# classes for managing FLUKA regions
# A. Mereghetti, 2024/04/16
# python version: >= 3.8.10

import numpy as np
from FLUKA import GeoObject, assembleLine, MaxLineLength, LineHeader, cleanRegLine, remSpaceAfters, addSpaceAfters

class Region(GeoObject):
    '''
    - no parsing/splitting of zones;
    - support only for preceding comments or commented lines between region
      definition lines, NO comments headed by !
    - a region definition always starts at column 1;
    - no check of length of lines for region definition;
    - no support for LATTICE cards at parsing;
    - ASSIGNMA cards: only 1:1 material:region assignment, NO material
      assignment for decay radiation simulation, NO magnetic/electric fields;
    '''

    def __init__(self,myName="",myComment=""):
        GeoObject.__init__(self,myName=myName,myComment=myComment)
        self.neigh=5
        self.definition=""
        self.material="BLACKHOLE"
        self.LatticeName=None
        self.TransformName=None
        # additional fields
        self.initCont()

    def initCont(self,rCont=0,rCent=np.array([0.0,0.0,0.0]),rMaxLen=0.0):
        self.rCont=rCont # containment indicator (1 per region):
                         # -1: region to be contained into/sized by another one
                         #  0: regular region (neither contains nor it is contained)
                         #  1: cell region (it contains another region)
        self.rCent=rCent # central point of one or more region (3D array)
        self.rMaxLen=rMaxLen # max length from hive cell

    @staticmethod
    def fromBuf(myBuffer):
        newReg=Region()
        lHeadParsed=False
        for tmpLine in myBuffer.splitlines():
            if (not lHeadParsed):
                if (tmpLine.startswith("*")):
                    newReg.tailMe(tmpLine)
                else:
                    data=tmpLine.split()
                    newReg.rename(data[0],lNotify=False)
                    newReg.neigh=float(data[1])
                    newReg.definition=tmpLine
                    # remove region name and number of neighbour regions from
                    #       definition
                    newReg.definition=newReg.definition.replace(data[0],"",1)
                    newReg.definition=newReg.definition.replace(data[1],"",1)
                    # remove heading/trailing empty spaces
                    newReg.definition=newReg.definition.strip() 
                    lHeadParsed=True
            else:
                newReg.definition=newReg.definition+"\n"+tmpLine
        return newReg
        
    def BodyNameReplaceInDef(self,oldNames,newNames):
        'query-replace body names in region definition'
        for oName,nName in zip(oldNames,newNames):
            self.definition=self.definition.replace(oName,nName)

    def assignMat(self,myMaterial):
        self.material=myMaterial

    def assignLat(self,myLatName=None):
        if (myLatName is None): myLatName=self.echoName()
        self.LatticeName=myLatName

    def assignTrasf(self,myTrasfName=None):
        if (myTrasfName is None): myTrasfName=self.echoName()
        self.TransformName=myTrasfName

    def addZone(self,myDef,nSpace=16):
        if (self.isNonEmpty()):
            if ( '|' not in self.definition ):
                # if more than a zone, be sure that the first zone has the "|" preceeding char,
                #    just for the sake of proper alignment
                self.definition="| %s"%(self.definition)
            self.definition=self.definition+'\n'+" "*nSpace+'| %s'%(myDef)
        else:
            self.definition=myDef

    def isNonEmpty(self):
        return len(self.definition)>0

    def isLattice(self):
        return self.LatticeName is not None
        
    def isLinkedToTransform(self):
        return self.TransformName is not None

    def echoLatticeName(self):
        return self.LatticeName
    
    def echoTransformName(self):
        return self.TransformName
    
    def echo(self,lMat=False,lLat=False,lFree=True,maxLen=MaxLineLength,header=LineHeader,lMultiLine=True,\
             remSpaceAfters=remSpaceAfters,addSpaceAfters=addSpaceAfters):
        '''take into account comment'''
        if (lMat):
            # echo ASSIGNMA card
            return GeoObject.echoComm(self)+"ASSIGNMA  %10s%10s" % ( self.material, self.echoName() )
        elif (lLat):
            # echo LATTICE card
            return GeoObject.echoComm(self)+"%-10s%10s%20s%10s%20s%-10s" % ( "LATTICE", self.echoName(), "", \
                                                        self.echoLatticeName(), "", self.echoTransformName() )
        else:
            myString=""
            myStrings=["%-8s"%(self.echoName()),"%4d"%(self.neigh)]
            lFirst=True
            for tmpString in self.definition.splitlines():
                if (tmpString.startswith("*")):
                    # print comment lines as they are
                    addLine=tmpString
                    if (lFirst):
                        myString=assembleLine(myStrings,maxLen=maxLen,header=header,lMultiLine=lMultiLine)
                else:
                    # print definition line respecting all FLUKA constraints
                    ttString=cleanRegLine(tmpString.strip(),remSpaceAfters=remSpaceAfters,addSpaceAfters=addSpaceAfters)
                    if (lFirst):
                        addLine=assembleLine(myStrings+ttString.split(),maxLen=maxLen,header=header,lMultiLine=lMultiLine)
                    else:
                        addLine=assembleLine([header]+ttString.split(),maxLen=maxLen,header=header,lMultiLine=lMultiLine)

                liasChar="\n"
                if (lFirst):
                    liasChar=""
                    lFirst=False

                myString=myString+liasChar+addLine
            return GeoObject.echoComm(self)+myString

    def retBodiesInDef(self):
        'return list of body names in region definition'
        bodiesInDef=[]
        if (self.isNonEmpty):
            for myLine in self.definition.splitlines():
                if (myLine.startswith("*")): continue # comment line
                tmpLine=myLine.replace("|"," ")
                tmpLine=tmpLine.replace("+"," ")
                tmpLine=tmpLine.replace("-"," ")
                currBodies=tmpLine.strip().split()
                for currBody in currBodies:
                    tmpBody=currBody.strip()
                    if (tmpBody not in bodiesInDef):
                        bodiesInDef.append(tmpBody)
        else:
            print("Region %s is empty!"%(self.echoName()))
            exit(1)
        return bodiesInDef

    def merge(self,newReg,spacing=" "*16):
        # warn user in comment
        self.tailMe("* --> merged with region %s <--"%(newReg.echoName()))
        # merge comments
        if (len(newReg.comment)>0):
            self.tailMe(newReg.comment)
        # merge definitions
        mergedDef=""
        newRegZones=newReg.definition.split("|")
        myRegZones=self.definition.split("|")
        if (len(myRegZones)==1 and len(myRegZones[0])==0):
            # empty definition of self
            mergedDef="%s"%(newReg.definition)
            tmpComment="* copying definition of %s into %s"%(\
                        newReg.echoName(),self.echoName())
            self.tailMe(tmpComment)
        else:
            lFirst=True
            for ii,myRegZone in enumerate(myRegZones,1):
                sMyRegZone=myRegZone.strip()
                if (len(sMyRegZone)>0):
                    for jj,newRegZone in enumerate(newRegZones,1):
                        sNewRegZone=newRegZone.strip()
                        if (len(sNewRegZone)>0):
                            tmpComment="* merging zone %s:%d into %s:%d"%(\
                                    newReg.echoName(),jj,self.echoName(),ii)
                            if (lFirst):
                                self.tailMe(tmpComment)
                                mergedDef="| %s %s"%(sMyRegZone,sNewRegZone)
                                lFirst=False
                            else:
                                mergedDef=mergedDef+"\n%s"%(tmpComment)+"\n%s| %s %s"%(spacing,sMyRegZone,sNewRegZone)
        self.definition=mergedDef
        # merge neighbours
        self.neigh=self.neigh+newReg.neigh
        # remove any sign of merge labelling
        self.initCont()

