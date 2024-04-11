# classes for managing FLUKA geometries
# A. Mereghetti, 2023/09/13
# python version: >= 3.8.10;

import numpy as np
from copy import deepcopy

import myMath
import grid

StringPrec=1.0E-14
MaxLineLength=132
LineHeader="%13s"%("")
cut0s="0"*4

class GeoObject():
    '''
    A very basic object, implementing simply a name and a comment.
    This object is supposed to collect all the stuff concerning
      handling of names and modification of comments of instances
      of the other classes;
    '''

    def __init__(self,myName="",myComment=""):
        self.name=myName
        self.comment=myComment

    def rename(self,newName,lNotify=True):
        'change instance name'
        if (lNotify):
            self.tailMe("* NAME CHANGE: FROM %s TO %s"%(self.name,newName))
        self.name=newName

    def headMe(self,myString):
        'head a string to the comment'
        if (len(self.comment)>0):
            self.comment=myString+"\n"+self.comment
        else:
            self.comment=myString
            
    def tailMe(self,myString):
        'tail a string to the comment'
        if (len(self.comment)>0):
            self.comment=self.comment+"\n"+myString
        else:
            self.comment=myString
                
    def echoComm(self):
        myBuf=""
        if (len(self.comment)>0):
            myBuf=self.comment+"\n"
        return myBuf

    def echoName(self):
        return self.name
            
class Body(GeoObject):
    '''
    - supported bodies: planes, spheres, and TRCs, nothing else for the time being;
    - body definition on ONE line only, and always starts at col 1;
    - support only for preceding comments, NO comments headed by "!" or after
      body definition;
    - no support for $start_* cards;
    '''

    def __init__(self,myName="",myComment=""):
        '''
        default body is a FLUKA PLA
        '''
        GeoObject.__init__(self,myName=myName,myComment=myComment)
        self.bType="PLA"
        self.P=np.zeros(3)
        self.V=np.array([0.0,0.0,1.0])
        self.Rs=np.zeros(2)
        self.lIsRotatable=True
        self.TransformName=None

    def echo(self,expDigits=None,numDigist=None,prec=StringPrec,cut0s=cut0s,maxLen=MaxLineLength,header=LineHeader,lMultiLine=True):
        '''take into account comment'''
        myStr="%-3s %8s"%(self.bType,self.echoName())
        if (self.bType=="PLA"):
            myFloats=list(self.V)+list(self.P)
        elif (self.bType=="YZP"): 
            myFloats=self.P[0]
        elif (self.bType=="XZP"):
            myFloats=self.P[1]
        elif (self.bType=="XYP"):
            myFloats=self.P[2]
        elif (self.bType=="SPH"):
            myFloats=list(self.P)+[self.Rs[0]]
        elif (self.bType=="TRC"):
            myFloats=list(self.P)+list(self.V)+list(self.Rs)
        elif (self.bType=="RCC"):
            myFloats=list(self.P)+list(self.V)+[self.Rs[0]]
        elif (self.bType=="XCC"):
            myFloats=list(self.P[1:])+[self.Rs[0]]
        elif (self.bType=="YCC"):
            myFloats=[self.P[2],self.P[0],self.Rs[0]]
        elif (self.bType=="ZCC"):
            myFloats=list(self.P[0:-1])+[self.Rs[0]]
        else:
            print("body %s NOT supported yet!"%(self.bType))
            exit(1)
        myStrings=[myStr]+echoFloats(myFloats,expDigits=expDigits,numDigist=numDigist,prec=prec,cut0s=cut0s)
        myStr=assembleLine(myStrings,maxLen=maxLen,header=header,lMultiLine=lMultiLine)
        if (self.TransformName is not None):
            myStr="$Start_transform -%s\n%s\n$end_transform"%(self.TransformName,myStr)
        return GeoObject.echoComm(self)+myStr

    def isRotatable(self):
        return self.lIsRotatable

    def linkTransformName(self,Tname=None):
        if (Tname is not None or Tname!=""):
            self.TransformName=Tname
            
    def retTransformName(self):
        return self.TransformName

    @staticmethod
    def fromBuf(tmpLines):
        newBody=Body()
        for tmpLine in tmpLines.splitlines():
            data=tmpLine.split()
            if (data[0]=="PLA"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.V=np.array(data[2:5]).astype(float)
                newBody.P=np.array(data[5:8]).astype(float)
            elif (data[0]=="YZP"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P[0]=float(data[2])
                newBody.lIsRotatable=False
            elif (data[0]=="XZP"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P[1]=float(data[2])
                newBody.lIsRotatable=False
            elif (data[0]=="XYP"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P[2]=float(data[2])
                newBody.lIsRotatable=False
            elif (data[0]=="SPH"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P=np.array(data[2:5]).astype(float)
                newBody.Rs[0]=float(data[5])
            elif (data[0]=="TRC"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P=np.array(data[2:5]).astype(float)
                newBody.V=np.array(data[5:8]).astype(float)
                newBody.Rs=np.array(data[8:]).astype(float)
            elif (data[0]=="RCC"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P=np.array(data[2:5]).astype(float)
                newBody.V=np.array(data[5:8]).astype(float)
                newBody.Rs[0]=float(data[8])
            elif (data[0]=="XCC"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P[1:]=np.array(data[2:4]).astype(float)
                newBody.Rs[0]=float(data[4])
                newBody.lIsRotatable=False
            elif (data[0]=="YCC"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P[0]=data[3]; newBody.P[2]=data[2]
                newBody.Rs[0]=float(data[4])
                newBody.lIsRotatable=False
            elif (data[0]=="ZCC"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P[0:-1]=np.array(data[2:4]).astype(float)
                newBody.Rs[0]=float(data[4])
                newBody.lIsRotatable=False
            elif (tmpLine.startswith("*")):
                newBody.tailMe(tmpLine)
            else:
                print("body %s NOT supported yet!"%(data[0]))
                exit(1)
        return newBody

    def makeRotatable(self,lDebug=False,infL=1000.0):
        if (not self.isRotatable()):
           myStr=GeoObject.echoComm(self)+"%-3s %8s"%(self.bType,self.echoName())
           if (self.bType=="YZP"):
               if (lDebug):
                   print("...converting body %s into a PLA"%(myStr))
               self.bType="PLA"
               self.V=np.array([1.0,0.0,0.0])
           elif (self.bType=="XZP"):
               if (lDebug):
                   print("...converting body %s into a PLA"%(myStr))
               self.bType="PLA"
               self.V=np.array([0.0,1.0,0.0])
           elif (self.bType=="XYP"):
               if (lDebug):
                   print("...converting body %s into a PLA"%(myStr))
               self.bType="PLA"
               self.V=np.array([0.0,0.0,1.0])
           elif (self.bType=="XCC"):
               if (lDebug):
                   print("...converting body %s into an RCC"%(myStr))
               self.bType="RCC"
               self.P[0]=-infL
               self.V=np.array([2*infL,0.0,0.0])
           elif (self.bType=="YCC"):
               if (lDebug):
                   print("...converting body %s into an RCC"%(myStr))
               self.bType="RCC"
               self.P[1]=-infL
               self.V=np.array([0.0,2*infL,0.0])
           elif (self.bType=="ZCC"):
               if (lDebug):
                   print("...converting body %s into an RCC"%(myStr))
               self.bType="RCC"
               self.P[2]=-infL
               self.V=np.array([0.0,0.0,2*infL])
           else:
               print("cannot make body %s rotatable!"%(myStr))
               exit(1)
           self.lIsRotatable=True

    def traslate(self,dd=None):
        if (dd is not None):
            self.P=self.P+dd
        else:
            print("please provide me with a 3D traslation array!")
            exit(1)
                
    def rotate(self,myMat=None,myTheta=None,myAxis=3,lDegs=True,lDebug=False):
        if (myMat is not None):
            self.P=myMat.mulArr(self.P,lDebug=lDebug)
            self.V=myMat.mulArr(self.V,lDebug=lDebug)
        elif (myTheta is not None):
            if (myAxis is None):
                print("you must specify a rotation axis together with an angle!")
                exit(1)
            else:
                myMat=myMath.RotMat(myAng=myTheta,myAxis=myAxis,lDegs=lDegs,lDebug=lDebug)
                self.rotate(myMat=myMat,myTheta=None,myAxis=myAxis,lDegs=lDegs,lDebug=lDebug)

    def resize(self,newL,lDebug=False):
        if (self.bType=="RCC"):
           if (lDebug):
               myStr=GeoObject.echoComm(self)+"%-3s %8s"%(self.bType,self.echoName())
               print("...resizing body %s to L=%g"%(myStr,newL))
           meanPoint=self.P+self.V/2.
           versor=self.V/np.linalg.norm(self.V)
           self.V=versor*newL
           self.P=meanPoint-self.V/2.
                
class Region(GeoObject):
    '''
    - no parsing/splitting of zones;
    - support only for preceding comments or commented lines between region
      definition lines, NO comments headed by !
    - a region definition always starts at column 1;
    - no check of length of lines for region definition;
    - no support for LATTICE;
    - ASSIGNMA cards: only 1:1 material:region assignment, NO material
      assignment for decay radiation simulation, NO magnetic/electric fields;
    '''

    def __init__(self,myName="",myComment=""):
        GeoObject.__init__(self,myName=myName,myComment=myComment)
        self.neigh=5
        self.definition=""
        self.material="BLACKHOLE"
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
        
    def echo(self,lMat=False):
        '''take into account comment'''
        if (lMat):
            # echo ASSIGNMA card
            return GeoObject.echoComm(self)+"ASSIGNMA  %10s%10s" % ( self.material, self.echoName() )
        else:
            return GeoObject.echoComm(self)+"%-8s   %4d %s" % \
                ( self.echoName(), self.neigh, self.definition )

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

class RotDefi(GeoObject):
    '''
    Class handling ROT-DEFI cards.
    Please have a look also to the Transformation() class: the ROT-DEFI
       class simply stores infos of a specific rotation/traslation;
    - echo in FREE format uses an empty space as field separator;
    - parsing in FREE format expects an empty space as field separator;
    '''

    def __init__(self,myAx=3,myPhi=0.0,myTh=0.0,myDD=np.zeros(3),myComment=""):
        # mainly for comments, name NOT used!
        GeoObject.__init__(self,myComment=myComment)
        self.axis=myAx           # FLUKA coding: 1=x, 2=y, 3=z
        self.phi=myPhi           # polar angle [degrees] (FLUKA units)
        self.theta=myTh          # azimuth angle [degrees] (FLUKA units)
        self.DD=myDD             # translation [cm] (x,y,z components)

    def echo(self,lFree=True,myID=0,myName="",\
             expDigits=None,numDigist=None,prec=StringPrec,cut0s=cut0s,\
             maxLen=MaxLineLength,header=LineHeader,lMultiLine=False):
        '''echo in FREE format uses an empty space as field separator'''
        if (myID<=0):
            print("...cannot dump a ROT-DEFI card without an index!")
            exit(1)
        if (len(myName)==0):
            print("...cannot dump a ROT-DEFI card without a name!")
            exit(1)
        if (myID<100):
            myIDecho=myID
            if (self.phi!=0.0 or self.theta!=0.0):
                myIDecho=myIDecho+self.axis*100
        else:
            myIDecho=myID*1000
            if (self.phi!=0.0 or self.theta!=0.0):
                myIDecho=myIDecho+self.axis
        if (lFree):
            myFloats=[self.phi,self.theta]+list(self.DD)
            myStrings=["ROT-DEFI  %10.1f"%(myIDecho)]+\
                echoFloats(myFloats,expDigits=expDigits,numDigist=numDigist,prec=prec,cut0s=cut0s)+\
                [" %-10s"%(myName)]
            return GeoObject.echoComm(self)+ \
                assembleLine(myStrings,maxLen=maxLen,header=header,lMultiLine=lMultiLine)
        else:
            return GeoObject.echoComm(self)+ \
                "ROT-DEFI  %10.1f% 10G% 10G% 10G% 10G% 10G%-10s"\
                %( myIDecho, self.phi, self.theta, self.DD[0], self.DD[1], self.DD[2], myName)

    @staticmethod
    def fromBuf(myBuffer,lFree=True):
        '''parsing in FREE format expects an empty space as field separator'''
        newRotDefi=RotDefi()
        for tmpLine in myBuffer.splitlines():
            if (tmpLine.startswith("*")):
                newRotDefi.tailMe(tmpLine.strip())
            else:
                # parse data
                if (lFree):
                    data=tmpLine.split()
                    pID=data[1]
                    phi=data[2]
                    tht=data[3]
                    DD=data[4:7]
                    myName=data[7]
                else:
                    pID=tmpLine[10:20].strip()
                    phi=tmpLine[20:30].strip()
                    tht=tmpLine[30:40].strip()
                    DD=[tmpLine[40:50].strip(),tmpLine[50:60].strip(),tmpLine[60:70].strip()]
                    myName=tmpLine[70:80].strip()
                # store data
                pID=float(pID)
                if (pID>1000):
                    myID=int(pID/1000)
                    myAxis=int(pID%1000)
                else:
                    myAxis=int(pID/100)
                    myID=int(pID%100)
                if (myAxis>0):
                    newRotDefi.axis=myAxis
                newRotDefi.phi=float(phi)
                newRotDefi.theta=float(tht)
                newRotDefi.DD=np.array([float(myDD) for myDD in DD])
        return newRotDefi, myID, myName

class Transformation(GeoObject):
    '''
    This is the actual class to be used when dealing with ROT-DEFI cards.
    The class implements an array of ROT-DEFI cards - hence, at least one card.
    '''

    def __init__(self,myID=0,myName=""):
        GeoObject.__init__(self,myName=myName) # mainly for names, comments NOT used!
        self.RotDefis=[] # list of ROT-DEFI cards
        self.ID=myID
        
    def __len__(self):
        return len(self.RotDefis)

    def __iter__(self):
        self.current_index=0
        return self
    
    def __next__(self):
        '''
        from https://blog.finxter.com/python-__iter__-magic-method/
        '''
        if self.current_index < len(self):
            x = self.RotDefis[self.current_index]
            self.current_index += 1
            return x
        raise StopIteration

    # override methods of GeoObject, since comments are of the single
    #   ROT-DEFI cards
    def headMe(self,myString):
        self.RotDefis[0].headMe(myString)
    def tailMe(self,myString):
        self.RotDefis[0].tailMe(myString)
    def echoComm(self):
        myComm=""
        if (len(self.RotDefis)>0):
            myComm=self.RotDefis[0].echoComm()
        return myComm

    def AddRotDefi(self,newRotDefi,iAdd=-1,lMoveComment=True):
        if (iAdd==-1 or iAdd==len(self.RotDefis)):
            # append ROT-DEFI card to list of ROT-DEFIs
            self.RotDefis.append(newRotDefi)
        elif (iAdd==0):
            # head ROT-DEFI card to list of ROT-DEFIs
            if (lMoveComment):
                myComm=self.echoComm().strip()
                if (len(myComm)>0):
                    # take care of main comment
                    newRotDefi.headMe(myComm)
                    self.RotDefis[0].comment=""
            self.RotDefis=[newRotDefi]+self.RotDefis
        else:
            # insert ROT-DEFI card at specific position in the list of ROT-DEFIs
            self.RotDefis=self.RotDefis[0:iAdd]+[newRotDefi]+self.RotDefis[iAdd:]
    def AddRotDefis(self,newRotDefis,iAdd=-1,lMoveComment=True):
        for myRotDefi in newRotDefis:
            self.AddRotDefi(myRotDefi,iAdd=iAdd,lMoveComment=lMoveComment)

    @staticmethod
    def fromBuf(self,myBuffer,lFree=True):
        newTransf=Transformation()
        myBuf=""
        for tmpLine in myBuffer.splitlines():
            if (tmpLine.startswith("*")):
                if (len(myBuf)==0):
                    myBuf=tmpLine
                else:
                    myBuf=myBuf+"\n"+tmpLine
            else:
                tmpRotDefi=RotDefi()
                myID,myName=tmpRotDefi.fromBuf(myBuf)
                if (myID!=0):
                    if (newTransf.ID==0):
                        # first ROT-DEFI card of transformation
                        newTransf.ID=myID
                    elif(newTransf.ID!=myID):
                        print("...cannot add ROT-DEFI with ID=%d to tranformation ID=%d!"%(myID,newTransf.ID))
                        exit(1)
                if (len(myName)>0):
                    if (len(newTransf.echoName())==0):
                        # first ROT-DEFI card of transformation
                        newTransf.rename(myName,lNotify=False)
                    elif (newTransf.echoName()!=myName):
                        print("...cannot add ROT-DEFI named %s to tranformation named %s !"%(myName,newTransf.echoName()))
                        exit(1)
                newTransf.AddRotDefi(tmpRotDefi)
                myBuf=""
        return newTransf

    def echo(self,lFree=True):
        buf=""
        for myRotDefi in self.RotDefis:
            if (len(buf)>0):
                buf=buf+"\n"
            buf=buf+myRotDefi.echo(myID=self.ID,myName=self.name)
        return buf

class Usrbin(GeoObject):
    '''
    A very basic class for handling USRBIN cards.
    For the time being, there is support only for:
    - parsing/echoing;
    - always 2 cards to fully define scoring;
    - a single ROTPRBIN card per transformation (1:1 mapping between ROTPRBIN and
      USRBIN cards); the mapping to the transformation is ALWAYS name-based;
      the ROTPRBIN card always PRECEEDs the respective USRBIN card;
    - NO comments between USRBIN cards defining the same scoring detector;
    - parsing ALWAYS in FIXED format;
    '''
    def __init__(self,myName="",myComment=""):
        GeoObject.__init__(self,myName=myName,myComment=myComment)
        self.definition=[] # array of strings storing the lines defining the
                           #   USRBIN. NB: USRBIN tags, & char and bin name
                           #   are NOT stored.
        self.TransfName=""

    def echo(self):
        '''take into account comment'''
        tmpBuf=GeoObject.echoComm(self)
        if (len(self.TransfName)>0):
            tmpBuf=tmpBuf+"%-10s%10s%10s%10s%10s\n"\
                %("ROTPRBIN","",self.TransfName,"",self.echoName())
        for myDef,mySdum,myEoL in zip(self.definition,[self.echoName(),"&"],["\n",""]):
            tmpBuf=tmpBuf+"%-10s%60s%-10s%-s"%("USRBIN",myDef,mySdum,myEoL)
        return tmpBuf
                        
    @staticmethod
    def fromBuf(myBuffer):
        newUsrBin=Usrbin()
        for myLine in myBuffer.split("\n"):
            myLine=myLine.strip()
            if (myLine.startswith("*")):
                newUsrBin.tailMe(myLine)
            elif (myLine.startswith("ROTPRBIN")):
                newUsrBin.assignTrasf(myLine[20:30].strip())
                if (len(newUsrBin.TransfName)==0):
                    print("...something wrong when parsing USRBIN cards: no ROT-DEFI card name!")
                    exit(1)
            else:
                newUsrBin.definition.append(myLine[10:70])
                if (len(newUsrBin.definition)==1):
                    # from first line defining USRBIN, get USRBIN name
                    newUsrBin.rename(myLine[70:].strip(),lNotify=False)
        if (len(newUsrBin.definition)!=2):
            print("...something wrong when parsing USRBIN cards: got %d lines!"\
                  %(len(newUsrBin.definition)))
            exit(1)
        if (len(newUsrBin.name)==0):
            print("...something wrong when parsing USRBIN cards: no name!")
            exit(1)
        return newUsrBin

    def getUnit(self):
        return int(abs(float(self.definition[0][20:30]))+1E-4)

    def setUnit(self,myUnit):
        self.definition[0]=self.definition[0][0:20]+"% 10.1f"%(myUnit)+self.definition[0][30:]

    def isSpecialBinning(self):
        binType=int(abs(float(self.definition[0][0:10]))+1E-4)
        return binType==8.0 or binType==18.0

    def getNbins(self):
        nBins=None
        if (not self.isSpecialBinning() ):
            nBins=int(abs(float(self.definition[1][30:40]))+1E-4)*\
                  int(abs(float(self.definition[1][40:50]))+1E-4)*\
                  int(abs(float(self.definition[1][50:  ]))+1E-4)
        return nBins

    def assignTrasf(self,trasfName):
        if (len(trasfName)>0):
            self.TransfName=trasfName
    def returnTrasf(self):
        return self.TransfName
    
class Geometry():
    '''
    - name-based FLUKA geometry defition;
    - NO support of LATTICE cards;
    - NO support for #include cards or geo defitions in files other than that
       being parsed;
    - comments:
      . in general, only comments preceeding a specific card are retained;
        --> trailing comments to body declaration section will disappear!
      . body: commented lines are considered always before the body;
      . region: commented lines are kept where they are found, if they are
                found before or along the region declaration;
      . rot-defi: commented lines are considered always before the rot-defi card;
      . usrbin: commented lines are considered always before the URSBIN card;
    '''
    def __init__(self):
        self.bods=[]
        self.regs=[]
        self.tras=[]
        self.bins=[]
        self.title=""

    def add(self,tmpEl,what="body"):
        '''
        tmpEl is already the instance to be added, NOT the string buffer to be parsed
        '''
        if (what.upper().startswith("BOD")):
            self.bods.append(tmpEl)
        elif (what.upper().startswith("REG")):
            self.regs.append(tmpEl)
        elif (what.upper().startswith("TRAS") or what.upper().startswith("TRANSF")):
            self.tras.append(tmpEl)
        elif (what.upper().startswith("BIN") or what.upper().startswith("USRBIN")):
            self.bins.append(tmpEl)
        else:
            print("...what do you want to add to geometry? %s NOT recognised!"%(what))
            exit(1)

    def setTitle(self,tmpTitle="My Geometry"):
        self.title=tmpTitle
        
    def assignma(self,tmpLine):
        'crunch info by ASSIGNMA card'
        data=tmpLine.split()
        myMatName=data[1]
        myRegName=data[2]
        print("...assigning material %s to region %s..."%(myMatName,myRegName))
        myEntry,iEntry=self.ret("REG",myRegName)
        self.regs[iEntry].assignMat(myMatName)

    def rotdefi(self,tmpBuf,lFree=True):
        'crunch info by ROT-DEFI card'
        tmpRotDefi,myID,myName=RotDefi.fromBuf(tmpBuf,lFree=lFree)
        myEntry,iEntry=self.ret("TRANSF",myName)
        if (iEntry==-1):
            # brand new transformation
            myTrans=Transformation(myID=myID,myName=myName)
            self.add(myTrans,"TRANSF")
            iEntry=0
        self.tras[iEntry].AddRotDefi(tmpRotDefi)

    def headMe(self,myString,begLine="* \n"+"* "+"="*130,endLine="* "+"-"*130+"\n* "):
        '''
        simple method to head a string to the geometry declaration (bodies,
           regions, assignma, usrbin, rot-defi cards)
        '''
        actualString=begLine+"\n* "+myString+" \n"+endLine
        if (len(self.bods)>0):
            self.bods[0].headMe(actualString)
        if (len(self.regs)>0):
            self.regs[0].headMe(actualString)
        if (len(self.tras)>0):
            self.tras[0].headMe(actualString)
        if (len(self.bins)>0):
            self.bins[0].headMe(actualString)
            
    def ret(self,what,myName):
        lFound=False
        if (what.upper().startswith("BODSINREG")):
            if (myName.upper()=="ALL"):
                myReg,iReg=self.ret("reg","ALL")
                myEntry=[]; iEntry=[]
                for findReg in myReg:
                    tmpEl,tmpInd=self.ret("BodsInReg",findReg)
                    myEntry=myEntry+tmpEl
                    iEntry=iEntry+tmpInd
            else:
                myReg,iReg=self.ret("reg",myName)
                myEntry=myReg.retBodiesInDef()
                iEntry=[]
                for findBod in myEntry:
                    tmpEl,tmpInd=self.ret("bod",findBod)
                    iEntry.append(tmpInd)
                lFound=True
        elif (what.upper().startswith("BOD")):
            if (myName.upper()=="ALL"):
                myEntry=[ body.echoName() for body in self.bods ]
                iEntry=[ ii for ii in range(len(self.bods)) ]
                lFound=True
            else:
                for iEntry,myEntry in enumerate(self.bods):
                    if (myEntry.echoName()==myName):
                        lFound=True
                        break
        elif (what.upper().startswith("REG")):
            if (myName.upper()=="ALL"):
                myEntry=[ reg.echoName() for reg in self.regs ]
                iEntry=[ ii for ii in range(len(self.regs)) ]
                lFound=True
            else:
                for iEntry,myEntry in enumerate(self.regs):
                    if (myEntry.echoName()==myName):
                        lFound=True
                        break
        elif (what.upper().startswith("TRANSF")):
            if (myName.upper()=="ALL"):
                myEntry=[ tras.echoName() for tras in self.tras ]
                iEntry=[ ii for ii in range(len(self.tras)) ]
                lFound=True
            else:
                for iEntry,myEntry in enumerate(self.tras):
                    if (myEntry.echoName()==myName):
                        lFound=True
                        break
        elif (what.upper().startswith("BIN") or what.upper().startswith("USRBIN")):
            if (myName.upper()=="ALL"):
                myEntry=[ myBin.echoName() for myBin in self.bins ]
                iEntry=[ ii for ii in range(len(self.bins)) ]
                lFound=True
            else:
                for iEntry,myEntry in enumerate(self.bins):
                    if (myEntry.echoName()==myName):
                        lFound=True
                        break
        else:
            print("%s not recognised! What should I look for in the geometry?"%(what))
            exit(1)
        if (not lFound):
            myEntry=None
            iEntry=-1
        return myEntry, iEntry

    @staticmethod
    def fromInp(myInpName):
        '''
        FREE format used only for parsing ROT-DEFI cards (blank space as separator)!
        '''
        newGeom=Geometry()
        print("parsing file %s..."%(myInpName))
        ff=open(myInpName,'r')
        lRead=0
        lFree=False
        tmpBuf=""
        regBuf=""
        for tmpLine in ff.readlines():
            
            # OUTSIDE GEOBEGIN-GEOEND
            if (lRead==0):
                if (tmpLine.startswith("GEOBEGIN")):
                    lRead=1
                elif (tmpLine.startswith("ASSIGNMA")):
                    newGeom.assignma(tmpLine)
                    tmpBuf="" # flush buffer
                elif (tmpLine.startswith("ROT-DEFI")):
                    tmpBuf=tmpBuf+tmpLine
                    newGeom.rotdefi(tmpBuf,lFree=lFree)
                    tmpBuf="" # flush buffer
                elif (tmpLine.startswith("ROTPRBIN")):
                    tmpBuf=tmpBuf+tmpLine
                elif (tmpLine.startswith("USRBIN")):
                    tmpBuf=tmpBuf+tmpLine
                    if(tmpLine.strip().endswith("&")):
                       newGeom.add(Usrbin.fromBuf(tmpBuf.strip()),what="BIN")
                       tmpBuf="" # flush buffer
                elif (tmpLine.startswith("*")):
                    tmpBuf=tmpBuf+tmpLine
                elif (tmpLine.startswith("FREE")):
                    lFree=True
                elif (tmpLine.startswith("FIXED")):
                    lFree=False
                else:
                    # another card, to be skipped
                    tmpBuf="" # flush buffer
                    
            # INSIDE GEOBEGIN-GEOEND
            elif (lRead==1):
                # title after GEOBEGIN
                newGeom.title=tmpLine[20:].strip()
                lRead=2
                tmpBuf="" # flush buffer
            elif (lRead==2):
                # definition of FLUKA bodies
                if (tmpLine.startswith("END")):
                    print("...acquired %d bodies;"%(len(newGeom.bods)))
                    lRead=3
                    tmpBuf="" # flush buffer
                elif (tmpLine.startswith("$start")):
                    print("$start* cards NOT supported!")
                    exit(1)
                else:
                    tmpBuf=tmpBuf+tmpLine
                    if (not tmpLine.startswith("*")):
                        # acquire body
                        newGeom.add(Body.fromBuf(tmpBuf.strip()),what="BODY")
                        tmpBuf="" # flush buffer
            elif (lRead==3):
                # definition of FLUKA regions
                if (tmpLine.startswith("END")):
                    if (len(regBuf)>0):
                        # acquire region
                        newGeom.add(Region.fromBuf(regBuf),what="reg")
                        regBuf="" # flush region def buffer
                    print("...acquired %d regions;"%(len(newGeom.regs)))
                    tmpBuf="" # flush buffer
                    lRead=0 # ignore LATTICE cards
                else:
                    if (tmpLine.startswith("*")):
                        # comment line: store in buffer
                        tmpBuf=tmpBuf+tmpLine
                        continue
                    if (tmpLine.startswith(" ")):
                        # region definition continues
                        regBuf=regBuf+tmpBuf+tmpLine
                        tmpBuf="" # flush buffer
                    else:
                        # a new region
                        if (len(regBuf)>0):
                            # acquire region previously read (if any)
                            newGeom.add(Region.fromBuf(regBuf),what="reg")
                            regBuf="" # flush region def buffer
                        regBuf=tmpBuf+tmpLine
                        tmpBuf="" # flush buffer
                    
        ff.close()
        print("...acquired %d USRBIN(s);"%(len(newGeom.bins)))
        print("...acquired %d transformation(s);"%(len(newGeom.tras)))
        print("...done;")
        return newGeom

    @staticmethod
    def appendGeometries(myGeos,myTitle=None):
        '''
        Barely appending geometries one to another;
        '''
        new=Geometry()
        for myGeo in myGeos:
            new.bods=new.bods+myGeo.bods
            new.regs=new.regs+myGeo.regs
            for ii,myTras in enumerate(myGeo.tras,1):
                myTras.ID=ii+len(new.tras)
            new.tras=new.tras+myGeo.tras
            new.bins=new.bins+myGeo.bins
        if (myTitle is None):
            myTitle="appended geometries"
        new.title=myTitle
        return new

    def echo(self,oFileName,lSplit=False,what="all",dMode="w"):
        '''
        - what="all"/"bodies"/"regions"/"materials"/"transformation"/"bins"
        In case oFileName ends with ".geo", lSplit is activated,
           overriding the user request, and bodies/regions are
           dumped in the .geo file, whereas assignmat cards are
           dumped in the _assignmat.inp file.
        '''
        import os
        
        if (not oFileName.endswith(".inp") and not oFileName.endswith(".geo")):
            oFileName=oFileName+".inp"
        if (oFileName.endswith(".geo")):
            lSplit=True

        if (what.upper().startswith("BOD")):
            print("saving bodies in file %s..."%(oFileName))
            ff=open(oFileName,dMode)
            ff.write("* \n")
            for tmpBody in self.bods:
                ff.write("%s\n"%(tmpBody.echo()))
            ff.write("* \n")
            ff.close()
            print("...done;")
        elif (what.upper().startswith("REG")):
            print("saving regions in file %s..."%(oFileName))
            ff=open(oFileName,dMode)
            ff.write("* \n")
            for tmpReg in self.regs:
                ff.write("%s\n"%(tmpReg.echo()))
            ff.write("* \n")
            ff.close()
            print("...done;")
        elif (what.upper().startswith("TRANSF")):
            print("saving ROT-DEFI cards in file %s..."%(oFileName))
            ff=open(oFileName,dMode)
            ff.write("* \n")
            ff.write("FREE \n")
            for tmpTrasf in self.tras:
                ff.write("%s\n"%(tmpTrasf.echo())) # FREE format by default
            ff.write("FIXED \n")
            ff.write("* \n")
            ff.close()
            print("...done;")
        elif (what.upper().startswith("MAT")):
            print("saving ASSIGNMA cards in file %s..."%(oFileName))
            ff=open(oFileName,dMode)
            ff.write("* \n")
            for tmpReg in self.regs:
                ff.write("%s\n"%(tmpReg.echo(lMat=True)))
            ff.write("* \n")
            ff.close()
            print("...done;")
        elif (what.upper().startswith("BIN")):
            print("saving USRBINs in file %s..."%(oFileName))
            ff=open(oFileName,dMode)
            ff.write("* \n")
            for tmpBin in self.bins:
                ff.write("%s\n"%(tmpBin.echo()))
            ff.write("* \n")
            ff.close()
            print("...done;")
        elif (what.upper()=="ALL"):
            if (lSplit):
                if (oFileName.endswith(".inp")):
                    # split geometry definition into bodies,
                    #   regions, assignmat, rotdefi and usrbin files,
                    #   to be used with pure #include statements;
                    self.echo(oFileName.replace(".inp","_bodies.inp",1),\
                              lSplit=False,what="bodies",dMode="w")
                    self.echo(oFileName.replace(".inp","_regions.inp",1),\
                              lSplit=False,what="regions",dMode="w")
                    self.echo(oFileName.replace(".inp","_assignmats.inp",1),\
                              lSplit=False,what="materials",dMode="w")
                    if (len(self.tras)>0):
                        self.echo(oFileName.replace(".inp","_rotdefis.inp",1),\
                              lSplit=False,what="transf",dMode="w")
                    if (len(self.bins)>0):
                        self.echo(oFileName.replace(".inp","_usrbins.inp",1),\
                              lSplit=False,what="bin",dMode="w")
                else:
                    # split geometry definition into a .geo file
                    #   and an assignmat, rotdefi and usrbin files;
                    #   the former is encapsulated between GEOBEGIN and GEOEND cards,
                    #   the other is imported via an #include statement
                    ff=open(oFileName,"w")
                    ff.write("% 5d% 5d%10s%s\n"%(0,0,"",self.title))
                    ff.close()
                    self.echo(oFileName,lSplit=False,what="bodies",dMode="a")
                    ff=open(oFileName,"a")
                    ff.write("%-10s\n"%("END"))
                    ff.close()
                    self.echo(oFileName,lSplit=False,what="regions",dMode="a")
                    ff=open(oFileName,"a")
                    ff.write("%-10s\n"%("END"))
                    ff.close()
                    self.echo(oFileName.replace(".geo","_assignmats.inp",1),\
                              lSplit=False,what="materials",dMode="w")
                    if (len(self.tras)>0):
                        self.echo(oFileName.replace(".geo","_rotdefis.inp",1),\
                              lSplit=False,what="transf",dMode="w")
                    if (len(self.bins)>0):
                        self.echo(oFileName.replace(".geo","_usrbins.inp",1),\
                              lSplit=False,what="bin",dMode="w")
            else:
                ff=open(oFileName,"w")
                ff.write("%-10s%60s%-10s\n"%("GEOBEGIN","","COMBNAME"))
                ff.write("% 5d% 5d%10s%s\n"%(0,0,"",self.title))
                ff.close()
                self.echo(oFileName,lSplit=False,what="bodies",dMode="a")
                ff=open(oFileName,"a")
                ff.write("%-10s\n"%("END"))
                ff.close()
                self.echo(oFileName,lSplit=False,what="regions",dMode="a")
                ff=open(oFileName,"a")
                ff.write("%-10s\n"%("END"))
                ff.write("%-10s\n"%("GEOEND"))
                ff.close()
                if (len(self.tras)>0):
                    self.echo(oFileName,lSplit=False,what="transf",dMode="a")
                self.echo(oFileName,lSplit=False,what="materials",dMode="a")
                if (len(self.bins)>0):
                    self.echo(oFileName,lSplit=False,what="bin",dMode="a")
        else:
            print("...what should I echo? %s NOT reconised!"%(what))

    def solidTrasform(self,dd=None,myMat=None,myTheta=None,myAxis=3,lDegs=True,lGeoDirs=False,lDebug=False):
        '''
        Bodies are applied the transformation as requested by the user:
        * first, rotation (meant to be expressed in ref sys of the object to rotate);
        * second, translation.

        ROT-DEFI cards express the opposite transformation (it is used runtime
          in particle tracking, so it is used to track particles in the original
          position). In addition, in ROT-DEFI cards:
        * the translation is applied first, then the rotation;
        * the rotation angles (azimuth) have opposite sign wrt a rotation angle
          in a right-handed system.
        '''
        print("applying solid transformation(s)...")
        if (myMat is None and myTheta is None and dd is None):
            print("...no transformation provided!")
        else:
            ROTDEFIlist=[]
            myName="ToFinPos"
            if (myMat is not None):
                print("...applying transformation expressed by matrix to geometry...")
                if (not lGeoDirs):
                    for ii in range(len(self.bods)):
                        self.bods[ii].rotate(myMat=myMat,myTheta=None,myAxis=None,\
                                             lDegs=lDegs,lDebug=lDebug)
                thetas=myMat.GetGimbalAngles() # [degs]
                for myAx,myTh in enumerate(thetas,1):
                    if (myTh!=0.0):
                        ROTDEFIlist.append(RotDefi(myAx=myAx,myTh=myTh))
            elif (myTheta is not None):
                if (type(myTheta) is list and type(myAxis) is list):
                    if ( len(myTheta)!=len(myAxis) ):
                        print("...inconsistent number of angles (%d) and axes (%d)!"\
                              %(len(myTheta),len(myAxis)))
                        exit(1)
                    for myT,myAx in zip(myTheta,myAxis):
                        # iteratively call this same method, on a single angle
                        self.solidTrasform(myTheta=myT,myAxis=myAx,lDegs=lDegs,lGeoDirs=lGeoDirs,lDebug=lDebug)
                elif (myTheta!=0.0):
                    print("...applying rotation by %f degs around axis %d..."%\
                          (myTheta,myAxis))
                    if (not lGeoDirs):
                        for ii in range(len(self.bods)):
                            self.bods[ii].rotate(myMat=None,myTheta=myTheta,myAxis=myAxis,\
                                                 lDegs=lDegs,lDebug=lDebug)
                    ROTDEFIlist.append(RotDefi(myAx=myAxis,myTh=myTheta))
            if (dd is not None):
                print("...applying traslation array [%f,%f,%f] cm..."%\
                      (dd[0],dd[1],dd[2]))
                if (not lGeoDirs):
                    for ii in range(len(self.bods)):
                        self.bods[ii].traslate(dd=dd)
                myDD=-dd
                if (isinstance(myDD,list)):
                    myDD=np.array(DD)
                ROTDEFIlist.append(RotDefi(myDD=myDD))
            if (len(ROTDEFIlist)>0):
                # add transformation to actual position...
                lCreate=False
                if ( lGeoDirs and not lCreate ):
                    # ...if $start_tranform directives are used
                    myEntry, iEntry = self.ret("TRANSF",myName)
                    lCreate=myEntry is None
                if ( not lCreate ):
                    # ...if USRBINs are present without ROTPRBIN cards
                    for myBin in self.bins:
                        if ( myBin.returnTrasf()=="" ):
                            lCreate=True
                            myBin.assignTrasf(myName)
                if (lCreate):
                    print("...adding the final transformation to (existing) list of transformations...")
                    self.add(Transformation(myID=len(self.tras),myName=myName),what="tras")
                # - add list of ROT-DEFIs for the solid transformation to the list
                #   of existing transformation
                for myTras in self.tras:
                    myTras.AddRotDefis(reversed(ROTDEFIlist),iAdd=0)
                # - link solid trasformation
                if (lGeoDirs):
                    for ii in range(len(self.bods)):
                        self.bods[ii].linkTransformName(Tname=myName)
        print("...done.")

    def rename(self,newName,lNotify=True):
        print("renaming geometry...")
        maxLenName=8
        if (len(newName)>=maxLenName):
            print("Geometry.rename(): cannot rename entities with len(%s)>=%d!"%(newName,maxLenName))
            exit(1)
        newNameFmt=newName+"%0"+"%d"%(maxLenName-len(newName))+"d"
        oldBodyNames=[]; newBodyNames=[]
        oldTrasNames=[]; newTrasNames=[]
        for iBody in range(len(self.bods)):
            oldBodyNames.append(self.bods[iBody].echoName())
            newBodyNames.append(newNameFmt%(iBody+1))
            self.bods[iBody].rename(newBodyNames[-1],lNotify=lNotify)
        for iReg in range(len(self.regs)):
            self.regs[iReg].rename(newNameFmt%(iReg+1),lNotify=lNotify)
            self.regs[iReg].BodyNameReplaceInDef(oldBodyNames,newBodyNames)
        for iTras in range(len(self.tras)):
            oldTrasNames.append(self.tras[iTras].echoName())
            newTrasNames.append(newNameFmt%(iTras+1))
            self.tras[iTras].rename(newTrasNames[-1],lNotify=lNotify)
        for iBin in range(len(self.bins)):
            self.bins[iBin].rename(newNameFmt%(iBin+1),lNotify=lNotify)
            trName=self.bins[iBin].returnTrasf()
            if (len(trName)>0):
                lFound=False
                for oldTrasName,newTrasName in zip(oldTrasNames,newTrasNames):
                    if (trName==oldTrasName):
                        self.bins[iBin].assignTrasf(newTrasName)
                        lFound=True
                        break
                if (not lFound):
                    print("cannot find transformation named %s!"%(trName))
                    exit(1)
            else:
                print("Geometry.rename(): USRBIN with no original name in geometry!")
                exit(1)
        for iBody in range(len(self.bods)):
            if (self.bods[iBody].retTransformName() is not None):
                if (self.bods[iBody].retTransformName() not in oldTrasNames):
                    print("cannot find name of transformation for moving geo: %s!"%(self.bods[iBody].retTransformName()))
                    exit(1)
                ii=oldTrasNames.index(self.bods[iBody].retTransformName())
                self.bods[iBody].linkTransformName(Tname=newTrasNames[ii])
        print("...done.")

    def flagRegs(self,whichRegs,rCont,rCent):
        if (whichRegs is str):
            if (whichRegs.upper()=="ALL"):
                regs2mod,iRegs2mod=self.ret("reg","ALL")
            else:
                regs2mod=[whichRegs]
        else:
            if ( "ALL" in [tmpStr.upper() for tmpStr in whichRegs] ):
                regs2mod,iRegs2mod=self.ret("reg","ALL")
            else:
                regs2mod=whichRegs
        for whichReg in regs2mod:
            outReg,iReg=self.ret("reg",whichReg)
            outReg.initCont(rCont=rCont,rCent=rCent)

    def resizeBodies(self,newL,whichBods):
        '''
        input:
        - newL: new length [cm];
        - whichBods: list of body names to be updated (if infinite)
        '''
        if (whichBods is str):
            if (whichBods.upper()=="ALL"):
                bods2mod,iBods2mod=self.ret("bod","ALL")
            else:
                bods2mod=[whichBods]
        else:
            if ( "ALL" in [tmpStr.upper() for tmpStr in whichBods] ):
                bods2mod,iBods2mod=self.ret("bod","ALL")
            else:
                bods2mod=whichBods
        for bod2mod in bods2mod:
            whichBod,iBod=self.ret("bod",bod2mod)
            whichBod.resize(newL)

    def makeBodiesRotatable(self,lDebug=False,infL=1000.):
        for myBod in self.bods:
            myBod.makeRotatable(lDebug=lDebug,infL=infL)

    def reAssiginUSRBINunits(self,nMaxBins=None,nUSRBINs=None,usedUnits=None):
        print("re-assigning USRBIN units...")
        if (nMaxBins is not None and nUSRBINs is not None):
            print("Please tell me if I have to re-assign USRBIN units based on:")
            print("- max number of bins in a unit;")
            print("- max number of USRBINs in a unit;")
            exit(1)
        if (usedUnits is None): usedUnits=[]
        if (type(usedUnits) is not list): usedUnits=[usedUnits]
        uniqueUnits=list(set([ myBin.getUnit() for myBin in self.bins ]))
        print("...%d original units:"%(len(uniqueUnits)),uniqueUnits)
        if (len(usedUnits)>0):
            print("...units NOT to be used:",usedUnits)
        currUnits=deepcopy(uniqueUnits)
        for ii in range(len(currUnits)):
            while (currUnits[ii] in usedUnits):
                currUnits[ii]=currUnits[ii]+1
                if (currUnits[ii]>99):
                    print("...exceeding max number of supported units!")
                    exit(1)
        myN=[ 0 for ii in range(len(currUnits)) ]
        if (nMaxBins is not None):
            nMax=nMaxBins
        elif(nUSRBINs is not None):
            nMax=nUSRBINs
        else:
            print("Please tell me if I have to re-assign USRBIN units based on:")
            print("- max number of bins in a unit;")
            print("- max number of USRBINs in a unit;")
            exit(1)
        for myBin in self.bins:
            iUnit=uniqueUnits.index(myBin.getUnit())
            if (nMaxBins is not None):
                nAdd=myBin.getNbins()
            elif(nUSRBINs is not None):
                nAdd=1
            if (myN[iUnit]+nAdd>nMax):
                myUnit=currUnits[iUnit]
                while(myUnit in currUnits or \
                      myUnit in usedUnits ):
                    myUnit=myUnit+1
                    if (myUnit>99):
                        print("...exceeding max number of supported units!")
                        exit(1)
                currUnits[iUnit]=myUnit
                myN[iUnit]=nAdd
            else:
                myN[iUnit]=myN[iUnit]+nAdd
            myBin.setUnit(currUnits[iUnit])
        print("...done;")

    @staticmethod
    def DefineHive_SphericalShell(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,defMat="VACUUM",tmpTitle="Hive for a spherical shell",lWrapBHaround=False,lDebug=True):
        '''
        This method defines the hive for a grid on a spherical shell.

        The grid is described in spherical coordinates:
        - r [cm];
        - theta [degs]: angle in yz-plane (positive when pointing towards y>0);
        - phi [degs]: angle in xz-plane (positive when pointing towards x>0).
        The grid is centred around the z-axis.

        NR, NT, NP are number of points (cell centers) along r, theta and phi
          in the grid; therefore, the number of steps in the grid are NR-1,
          NT-1, NP-1. Similarly, the number of bodies that limit the grid
          along r, theta and phy are NR+1, NT+1 and NP+1

        The hive is delimited:
        - radially, by spheres;
        - on theta, by TRCs (i.e. identifying circles of latitude on the spherical shell);
        - on phi, by rotated YZPs (i.e. identifying meridians on the spherical shell);
        '''
        
        print("Preparing the hive for a spherical shell...")
        
        print("...generating the grid of cells...")
        cellGrid=grid.SphericalShell(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,lDebug=lDebug)
        print("...defining hive boundaries...")
        myHive=grid.SphericalHive(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,lDebug=lDebug)
        RRs,TTs,PPs=myHive.ret(what="all")

        print("...generating FLUKA geometry...")
        newGeom=Geometry()
        
        print("   ...bodies...")
        # - concentric spherical shells
        spheres=[]
        for ii,RR in enumerate(RRs,1):
            tmpBD=Body()
            tmpBD.rename("HVRAD%03i"%(ii),lNotify=False)
            tmpBD.bType="SPH"
            tmpBD.Rs[0]=RR
            tmpBD.comment="* Hive radial boundary at R[cm]=%g"%(RR)
            spheres.append(tmpBD)
        spheres[0].comment="* \n"+spheres[0].comment
        # - circles of latitude on the spherical shell
        thetas=[]
        for ii,TT in enumerate(TTs,1):
            tmpBD=Body()
            tmpBD.rename("HVTHT%03i"%(ii),lNotify=False)
            if (TT==0.0):
                tmpBD.V=np.array([0.0,1.0,0.0])
            else:
                tmpBD.bType="TRC"
                hh=RRs[-1]+10
                if ( myHive.SphericalHive_ThetaCoversPi() and ( ii==1 or ii==len(TTs) ) ):
                    # extremely small angle, almost 90degs,
                    #  almost negligible, but necessary for a more
                    #  robust geometry
                    tmpBD.Rs[0]=hh*1E-4
                else:
                    tmpBD.Rs[0]=hh*np.tan(np.radians(90-abs(TT)))
                if (TT>0):
                    hh=-hh # TT>0: angle towards y>0 --> TRC points downwards
                tmpBD.V=np.array([0.0, hh,0.0])
                tmpBD.P=np.array([0.0,-hh,0.0])
            tmpBD.comment="* Hive theta boundary at theta[deg]=%g"%(TT)
            thetas.append(tmpBD)
        thetas[0].comment="* \n"+thetas[0].comment
        # - meridians on the spherical shell
        phis=[]
        for ii,PP in enumerate(PPs,1):
            tmpBD=Body()
            tmpBD.rename("HVPHI%03i"%(ii),lNotify=False)
            tmpBD.V=np.array([1.0,0.0,0.0])
            tmpBD.rotate(myMat=None,myTheta=PP,myAxis=2,lDegs=True,lDebug=lDebug)
            tmpBD.comment="* Hive phi boundary at phi[deg]=%g"%(PP)
            phis.append(tmpBD)
        phis[0].comment="* \n"+phis[0].comment
        newGeom.bods=spheres+thetas+phis

        print("   ...regions...")
        # - outside hive
        tmpReg=Region()
        tmpReg.rename("HV_OUTER",lNotify=False)
        tmpReg.material=defMat
        tmpReg.addZone('-%-8s'%(spheres[-1].name))
        tmpReg.comment="* region outside hive"
        tmpReg.initCont(rCont=-1)
        newGeom.add(tmpReg,what="reg")
        # - inside hive
        tmpReg=Region()
        tmpReg.rename("HV_INNER",lNotify=False)
        tmpReg.material=defMat
        tmpReg.addZone('+%-8s'%(spheres[0].name))
        tmpReg.comment="* region inside hive"
        newGeom.add(tmpReg,what="reg")
        # - around hive
        tmpReg=Region()
        tmpReg.rename("HVAROUND",lNotify=False)
        tmpReg.material=defMat
        if (not cellGrid.HasPole("N")):
            tmpReg.addZone('+%-8s -%-8s +%-8s'%( \
                spheres[-1].echoName(),spheres[0].echoName(), thetas[-1].echoName() ))
        if (not cellGrid.HasPole("S")):
            tmpReg.addZone('+%-8s -%-8s +%-8s'%( \
                spheres[-1].echoName(),spheres[0].echoName(), thetas[ 0].echoName() ))
        if (not myHive.SphericalHive_PhiCovers2pi()):
            tmpReg.addZone('+%-8s -%-8s -%-8s -%-8s -%-8s +%-8s'%( \
                spheres[-1].echoName(),spheres[0].echoName(), thetas[-1].echoName(), thetas[ 0].echoName(), phis[-1].echoName(), phis[ 0].echoName() ))
        if (PPs[-1]-PPs[0]<180.0):
            tmpReg.addZone('+%-8s -%-8s -%-8s -%-8s -%-8s -%-8s'%( \
                spheres[-1].echoName(),spheres[0].echoName(), thetas[-1].echoName(), thetas[ 0].echoName(), phis[-1].echoName(), phis[ 0].echoName() ))
            tmpReg.addZone('+%-8s -%-8s -%-8s -%-8s +%-8s +%-8s'%( \
                spheres[-1].echoName(),spheres[0].echoName(), thetas[-1].echoName(), thetas[ 0].echoName(), phis[-1].echoName(), phis[ 0].echoName() ))
        if (tmpReg.isNonEmpty()):
            tmpReg.comment="* region around hive"
            newGeom.add(tmpReg,what="reg")
        # - actual hive
        iPs=[ iP for iP in range(1,len(phis)) ]
        if myHive.SphericalHive_PhiCovers2pi():
            # re-use first plane
            iPs=iPs+[0]
        iHive=0
        for iR in range(1,len(spheres)):
            if (cellGrid.HasPole("S")):
                tmpReg=Region()
                tmpReg.rename("HVCL%04i"%(iHive),lNotify=False)
                tmpReg.material=defMat
                tmpReg.addZone('+%-8s -%-8s +%-8s'%(spheres[iR].echoName(),spheres[iR-1].echoName(),thetas[0].echoName()))
                myCenter=cellGrid.ret(what="POINT",iEl=iHive)
                rMaxLen=1.1*max(RRs[iR]-RRs[iR-1],RRs[iR]*2*np.absolute(np.radians(TTs[0]-90)))
                tmpReg.tailMe("* - hive region %4d: SOUTH POLE! R[cm]=[%g:%g], theta[deg]=[-90:%g]"%(
                    iHive,RRs[iR-1],RRs[iR],TTs[0]))
                tmpReg.tailMe("*   center=[%g,%g,%g];"%(myCenter[0],myCenter[1],myCenter[2]))
                tmpReg.initCont(rCont=1,rCent=myCenter,rMaxLen=rMaxLen)
                newGeom.add(tmpReg,what="reg")
                iHive=iHive+1
            for iT in range(1,len(thetas)):
                for iP in iPs:
                    tmpReg=Region()
                    tmpReg.rename("HVCL%04i"%(iHive),lNotify=False)
                    tmpReg.material=defMat
                    rDef='+%-8s -%-8s'%(spheres[iR].echoName(),spheres[iR-1].echoName())
                    if (TTs[iT]<=0.0):
                        tDef='+%-8s -%-8s'%(thetas [iT].echoName(),thetas [iT-1].echoName())
                    elif (TTs[iT-1]==0.0 or (TTs[iT]>0.0 and TTs[iT-1]<0.0)):
                        tDef='-%-8s -%-8s'%(thetas [iT].echoName(),thetas [iT-1].echoName())
                    elif (TTs[iT-1]>0.0):
                        tDef='-%-8s +%-8s'%(thetas [iT].echoName(),thetas [iT-1].echoName())
                    pDef='+%-8s -%-8s'%(phis[iP].echoName(),phis[iP-1].echoName())
                    tmpReg.addZone('%s %s %s'%(rDef,tDef,pDef))
                    myCenter=cellGrid.ret(what="POINT",iEl=iHive)
                    rMaxLen=1.1*max(RRs[iR]-RRs[iR-1],RRs[iR]*np.radians(TTs[iT]-TTs[iT-1]),RRs[iR]*np.radians(PPs[iT]-PPs[iT-1]))
                    tmpReg.tailMe("* - hive region %4d: R[cm]=[%g:%g], theta[deg]=[%g:%g], phi[deg]=[%g:%g]"%(
                        iHive,RRs[iR-1],RRs[iR],TTs[iT-1],TTs[iT],PPs[iP-1],PPs[iP]))
                    tmpReg.tailMe("*   center=[%g,%g,%g];"%(myCenter[0],myCenter[1],myCenter[2]))
                    tmpReg.initCont(rCont=1,rCent=myCenter,rMaxLen=rMaxLen)
                    newGeom.add(tmpReg,what="reg")
                    iHive=iHive+1
            if (cellGrid.HasPole("N")):
                tmpReg=Region()
                tmpReg.rename("HVCL%04i"%(iHive),lNotify=False)
                tmpReg.material=defMat
                tmpReg.addZone('+%-8s -%-8s +%-8s'%(spheres[iR].echoName(),spheres[iR-1].echoName(),thetas[-1].echoName()))
                myCenter=cellGrid.ret(what="POINT",iEl=iHive)
                rMaxLen=1.1*max(RRs[iR]-RRs[iR-1],RRs[iR]*2*np.absolute(np.radians(90-TTs[-1])))
                tmpReg.tailMe("* - hive region %4d: NORTH POLE! R[cm]=[%g:%g], theta[deg]=[%g:90]"%(
                    iHive,RRs[iR-1],RRs[iR],TTs[-1]))
                tmpReg.tailMe("*   center=[%g,%g,%g];"%(myCenter[0],myCenter[1],myCenter[2]))
                tmpReg.initCont(rCont=1,rCent=myCenter,rMaxLen=rMaxLen)
                newGeom.add(tmpReg,what="reg")
                iHive=iHive+1

        newGeom.headMe(tmpTitle)
        newGeom.setTitle(tmpTitle=tmpTitle)

        if (lWrapBHaround):
            newGeom=Geometry.WrapBH_Sphere(newGeom,2*max(RRs),3*max(RRs))
        return newGeom

    @staticmethod
    def BuildGriddedGeo(myGrid,myProtoList,myProtoGeos,osRegNames=[],lGeoDirs=False,lDebug=True):
        '''
        This method defines a FLUKA geometry representing a grid of objects.

        input parameters:
        - grid: an instance of a Grid() class, i.e. a list of locations;
        - myProtoList: this list states which prototype should be used at each
                       location. NB: len(myProtoList)=len(myGrid);
        - myProtoGeos: dictionary of actual prototype geometries. The unique
                       entries of myProtoList are the full set or a subset of
                       the keys of this dictionary;
        - osRegNames:  list of regions of the prototypes expressing the outermost
                       part of the prototypes, to be 'subtracted' from the region
                       definition of the hive cells, that should be sized by the
                       regions of the hive cell;
        '''
        myGeos=[]
        # loop over locations, to clone prototypes
        for iLoc,myLoc in enumerate(myGrid):
            if (myProtoList[iLoc] not in myProtoGeos):
                print("Geometry.BuildGriddedGeo(): unknown prototype %s!"%(\
                        myProtoList[iLoc]))
                exit(1)
            # - clone prototype
            myGeo=deepcopy(myProtoGeos[myProtoList[iLoc]])
            # - move clone to requested location/orientation
            #   NB: give priority to angles/axis wrt matrices, for higher
            #       numerical accuracy in final .inp file
            if (len(myLoc.ret("ANGLE"))>0):
                myGeo.solidTrasform(dd=myLoc.ret("POINT"),myTheta=myLoc.ret("ANGLE"),myAxis=myLoc.ret("AXIS"),lGeoDirs=lGeoDirs,lDebug=lDebug)
            else:
                myGeo.solidTrasform(dd=myLoc.ret("POINT"),myMat=myLoc.ret("MATRIX"),lGeoDirs=lGeoDirs,lDebug=lDebug)
            # - flag the region(s) outside the prototypes or that should be sized
            #   by the hive cells
            myGeo.flagRegs(osRegNames,-1,myLoc.ret("POINT"))
            # - rename the clone
            baseName="GR%03d"%(iLoc)
            myGeo.rename(baseName)
            # - notify the user about original prototype and location
            myGeo.headMe("GRID cell # %3d - family name: %s - prototype: %s\n"%(\
                            iLoc,baseName,myProtoList[iLoc])+ \
                         "* "+myLoc.echo(myFmt="% 13.6E",mySep="\n* ") )
            # - append clone to list of geometries
            myGeos.append(myGeo)
        # return merged geometry
        return Geometry.appendGeometries(myGeos)

    @staticmethod
    def MergeGeos(hiveGeo,gridGeo,lDebug=True,myTitle=None,prec=0.001):
        '''
        This method merges one FLUKA geometry onto another one.

        input parameters:
        - hiveGeo: Geometry instance of the hive;
        - gridGeo: Geometry instance of the grid of objects.
        - prec: precision of identification of points proximity [cm]

        The two geometries must not have common names - no check is performed
          for the time being.

        All the regions of gridGeo with rCont==-1 will be matched with regions
          of hiveGeo with rCont==1; a one-to-one mapping is established based
          on the rCent arrays. The merged regions will still belong to gridGeo
          and the respective region in hiveGeo will disappear.
        '''
        print("merging geometries...")
        iRhg=[ ii for ii,mReg in enumerate(hiveGeo.regs) if mReg.rCont==1 ]
        cRhg=np.array([ mReg.rCent for ii,mReg in enumerate(hiveGeo.regs) if mReg.rCont==1 ])
        ucRhg=np.unique(cRhg,axis=0)
        iRgg=[ ii for ii,mReg in enumerate(gridGeo.regs) if mReg.rCont==-1 ]
        cRgg=np.array([ mReg.rCent for ii,mReg in enumerate(gridGeo.regs) if mReg.rCont==-1 ])
        ucRgg=np.unique(cRgg,axis=0)
        if (lDebug):
            print("...found %d containing regions in hive:"%(len(iRhg)))
            print("   ...with %d unique centers!"%(len(ucRhg)))
            print(iRhg)
            print(cRhg)
            print(ucRhg)
            print("...found %d contained/sized regions in grid:"%(len(iRgg)))
            print("   ...with %d unique centers!"%(len(ucRgg)))
            print(iRgg)
            print(cRgg)
            print(ucRgg)
        if (len(iRhg)==len(ucRhg)):
            # for each location, only one hive region is concerned:
            #   merge each containing hive region into the concerned
            #   grid regions, and then remove the hive region
            # - merge region defs
            removeRegs=[]
            for jRhg in iRhg:
                for jRgg in iRgg:
                    if (np.linalg.norm(gridGeo.regs[jRgg].rCent-hiveGeo.regs[jRhg].rCent)<prec):
                        # resize infinite bodies
                        newL=hiveGeo.regs[jRhg].rMaxLen
                        whichBods=gridGeo.regs[jRgg].retBodiesInDef()
                        gridGeo.resizeBodies(newL,whichBods=whichBods)
                        # actually merge
                        gridGeo.regs[jRgg].merge(hiveGeo.regs[jRhg])
                        if (hiveGeo.regs[jRhg].echoName() not in removeRegs):
                            removeRegs.append(hiveGeo.regs[jRhg].echoName())
            # - remove merged regs
            for removeReg in removeRegs:
                myReg,iReg=hiveGeo.ret("REG",removeReg)
                hiveGeo.regs.pop(iReg)
        elif(len(iRgg)==len(ucRgg)):
            # for each location, only one grid region is concerned:
            #   merge each contained grid region into the concerned
            #   hive regions, and then remove the grid region
            # - merge region defs
            removeRegs=[]
            for jRgg in iRgg:
                for jRhg in iRhg:
                    if (np.linalg.norm(hiveGeo.regs[jRhg].rCent-gridGeo.regs[jRgg].rCent)<prec):
                        hiveGeo.regs[jRhg].merge(gridGeo.regs[jRgg])
                        if (gridGeo.regs[jRgg].echoName() not in removeRegs):
                            removeRegs.append(gridGeo.regs[jRgg].echoName())
            # - remove merged regs
            for removeReg in removeRegs:
                myReg,iReg=gridGeo.ret("REG",removeReg)
                gridGeo.regs.pop(iReg)
        else:
            print("...cannot merge more than a region of the hive and more than a region of the grid for a single location!")
            exit(1)
        print("...done.")
            
        return Geometry.appendGeometries([hiveGeo,gridGeo],myTitle=myTitle)

    @staticmethod
    def WrapBH_Sphere(myGeo,Rmin,Rmax,defMat="VACUUM",lDebug=True):
        '''
        Method for wrapping a spherical layer of blackhole around a given geometry
        '''
        print('wrapping a spherical layer of blackhole around geometry: Rmin=%g; Rmax=%g'%(Rmin,Rmax))
        newGeom=Geometry()

        print("...bodies...")
        bodies=[]
        for RR,tagName in zip([Rmin,Rmax],["inner","outer"]):
            tmpBD=Body()
            tmpBD.rename("BLK%s"%(tagName.upper()),lNotify=False)
            tmpBD.bType="SPH"
            tmpBD.Rs[0]=RR
            tmpBD.comment="* blackhole: %s radial boundary at R[cm]=%g"%(tagName,RR)
            bodies.append(tmpBD)
            
        print("...regions...")
        regions=[]
        # - regions outside / inside layer
        for iBod, (tagName,mySig,myPos) in enumerate(zip(["inner","outer"],["+","-"],["inside","outside"])):
            tmpReg=Region()
            tmpReg.rename("BLK%s"%(tagName.upper()),lNotify=False)
            tmpReg.material=defMat
            tmpReg.definition='''%s%-8s'''%(mySig,bodies[iBod].name)
            tmpReg.comment="* region %s blakchole layer"%(myPos)
            if (iBod==0):
                tmpReg.initCont(rCont=1)
            regions.append(tmpReg)
        # - actual layer
        tmpReg=Region()
        tmpReg.rename("BLKLAYER",lNotify=False)
        tmpReg.material="BLCKHOLE"
        tmpReg.definition='''+%-8s -%-8s'''%(bodies[-1].echoName(),bodies[0].echoName())
        tmpReg.comment="* blackhole layer"
        regions.append(tmpReg)

        newGeom.bods=bodies
        newGeom.regs=regions
        
        print('...done.')
        return Geometry.MergeGeos(newGeom,myGeo,lDebug=lDebug,myTitle=myGeo.title)

def acquireGeometries(fileNames,geoNames=None,lMakeRotatable=False):
    '''
    A simple function to parse a series of geometry files and
      store them in a dictionary of geometries.
    This function can be used to parse the database of prototype geometries.
    '''
    import os.path
    # check user input
    if (geoNames is None):
        geoNames=fileNames
    elif (len(fileNames)!=len(geoNames)):
        print("Number of items (%d) and names (%d) of geometries to acquire do not coincide!"%\
              (len(fileNames),len(geoNames)))
        exit(1)
        
    # acquire geometries:
    print("acquiring geometries...")
    myGeos={}
    for ii in range(len(fileNames)):
        if (not os.path.isfile(fileNames[ii]) or not os.path.exists(fileNames[ii])):
            print("something wrong with file %s! please check path, existence, access rights!"%(fileNames[ii]))
            exit(1)
        myGeos[geoNames[ii]]=Geometry.fromInp(fileNames[ii])
        if (geoNames[ii]!=fileNames[ii]):
            print("--> geometry saved in DB as %s;"%(geoNames[ii]))
    print("...acquired %d/%d geometries;"%(len(myGeos),len(fileNames)))

    # make geometries rotatable
    if (lMakeRotatable):
        print("making acquired geometries rotatable...")
        for myProtoName,myProtoGeo in myGeos.items():
            myProtoGeo.makeBodiesRotatable()
        print("...done;")
        
    return myGeos

def echoFloats(myFloats,expDigits=None,numDigist=None,prec=StringPrec,cut0s=cut0s):
    if (expDigits is None): expDigits=22
    if (numDigist is None): numDigist=expDigits-7
    numFmt="%% %d.%dE"%(expDigits,numDigist)
    if (type(myFloats) is list):
        floats2convert=myFloats
    else:
        floats2convert=[myFloats]
    myStrs=[]
    for float2convert in floats2convert:
        if (abs(float2convert)%1<prec): # integer
            myStrs.append("% .1f"%(float2convert))
        else: # actual floating point number
            myStr=numFmt%(float2convert)
            data=myStr.split("E")
            mantissa=data[0]; exponent=data[1]
            if (cut0s is not None): # truncate numbers if a series of "0" is found
                if (mantissa.find(cut0s)>-1):
                    mantissa=mantissa[:mantissa.index(cut0s)]
            if (exponent=="+00"):
                exponent=None
            elif (exponent=="+01"):
                exponent=None
                mantissa=mantissa.replace(".","")
                mantissa=mantissa[0:3]+"."+mantissa[3:]
            elif (exponent=="+02"):
                exponent=None
                mantissa=mantissa.replace(".","")
                mantissa=mantissa[0:4]+"."+mantissa[4:]
            elif (exponent=="-01"):
                exponent=None
                mantissa=mantissa.replace(".","")
                mantissa=mantissa[0]+"0."+mantissa[1:]
            elif (exponent=="-02"):
                exponent=None
                mantissa=mantissa.replace(".","")
                mantissa=mantissa[0]+"0.0"+mantissa[1:]
            if (mantissa=="-0.0"): mantissa=" 0.0"
            if (exponent is None):
                myStr=mantissa
            else:
                myStr="%sE%s"%(mantissa,exponent)
            myStrs.append(myStr)
    return myStrs

def assembleLine(myStrings,maxLen=MaxLineLength,header=LineHeader,lMultiLine=True):
    if (len(myStrings)==0):
        finStr=""
    else:
        finStr=myStrings[0]; curLen=len(finStr)
        for myStr in myStrings[1:]:
            if (curLen+len(myStr)+1>=maxLen):
                if (lMultiLine):
                    finStr=finStr+"\n"+header+myStr
                    curLen=len(header+myStr)
                else:
                    print("Cannot prepare a line longer than %d chars!"%(maxLen))
                    print("%s %s"%(finStr,myStr))
                    exit(1)
            else:
                finStr=finStr+" "+myStr
                curLen=curLen+len(myStr)+1
    return finStr

if (__name__=="__main__"):
    lDebug=False
    lGeoDirs=True
    # # - manipulate a geometry
    # caloCrysGeo=Geometry.fromInp("caloCrys.inp")
    # myMat=RotMat(myAng=60,myAxis=3,lDegs=True,lDebug=lDebug)
    # caloCrysGeo.solidTrasform(dd=[0,10,-20],myMat=myMat)
    # caloCrysGeo.echo("pippo.inp")
    
    # - test generation of geometry
    R=75
    dR=50
    NR=2
    Tmin=-20.0 # theta [degs] --> range: -Tmax:Tmax
    Tmax=20.0  # theta [degs] --> range: -Tmax:Tmax
    NT=4       # number of steps (i.e. grid cells)
    Pmin=-7.5
    Pmax=7.5   # phi [degs] --> range: -Pmax:Pmax
    NP=5       # number of steps (i.e. grid cells)
    # Tmax=3.0  # theta [degs] --> range: -Tmax:Tmax
    # NT=4      # number of steps (i.e. grid cells)
    # Pmax=2.0  # phi [degs] --> range: -Pmax:Pmax
    # NP=3      # number of steps (i.e. grid cells)

    # - hive geometry
    HiveGeo=Geometry.DefineHive_SphericalShell(R,R+dR,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,lDebug=lDebug)
    # HiveGeo.echo("hive.inp")

    # - gridded crystals
    #   acquire geometries
    fileNames=[ "caloCrys_02.inp" ] ; geoNames=fileNames
    myProtoGeos=acquireGeometries(fileNames,geoNames=geoNames,lMakeRotatable=not lGeoDirs);
    cellGrid=grid.SphericalShell(R,R+dR,NR,-Tmax,Tmax,NT,Pmin,Pmax,NP,lDebug=lDebug)
    myProtoList=[ "caloCrys_02.inp" for ii in range(len(cellGrid)) ]
    GridGeo=Geometry.BuildGriddedGeo(cellGrid,myProtoList,myProtoGeos,osRegNames=["OUTER"],lGeoDirs=lGeoDirs,lDebug=lDebug)
    # GridGeo.echo("grid.inp")

    # - merge geometries
    mergedGeo=Geometry.MergeGeos(HiveGeo,GridGeo,lDebug=lDebug)
    mergedGeo.reAssiginUSRBINunits(nMaxBins=35*35*10,usedUnits=26)
    mergedGeo.echo("merged.inp")
    
