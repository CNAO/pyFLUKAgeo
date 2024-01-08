# classes for managing FLUKA geometries
# A. Mereghetti, 2023/09/13
# python version: >= 3.8.10;

import numpy as np
from copy import deepcopy

import myMath
import grid

class GeoObject():
    '''
    A very basic object, implementing simply a name and a comment.
    This object is supposed to collect all the stuff concerning
      handling of names and modification of comments
    '''

    def __init__(self):
        self.name=""
        self.comment=""

    def rename(self,newName,lNotify=True):
        '''change instance name'''
        if (lNotify):
            self.tailMe("* NAME CHANGE: FROM %s TO %s"%(self.name,newName))
        self.name=newName

    def headMe(self,myString):
        '''head a string to the comment'''
        if (len(self.comment)>0):
            self.comment=myString+"\n"+self.comment
        else:
            self.comment=myString
            
    def tailMe(self,myString):
        '''tail a string to the comment'''
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
    - supported bodies: all planes and spheres, nothing else for the time being;
    - body definition on ONE line only, and always starts at col 1;
    - support only for preceding comments, NO comments headed by "!" or after
      body definition;
    - no support for $start_* cards;
    '''

    def __init__(self):
        '''
        default body is a FLUKA PLA
        '''
        GeoObject.__init__(self)
        self.bType="PLA"
        self.P=np.zeros(3)
        self.V=np.array([0.0,0.0,1.0])
        self.Rs=np.zeros(2)

    def echo(self):
        '''take into account comment'''
        if (self.bType=="PLA"):
            return GeoObject.echoComm(self)+ \
                "PLA %8s % 15.8E % 15.8E % 15.8E % 15.8E % 15.8E % 15.8E" % \
                ( self.echoName(), self.V[0], self.V[1], self.V[2], \
                             self.P[0], self.P[1], self.P[2] )
        elif (self.bType=="SPH"):
            return GeoObject.echoComm(self)+"SPH %8s % 15.8E % 15.8E % 15.8E % 15.8E" \
             % ( self.echoName(), self.P[0], self.P[1], self.P[2], self.Rs[0] )
        else:
            print("body %s NOT supported yet!"%(self.bType))
            exit(1)

    @staticmethod
    def fromBuf(tmpLines):
        newBody=Body()
        for tmpLine in tmpLines.splitlines():
            data=tmpLine.split()
            if (data[0]=="PLA"):
                newBody.bType="PLA"
                newBody.rename(data[1],lNotify=False)
                newBody.V=np.array(data[2:5]).astype(float)
                newBody.P=np.array(data[5:8]).astype(float)
            elif (data[0]=="YZP"):
                newBody.bType="PLA"
                newBody.rename(data[1],lNotify=False)
                newBody.V=np.array([1.0,0.0,0.0])
                newBody.P[0]=data[2].astype(float)
            elif (data[0]=="XZP"):
                newBody.bType="PLA"
                newBody.rename(data[1],lNotify=False)
                newBody.V=np.array([0.0,1.0,0.0])
                newBody.P[1]=data[2].astype(float)
            elif (data[0]=="XYP"):
                newBody.bType="PLA"
                newBody.rename(data[1],lNotify=False)
                newBody.V=np.array([0.0,0.0,1.0])
                newBody.P[2]=data[2].astype(float)
            elif (data[0]=="SPH"):
                newBody.bType="SPH"
                newBody.rename(data[1],lNotify=False)
                newBody.P=np.array(data[2:5]).astype(float)
                newBody.Rs[0]=data[5].astype(float)
            elif (tmpLine.startswith("*")):
                newBody.tailMe(tmpLine)
            else:
                print("body %s NOT supported yet!"%(data[0]))
                exit(1)
        return newBody

    def traslate(self,dd=None):
        if (dd is not None):
            if (self.bType=="PLA" or self.bType=="SPH" ):
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

    def __init__(self):
        GeoObject.__init__(self)
        self.neigh=5
        self.definition=""
        self.material="BLACKHOLE"
        # additional fields
        self.initCont()

    def initCont(self,rCont=0,rCent=np.array([0.0,0.0,0.0])):
        self.rCont=rCont # containment indicator (1 per region):
                         # -1: region to be contained into/sized by another one
                         #  0: regular region (neither contains nor it is contained)
                         #  1: cell region (it contains another region)
        self.rCent=rCent # central point of one or more region (3D array)

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
        '''
        simple method to query-replace body names in region definition
        '''
        for oName,nName in zip(oldNames,newNames):
            self.definition=self.definition.replace(oName,nName)

    def assignMat(self,myMaterial):
        self.material=myMaterial
        
    def echo(self,lMat=False):
        '''take into account comment'''
        if (lMat):
            # echo ASSIGNMA card
            return GeoObject.echoComm(self)+"ASSIGNMA  %10s%10s" % ( self.material, self.echoName() )
        else:
            return GeoObject.echoComm(self)+"%-8s   %4d %s" % \
                ( self.echoName(), self.neigh, self.definition )

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
    Please have a look also to the Transformation class: the ROT-DEFI
       class simply stores infos of a specific rotation/traslation;
    - echo in FREE format uses an empty space as field separator;
    - parsing in FREE format expects an empty space as field separator;
    '''

    def __init__(self):
        GeoObject.__init__(self) # mainly for comments, name NOT used!
        self.axis=3    # FLUKA coding: 1=x, 2=y, 3=z
        self.phi=0.0   # polar angle [degrees] (FLUKA units)
        self.theta=0.0 # azimuth angle [degrees] (FLUKA units)
        self.DD=np.zeros(3) # translation [cm] (x,y,z components)

    def echo(self,lFree=True,myID=0,myName=""):
        '''echo in FREE format uses an empty space as field separator'''
        if (myID<=0):
            print("...cannot dump a ROT-DEFI card without an index!")
            exit(1)
        if (len(myName)==0):
            print("...cannot dump a ROT-DEFI card without a name!")
            exit(1)
        if (myID<100):
            myIDecho=self.axis*100+myID
        else:
            myIDecho=self.axis+myID*1000
        if (lFree):
            return GeoObject.echoComm(self)+ \
                "ROT-DEFI  %10.1f % 15.8E % 15.8E % 15.8E % 15.8E % 15.8E %-10s"\
                %( myIDecho, self.phi, self.theta, self.DD[0], self.DD[1], self.DD[2], myName)
        else:
            return GeoObject.echoComm(self)+ \
                "ROT-DEFI  %10.1f% 10G% 10G% 10G% 10G% 10G%-10s"\
                %( myIDecho, self.phi, self.theta, self.DD[0], self.DD[1], self.DD[2], myName)
            
    def fromBuf(self,myBuffer,lFree=True):
        '''parsing in FREE format expects an empty space as field separator'''
        for tmpLine in myBuffer.splitlines():
            if (tmpLine.startswith("*")):
                self.tailMe(tmpLine)
            else:
                # parse data
                if (lFree):
                    data=tmpLine.split()
                    pID=float(data[1])
                    phi=float(data[2])
                    tht=float(data[3])
                    DD=data[4:7]
                    myName=data[7]
                else:
                    pID=float(tmpLine[10:20].strip())
                    phi=float(tmpLine[20:30].strip())
                    tht=float(tmpLine[30:40].strip())
                    DD=[tmpLine[40:50].strip(),tmpLine[50:60].strip(),tmpLine[60:70].strip()]
                    myName=tmpLine[70:80].strip()
                # store data
                if (pID>1000):
                    myID=int(pID/1000)
                    self.axis=pID%1000
                else:
                    self.axis=int(pID/100)
                    myID=pID%100
                self.phi=phi
                self.theta=tht
                self.DD=np.array(DD)
        return myID, myName

class Transformation(GeoObject):
    '''
    This is the actual class to be used when dealing with ROT-DEFI cards.
    The class implements an array of ROT-DEFI cards - hence, at least one card.
    '''

    def __init__(self):
        GeoObject.__init__(self) # mainly for names, comments NOT used!
        self.RotDefis=[] # list of ROT-DEFI cards
        self.ID=0
        
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
        self.RotDefis[0].echoComm()

    def AddRotDefi(self,newRotDefi):
        self.RotDefis.append(newRotDefi)

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
                        newTransf.ID=myID
                    elif(newTransf.ID!=myID):
                        print("...cannot add ROT-DEFI with ID=%d to tranformation ID=%d!"%(myID,newTransf.ID))
                        exit(1)
                if (len(myName)>0):
                    if (len(newTransf.echoName())==0):
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
            buf=buf+myRotDefi.echo(myID=self.ID,myName=self.name)+"\n"
        return buf
                        
class Geometry():
    '''
    - name-based FLUKA geometry defition;
    - NO support of LATTICE cards;
    - NO support for #include cards or geo defitions in files other than that
       being parsed;
    - comments:
      . body: commented lines are considered always before the comment, and only
              commented lines before the body will be retained;
              --> trailing comments to body declaration section will disappear!
      . region: commented lines are kept where they are found, if they are
                found before or along the region declaration;
              --> trailing comments to region declaration will disappear!
    '''
    def __init__(self):
        self.bods=[]
        self.regs=[]
        self.tras=[]
        self.title=""

    def add(self,tmpEl,what="body"):
        '''
        tmpEl is already the instance to be added, NOT the string buffer to be parsed
        '''
        if (what.upper().startswith("BOD")):
            self.bods.append(tmpEl)
        elif (what.upper().startswith("REG")):
            self.regs.append(tmpEl)
        elif (what.upper().startswith("TRANSF")):
            self.tras.append(tmpEl)
        else:
            print("...what do you want to add to geometry? %s NOT recognised!"%(what))
            exit(1)

    def setTitle(self,tmpTitle="My Geometry"):
        self.title=tmpTitle
        
    def assignma(self,tmpLine):
        data=tmpLine.split()
        myMat=data[1]
        myReg=data[2]
        print("...assigning material %s to region %s..."%(myMat,myReg))
        myEntry,iEntry=self.ret(myReg,what="REG")
        self.regs[iEntry].assignMat(myMat)

    def headMe(self,myString,begLine="* \n"+"* "+"="*108,endLine="* "+"-"*108+"\n* "):
        '''
        simple method to head a string to the geometry declaration (bodies,
           regions, assignma cards)
        '''
        actualString=begLine+"\n* "+myString+" \n"+endLine
        if (len(self.bods)>0):
            self.bods[0].headMe(actualString)
        if (len(self.regs)>0):
            self.regs[0].headMe(actualString)
        if (len(self.tras)>0):
            self.tras[0].headMe(actualString)
            
    def ret(self,myWhat,myName):
        lFound=False
        if (myWhat.upper().startswith("BOD")):
            if (myName.upper()=="ALL"):
                myEntry=[ body.echoName() for body in self.bods ]
                iEntry=[ ii for ii in range(len(self.bods)) ]
            else:
                for iEntry,myEntry in enumerate(self.bods):
                    if (myEntry.echoName()==myName):
                        lFound=True
                        break
        elif (myWhat.upper().startswith("REG")):
            if (myName.upper()=="ALL"):
                myEntry=[ reg.echoName() for reg in self.regs ]
                iEntry=[ ii for ii in range(len(self.regs)) ]
            else:
                for iEntry,myEntry in enumerate(self.regs):
                    if (myEntry.echoName()==myName):
                        lFound=True
                        break
        elif (what.upper().startswith("TRANSF")):
            if (myName.upper()=="ALL"):
                myEntry=[ tras.echoName() for tras in self.tras ]
                iEntry=[ ii for ii in range(len(self.tras)) ]
            else:
                for iEntry,myEntry in enumerate(self.tras):
                    if (myEntry.echoName()==myName):
                        lFound=True
                        break
        else:
            print("%s not recognised! What should I look for in the geometry?"%(myWhat))
            exit(1)
        if (not lFound):
            print("unable to find %s named %s in geometry..."%(myWhat,myName))
            exit(1)
        return myEntry, iEntry

    @staticmethod
    def fromInp(myInpName):
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
                elif (tmpLine.startswith("ROT-DEFI")):
                    # acquire ROT-DEFI card
                    tmpRotDefi=RotDefi()
                    myID,myName=tmpRotDefi.fromBuf(tmpBuf,lFree=lFree)
                    myEntry,iEntry=self.ret("TRANSF",myName)
                    if (iEntry==-1):
                        # brand new transformation
                        self.add(tmpRotDefi,"TRANSF")
                    else:
                        self.tras[iEntry].AddRotDefi(tmpRotDefi)
                    tmpBuf="" # flush buffer
                elif (tmpLine.startswith("*")):
                    if (len(tmpBuf)>0):
                        tmpBuf=tmpBuf+"\n"+tmpLine
                    else:
                        tmpBuf=tmpLine
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
                        newGeom.add(Body.fromBuf(tmpBuf.strip()))
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
                            newGeom.add(Region.fromBuf(regBuf),what="reg"))
                            regBuf="" # flush region def buffer
                        regBuf=tmpBuf+tmpLine
                        tmpBuf="" # flush buffer
                    
        ff.close()
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
            new.tras=new.regs+myGeo.tras
        if (myTitle is None):
            myTitle="appended geometries"
        new.title=myTitle
        return new

    def echo(self,oFileName,lSplit=False,what="all",dMode="w"):
        '''
        - what="all"/"bodies"/"regions"/"materials"
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
        elif (what.upper()=="ALL"):
            if (lSplit):
                if (oFileName.endswith(".inp")):
                    # split geometry definition into bodies,
                    #   regions and assignmat files, to be used with
                    #   pure #include statements;
                    self.echo(oFileName.replace(".inp","_bodies.inp",1),\
                              lSplit=False,what="bodies",dMode="w")
                    self.echo(oFileName.replace(".inp","_regions.inp",1),\
                              lSplit=False,what="regions",dMode="w")
                    self.echo(oFileName.replace(".inp","_assignmats.inp",1),\
                              lSplit=False,what="materials",dMode="w")
                    if (len(self.tras)>0):
                        self.echo(oFileName.replace(".inp","_rotdefis.inp",1),\
                              lSplit=False,what="transf",dMode="w")
                else:
                    # split geometry definition into a .geo file
                    #   and an assignmat file; the former is encapsulated
                    #   between GEOBEGIN and GEOEND cards, the other is
                    #   imported via an #include statement
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
        else:
            print("...what should I echo? %s NOT reconised!"%(what))

    def solidTrasform(self,dd=None,myMat=None,myTheta=None,myAxis=3,lDegs=True,lDebug=False):
        print("applying solid transformation(s)...")
        if (myMat is not None):
            print("...applying transformation expressed by matrix to geometry...")
            for ii in range(len(self.bods)):
                self.bods[ii].rotate(myMat=myMat,myTheta=None,myAxis=None,\
                                     lDegs=lDegs,lDebug=lDebug)
        elif (myTheta is not None):
            print("...applying rotation by %f degs around axis %d..."%\
                  (myTheta,myAxis))
            for ii in range(len(self.bods)):
                self.bods[ii].rotate(myMat=None,myTheta=myTheta,myAxis=myAxis,\
                                     lDegs=lDegs,lDebug=lDebug)
        if (dd is not None):
            print("...applying traslation array [%f,%f,%f] cm..."%\
                  (dd[0],dd[1],dd[2]))
            for ii in range(len(self.bods)):
                self.bods[ii].traslate(dd=dd)
        if (myMat is None and myTheta is None and dd is None):
            print("...no transformation provided!")
        print("...done.")

    def rename(self,newName,lNotify=True):
        maxLenName=8
        if (len(newName)>=maxLenName):
            print("Geometry.rename(): cannot rename entities with len(%s)>=%d!"%(newName,maxLenName))
            exit(1)
        newNameFmt=newName+"%0"+"%d"%(maxLenName-len(newName))+"d"
        oldBodyNames=[]; newBodyNames=[]
        for iBody in range(len(self.bods)):
            oldBodyNames.append(self.bods[iBody].echoName())
            newBodyNames.append(newNameFmt%(iBody+1))
            self.bods[iBody].rename(newBodyNames[-1],lNotify=lNotify)
        for iReg in range(len(self.regs)):
            self.regs[iReg].rename(newNameFmt%(iReg+1),lNotify=lNotify)
            self.regs[iReg].BodyNameReplaceInDef(oldBodyNames,newBodyNames)
        for iTras in range(len(self.tras)):
            self.tras[iTras].rename(newNameFmt%(iTras+1),lNotify=lNotify)

    def flagRegs(self,whichRegs,rCont,rCent):
        if (whichRegs is str):
            if (whichRegs.upper()=="ALL"):
                regs2mod,iRegs2mod=myGeo.ret("reg","ALL")
            else:
                regs2mod=[whichRegs]
        else:
            if ( "ALL" in [tmpStr.upper() for tmpStr in whichRegs] ):
                regs2mod,iRegs2mod=myGeo.ret("reg","ALL")
            else:
                regs2mod=whichRegs
        for whichReg in whichRegs:
            outReg,iReg=self.ret("reg",whichReg)
            outReg.initCont(rCont=rCont,rCent=rCent)

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
        - on theta, by rotated XZPs;
        - on phi, by rotated YZPs;
        '''
        
        print("Preparing the hive for a spherical shell...")
        
        print("...generating the grid of cells...")
        cellGrid=grid.Grid.SphericalShell(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,lDebug=lDebug)
        print("...defining hive boundaries...")
        myHive=grid.SphericalHive(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,lDebug=lDebug)
        RRs,TTs,PPs=myHive.ret(myWhat="all")

        print("...generating FLUKA geometry...")
        newGeom=Geometry()
        
        print("   ...bodies...")
        spheres=[]
        for ii,RR in enumerate(RRs,1):
            tmpBD=Body()
            tmpBD.rename("HVRAD%03i"%(ii),lNotify=False)
            tmpBD.bType="SPH"
            tmpBD.Rs[0]=RR
            tmpBD.comment="* Hive radial boundary at R[cm]=%g"%(RR)
            spheres.append(tmpBD)
        spheres[0].comment="* \n"+spheres[0].comment
        thetas=[]
        for ii,TT in enumerate(TTs,1):
            tmpBD=Body()
            tmpBD.rename("HVTHT%03i"%(ii),lNotify=False)
            tmpBD.V=np.array([0.0,1.0,0.0])
            tmpBD.rotate(myMat=None,myTheta=-TT,myAxis=1,lDegs=True,lDebug=lDebug)
            tmpBD.comment="* Hive theta boundary at theta[deg]=%g"%(TT)
            thetas.append(tmpBD)
        thetas[0].comment="* \n"+thetas[0].comment
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
        tmpReg.definition='''-%-8s'''%(spheres[-1].name)
        tmpReg.comment="* region outside hive"
        tmpReg.initCont(rCont=-1)
        newGeom.add(tmpReg,what="reg"))
        # - inside hive
        tmpReg=Region()
        tmpReg.rename("HV_INNER",lNotify=False)
        tmpReg.material=defMat
        tmpReg.definition='''+%-8s'''%(spheres[0].name)
        tmpReg.comment="* region inside hive"
        newGeom.add(tmpReg,what="reg"))
        # - around hive
        tmpReg=Region()
        tmpReg.rename("HVAROUND",lNotify=False)
        tmpReg.material=defMat
        tmpReg.definition='''| +%-8s -%-8s -%-8s
                | +%-8s -%-8s +%-8s
                | +%-8s -%-8s +%-8s -%-8s +%-8s
                | +%-8s -%-8s +%-8s -%-8s -%-8s'''%( \
     spheres[-1].echoName(),spheres[0].echoName(), thetas[-1].echoName(), \
     spheres[-1].echoName(),spheres[0].echoName(), thetas[ 0].echoName(), \
     spheres[-1].echoName(),spheres[0].echoName(), thetas[-1].echoName(), thetas[ 0].echoName(), phis[ 0].echoName(), \
     spheres[-1].echoName(),spheres[0].echoName(), thetas[-1].echoName(), thetas[ 0].echoName(), phis[-1].echoName()  )
        tmpReg.comment="* region around hive"
        newGeom.add(tmpReg,what="reg"))
        # - actual hive
        iHive=0
        for iR in range(1,len(spheres)):
            for iT in range(1,len(thetas)):
                for iP in range(1,len(phis)):
                    tmpReg=Region()
                    tmpReg.rename("HVCL%04i"%(iHive),lNotify=False)
                    tmpReg.material=defMat
                    tmpReg.definition='+%-8s -%-8s +%-8s -%-8s +%-8s -%-8s'%\
                        (spheres[iR].echoName(),spheres[iR-1].echoName(),\
                         thetas [iT].echoName(),thetas [iT-1].echoName(),\
                         phis   [iP].echoName(),phis   [iP-1].echoName())
                    tmpComment=             "* - hive region %4d: R[cm]=[%g:%g], theta[deg]=[%g:%g], phi[deg]=[%g:%g]"%(
                        iHive,RRs[iR-1],RRs[iR],TTs[iT-1],TTs[iT],PPs[iP-1],PPs[iP])
                    myCenter=cellGrid.ret(myWhat="POINT",iEl=iHive)
                    tmpComment=tmpComment+"\n*   center=[%g,%g,%g];"%(myCenter[0],myCenter[1],myCenter[2])
                    tmpReg.comment=tmpComment
                    tmpReg.initCont(rCont=1,rCent=myCenter)
                    newGeom.add(tmpReg,what="reg"))
                    iHive=iHive+1

        newGeom.headMe(tmpTitle)
        newGeom.setTitle(tmpTitle=tmpTitle)

        if (lWrapBHaround):
            newGeom=Geometry.WrapBH_Sphere(newGeom,2*max(RRs),3*max(RRs))
        return newGeom

    @staticmethod
    def BuildGriddedGeo(myGrid,myProtoList,myProtoGeos,osRegNames=[],lDebug=True):
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
                print("Geometry.BuildGriddedGeo_SphericalShell(): unknown prototype %s!"%(\
                        myProtoList[iLoc]))
                exit(1)
            # - clone prototype
            myGeo=deepcopy(myProtoGeos[myProtoList[iLoc]])
            # - move clone to requested location/orientation
            myGeo.solidTrasform(dd=myLoc.ret("POINT"),myMat=myLoc.ret("MATRIX"),lDebug=lDebug)
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

def acquireGeometries(fileNames,geoNames=None):
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
    return myGeos

if (__name__=="__main__"):
    lDebug=False
    # # - manipulate a geometry
    # caloCrysGeo=Geometry.fromInp("caloCrys.inp")
    # myMat=RotMat(myAng=60,myAxis=3,lDegs=True,lDebug=lDebug)
    # caloCrysGeo.solidTrasform(dd=[0,10,-20],myMat=myMat)
    # caloCrysGeo.echo("pippo.inp")
    
    # - test generation of geometry
    R=100
    dR=50
    NR=1
    Tmax=14.0  # theta [degs] --> range: -Tmax:Tmax
    NT=15      # number of steps (i.e. grid cells)
    Pmax=23.0  # phi [degs] --> range: -Pmax:Pmax
    NP=24      # number of steps (i.e. grid cells)
    # Tmax=3.0  # theta [degs] --> range: -Tmax:Tmax
    # NT=4      # number of steps (i.e. grid cells)
    # Pmax=2.0  # phi [degs] --> range: -Pmax:Pmax
    # NP=3      # number of steps (i.e. grid cells)

    # - hive geometry
    HiveGeo=Geometry.DefineHive_SphericalShell(R,R+dR,NR,-Tmax,Tmax,NT,-Pmax,Pmax,NP,lDebug=lDebug)
    # HiveGeo.echo("hive.inp")

    # - gridded crystals
    #   acquire geometries
    fileNames=[ "caloCrys_01.inp" ] ; geoNames=fileNames
    myProtoGeos=acquireGeometries(fileNames,geoNames=geoNames);
    cellGrid=grid.Grid.SphericalShell(R,R+dR,NR,-Tmax,Tmax,NT,-Pmax,Pmax,NP,lDebug=lDebug)
    myProtoList=[ "caloCrys_01.inp" for ii in range(len(cellGrid)) ]
    GridGeo=Geometry.BuildGriddedGeo(cellGrid,myProtoList,myProtoGeos,osRegNames=["OUTER"],lDebug=lDebug)
    # GridGeo.echo("grid.inp")

    # - merge geometries
    mergedGeo=Geometry.MergeGeos(HiveGeo,GridGeo,lDebug=lDebug)
    mergedGeo.echo("merged.inp")
    
