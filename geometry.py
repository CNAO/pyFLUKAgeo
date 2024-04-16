# classes for managing FLUKA geometries
# A. Mereghetti, 2023/09/13
# python version: >= 3.8.10;

import numpy as np
from copy import deepcopy

import grid
from body import Body
from region import Region
from transformation import RotDefi, Transformation
from scorings import Usrbin

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
    - #define vars used only at parsing level, i.e. they are not stored
      in the parsed geometry;
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
            
    def ret(self,myKey,myValue):
        lFound=False
        if (myKey.upper().startswith("BODSINREG")):
            if (myValue.upper()=="ALL"):
                myReg,iReg=self.ret("reg","ALL")
                myEntry=[]; iEntry=[]
                for findReg in myReg:
                    tmpEl,tmpInd=self.ret("BodsInReg",findReg)
                    myEntry=myEntry+tmpEl
                    iEntry=iEntry+tmpInd
            else:
                myReg,iReg=self.ret("reg",myValue)
                myEntry=myReg.retBodiesInDef()
                iEntry=[]
                for findBod in myEntry:
                    tmpEl,tmpInd=self.ret("bod",findBod)
                    iEntry.append(tmpInd)
                lFound=True
        elif (myKey.upper().startswith("BOD")):
            if (myValue.upper()=="ALL"):
                myEntry=[ body.echoName() for body in self.bods ]
                iEntry=[ ii for ii in range(len(self.bods)) ]
                lFound=True
            else:
                for iEntry,myEntry in enumerate(self.bods):
                    if (myEntry.echoName()==myValue):
                        lFound=True
                        break
        elif (myKey.upper().startswith("REG")):
            if (myValue.upper()=="ALL"):
                myEntry=[ reg.echoName() for reg in self.regs ]
                iEntry=[ ii for ii in range(len(self.regs)) ]
                lFound=True
            else:
                for iEntry,myEntry in enumerate(self.regs):
                    if (myEntry.echoName()==myValue):
                        lFound=True
                        break
        elif (myKey.upper().startswith("TRANSF")):
            if (myValue.upper()=="ALL"):
                myEntry=[ tras.echoName() for tras in self.tras ]
                iEntry=[ ii for ii in range(len(self.tras)) ]
                lFound=True
            else:
                for iEntry,myEntry in enumerate(self.tras):
                    if (myEntry.echoName()==myValue):
                        lFound=True
                        break
        elif (myKey.upper().startswith("BININUNIT") or myKey.upper().startswith("USRBININUNIT")):
            myEntry=[] ; iEntry=[]
            if (isinstance(myValue,float) or isinstance(myValue,int)): myValue=[myValue]
            for ii,myBin in enumerate(self.bins):
                if (abs(myBin.getUnit()) in myValue):
                    myEntry.append(myBin)
                    iEntry.append(ii)
                    lFound=True
        elif (myKey.upper().startswith("BIN") or myKey.upper().startswith("USRBIN")):
            if (myValue.upper()=="ALL"):
                myEntry=[ myBin.echoName() for myBin in self.bins ]
                iEntry=[ ii for ii in range(len(self.bins)) ]
                lFound=True
            else:
                for iEntry,myEntry in enumerate(self.bins):
                    if (myEntry.echoName()==myValue):
                        lFound=True
                        break
        else:
            print("%s not recognised! What should I look for in the geometry?"%(myKey))
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
        defineFlags=[]
        print("parsing file %s..."%(myInpName))
        ff=open(myInpName,'r')
        iRead=0
        lFree=False
        lReads=[True]
        tmpBuf=""
        regBuf=""
        for tmpLine in ff.readlines():

            # PRE-PROCESSOR
            if (tmpLine.startswith("#define")):
                data=tmpLine.split()
                if (len(data)>2):
                    print("...cannot handle define vars with numerical value, only booleans!")
                    print(tmpLine)
                    exit(1)
                if (data[1] not in defineFlags):
                    defineFlags.append(data[1])
                continue
            elif (tmpLine.startswith("#if")):
                data=tmpLine.split()
                lReads.append(data[1] in defineFlags)
                continue
            elif (tmpLine.startswith("#elif")):
                data=tmpLine.split()
                lReads[-1]=(data[1] in defineFlags)
                continue
            elif (tmpLine.startswith("#else")):
                lReads[-1]=not lReads[-1]
                continue
            elif (tmpLine.startswith("#end")):
                lReads.pop()
                continue
            if (not all(lReads)): continue
            
            # OUTSIDE GEOBEGIN-GEOEND
            if (iRead==0):
                if (tmpLine.startswith("GEOBEGIN")):
                    iRead=1
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
            elif (iRead==1):
                # title after GEOBEGIN
                newGeom.title=tmpLine[20:].strip()
                iRead=2
                tmpBuf="" # flush buffer
            elif (iRead==2):
                # definition of FLUKA bodies
                if (tmpLine.startswith("END")):
                    print("...acquired %d bodies;"%(len(newGeom.bods)))
                    iRead=3
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
            elif (iRead==3):
                # definition of FLUKA regions
                if (tmpLine.startswith("END")):
                    if (len(regBuf)>0):
                        # acquire region
                        newGeom.add(Region.fromBuf(regBuf),what="reg")
                        regBuf="" # flush region def buffer
                    print("...acquired %d regions;"%(len(newGeom.regs)))
                    tmpBuf="" # flush buffer
                    iRead=0 # ignore LATTICE cards
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
                if (isinstance(myTheta,list) and isinstance(myAxis,list)):
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
        if (isinstance(whichRegs,str)):
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

    def resizeBodies(self,newL,whichBods="ALL",lDebug=False):
        '''
        input:
        - newL: new length [cm];
        - whichBods: list of body names to be updated (if infinite);
        '''
        if (isinstance(whichBods,str)):
            if (whichBods.upper()=="ALL"):
                bods2mod,iBods2mod=self.ret("bod","ALL")
                if (not isinstance(bods2mod,list)):
                    bods2mod, iBods2mod = [whichBods], [iBods2mod]
            else:
                bods2mod=[whichBods]
        else:
            if ( "ALL" in [tmpStr.upper() for tmpStr in whichBods] ):
                bods2mod,iBods2mod=self.ret("bod","ALL")
            else:
                bods2mod=whichBods
        if (lDebug): print("re-sizing %d bodies..."%(len(bods2mod)))
        for bod2mod in bods2mod:
            whichBod,iBod=self.ret("bod",bod2mod)
            whichBod.resize(newL)
        if (lDebug): print("...done;")

    def makeBodiesRotatable(self,lDebug=False,infL=1000.):
        for myBod in self.bods:
            myBod.makeRotatable(lDebug=lDebug,infL=infL)

    def reAssiginUSRBINunits(self,nMaxBins=None,nUSRBINs=None,usedUnits=None,lDebug=False):
        if (lDebug): print("re-assigning USRBIN units...")
        if (nMaxBins is not None and nUSRBINs is not None):
            print("Please tell me if I have to re-assign USRBIN units based on:")
            print("- max number of bins in a unit;")
            print("- max number of USRBINs in a unit;")
            exit(1)
        if (usedUnits is None): usedUnits=[]
        if (type(usedUnits) is not list): usedUnits=[usedUnits]
        uniqueUnits=list(set([ myBin.getUnit() for myBin in self.bins ]))
        uniqueUnits.sort(key=abs)
        if (lDebug):
            print("...%d original units:"%(len(uniqueUnits)),uniqueUnits)
            if (len(usedUnits)>0):
                print("...units NOT to be used:",usedUnits)
        currUnits=[21+ii for ii in range(len(uniqueUnits))]
        for ii in range(len(currUnits)):
            while (currUnits[ii] in [abs(iu) for iu in usedUnits]):
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
                while(myUnit<=max(currUnits) or \
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
        if (lDebug): print("...done;")

    def resizeUsrbins(self,newL,whichUnits=None,axis=3,lDebug=False):
        '''
        input:
        - newL: new length [cm];
        - whichUnits: units of USRBINs to be updated;
        '''
        if (whichUnits is None): whichUnits="ALL"
        if (isinstance(whichUnits,str)):
            if (whichUnits.upper()=="ALL"):
                bins2mod,iBins2mod=self.ret("bin","ALL")
                if (lDebug): print("re-sizing ALL USRBINs, i.e. %d ..."%(len(bins2mod)))
            else:
                print("Cannot specify USRBIN names for resizing!")
                exit(1)
        elif (isinstance(whichUnits,float) or isinstance(whichUnits,int)):
            bins2mod,iBins2mod=self.ret("BININUNIT",whichUnits)
            if (lDebug): print("re-sizing USRBINs in unit %d, i.e. %d ..."%(whichUnits,len(bins2mod)))
        elif (isinstance(whichUnits,list)):
            bins2mod=[]; iBins2mod=[]
            for whichUnit in whichUnits:
                tBins2mod,tIBins2mod=self.ret("BININUNIT",whichUnit)
                print("re-sizing USRBINs in unit %d, i.e. %d ..."%(whichUnit,len(tBins2mod)))
                if (isinstance(tBins2mod,list)):
                    bins2mod=bins2mod+tBins2mod
                else:
                    bins2mod.append(tBins2mod)
        else:
            print("Wrong indication of USRBINs for resizing!")
            print(whichUnits)
            exit(1)
        for bin2mod in bins2mod:
            bin2mod.resize(newL,axis=3)
        if (lDebug): print("...done;")

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
                rMaxLen=max(RRs[iR]-RRs[iR-1],RRs[iR]*2*np.absolute(np.radians(TTs[0]-90)))
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
                    rMaxLen=max(RRs[iR]-RRs[iR-1],RRs[iR]*np.radians(TTs[iT]-TTs[iT-1]),RRs[iR]*np.radians(PPs[iT]-PPs[iT-1]))
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
                rMaxLen=max(RRs[iR]-RRs[iR-1],RRs[iR]*2*np.absolute(np.radians(90-TTs[-1])))
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
        This method defines a list of FLUKA geometry representing a grid of objects.

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
        return myGeos

    @staticmethod
    def MergeGeos(hiveGeo,gridGeo,lDebug=True,myTitle=None,prec=0.001,enlargeFact=1.1,resBins="ALL"):
        '''
        This method merges one FLUKA geometry onto another one.

        input parameters:
        - hiveGeo: Geometry instance of the hive;
        - gridGeo: Geometry instance of the grid of objects;
        - prec: precision of identification of points proximity [cm];
        - enlargeFact: (safety) factor for scaling infinite bodies and USRBINs;
        - resBins: list of units of USRBINs that should be resized;

        The two geometries must not have common names - no check is performed
          for the time being.

        All the regions of gridGeo with rCont==-1 will be matched with regions
          of hiveGeo with rCont==1; a one-to-one mapping is established based
          on the rCent arrays. The merged regions will still belong to gridGeo
          and the respective region in hiveGeo will disappear.
        '''
        print("merging geometries...")
        if (isinstance(hiveGeo,Geometry)): hiveGeo=[hiveGeo]
        if (isinstance(gridGeo,Geometry)): gridGeo=[gridGeo]
        if (isinstance(enlargeFact,float) or isinstance(enlargeFact,int)):
            if (enlargeFact<=0): enlargeFact=None
        # hiveGeo
        iRhgs=[]; jRhgs=[]; cRhgs=[];
        for jj,myGeo in enumerate(hiveGeo):
            iRhgs=iRhgs+[ ii for ii,mReg in enumerate(myGeo.regs) if mReg.rCont==1 ]
            jRhgs=jRhgs+[ jj for mReg in myGeo.regs if mReg.rCont==1 ]
            cRhgs=cRhgs+[ mReg.rCent for mReg in myGeo.regs if mReg.rCont==1 ]
        cRhgs=np.array(cRhgs)
        ucRhgs=np.unique(cRhgs,axis=0)
        # gridGeo
        iRggs=[]; jRggs=[]; cRggs=[];
        for jj,myGeo in enumerate(gridGeo):
            iRggs=iRggs+[ ii for ii,mReg in enumerate(myGeo.regs) if mReg.rCont==-1 ]
            jRggs=jRggs+[ jj for mReg in myGeo.regs if mReg.rCont==-1 ]
            cRggs=cRggs+[ mReg.rCent for mReg in myGeo.regs if mReg.rCont==-1 ]
        cRggs=np.array(cRggs)
        ucRggs=np.unique(cRggs,axis=0)
        if (lDebug):
            print("...found %d containing regions in hive:"%(len(iRhgs)))
            print("   ...with %d unique centers!"%(len(ucRhgs)))
            print(iRhgs)
            print(jRhgs)
            print(cRhgs)
            print(ucRhgs)
            print("...found %d contained/sized regions in grid:"%(len(iRggs)))
            print("   ...with %d unique centers!"%(len(ucRggs)))
            print(iRggs)
            print(jRggs)
            print(cRggs)
            print(ucRggs)
        if (len(iRhgs)==len(ucRhgs)):
            # for each location, only one hive region is concerned:
            #   merge each containing hive region into the concerned
            #   grid regions, and then remove the hive region
            # - merge region defs
            removeRegs=[]; jRemoveRegs=[]
            for iRhg,jRhg in zip(iRhgs,jRhgs):
                for iRgg,jRgg in zip(iRggs,jRggs):
                    if (np.linalg.norm(gridGeo[jRgg].regs[iRgg].rCent-hiveGeo[jRhg].regs[iRhg].rCent)<prec):
                        # resize infinite bodies and USRBINs
                        if (enlargeFact is not None or resBins is not None):
                            newL=enlargeFact*hiveGeo[jRhg].regs[iRhg].rMaxLen
                        if (enlargeFact is not None): gridGeo[jRgg].resizeBodies(newL,lDebug=lDebug)
                        if (resBins is not None): gridGeo[jRgg].resizeUsrbins(newL,whichUnits=resBins,lDebug=lDebug)
                        # actually merge
                        gridGeo[jRgg].regs[iRgg].merge(hiveGeo[jRhg].regs[iRhg])
                        if (hiveGeo[jRhg].regs[iRhg].echoName() not in removeRegs):
                            removeRegs.append(hiveGeo[jRhg].regs[iRhg].echoName())
                            jRemoveRegs.append(jRhg)
            # - remove merged regs
            for removeReg,jRemoveReg in zip(removeRegs,jRemoveRegs):
                myReg,iReg=hiveGeo[jRemoveReg].ret("REG",removeReg)
                hiveGeo[jRemoveReg].regs.pop(iReg)
        elif(len(iRggs)==len(ucRggs)):
            # for each location, only one grid region is concerned:
            #   merge each contained grid region into the concerned
            #   hive regions, and then remove the grid region
            # - merge region defs
            removeRegs=[]; jRemoveRegs=[]
            for iRgg,jRgg in zip(iRggs,jRggs):
                for iRhg,jRhg in zip(iRhgs,jRhgs):
                    if (np.linalg.norm(gridGeo[jRgg].regs[iRgg].rCent-hiveGeo[jRhg].regs[iRhg].rCent)<prec):
                        hiveGeo[jRhg].regs[iRhg].merge(gridGeo[jRgg].regs[iRgg])
                        if (gridGeo[jRgg].regs[iRgg].echoName() not in removeRegs):
                            removeRegs.append(gridGeo[jRgg].regs[iRgg].echoName())
                            jRemoveRegs.append(jRgg)
            # - remove merged regs
            for removeReg in removeRegs:
                myReg,iReg=gridGeo[jRemoveReg].ret("REG",removeReg)
                gridGeo[jRemoveReg].regs.pop(iReg)
        else:
            print("...cannot merge more than a region of the hive and more than a region of the grid for a single location!")
            exit(1)
            
        return Geometry.appendGeometries(hiveGeo+gridGeo,myTitle=myTitle)

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

if (__name__=="__main__"):
    lDebug=False
    lGeoDirs=False
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
    
    (Geometry.appendGeometries(GridGeo)).echo("grid.inp")
    
    # - merge geometries
    mergedGeo=Geometry.MergeGeos(HiveGeo,GridGeo,lDebug=lDebug,resBins=25)
    mergedGeo.reAssiginUSRBINunits(nMaxBins=35*35*5000,usedUnits=26)
    mergedGeo.echo("merged.inp")
