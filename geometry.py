# classes for managing FLUKA geometries
# A. Mereghetti, 2023/09/13
# python version: >= 3.8.10;

import numpy as np
from copy import deepcopy

import myMath
import grid

class Body:
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
        self.bName=""
        self.bType="PLA"
        self.P=np.zeros(3)
        self.V=np.array([0.0,0.0,1.0])
        self.Rs=np.zeros(2)
        self.comment=""

    def echo(self):
        # take into account comment lines
        myBuf=""
        if (len(self.comment)>0):
            myBuf=self.comment+"\n"
        # actual body definition
        if (self.bType=="PLA"):
            return myBuf+ \
                "PLA %8s % 15.8E % 15.8E % 15.8E % 15.8E % 15.8E % 15.8E" % \
                ( self.bName, self.V[0], self.V[1], self.V[2], \
                              self.P[0], self.P[1], self.P[2] )
        elif (self.bType=="SPH"):
            return myBuf+"SPH %8s % 15.8E % 15.8E % 15.8E % 15.8E" \
             % ( self.bName, self.P[0], self.P[1], self.P[2], self.Rs[0] )
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
                newBody.bName=data[1]
                newBody.V=np.array(data[2:5]).astype(float)
                newBody.P=np.array(data[5:8]).astype(float)
            elif (data[0]=="YZP"):
                newBody.bType="PLA"
                newBody.bName=data[1]
                newBody.V=np.array([1.0,0.0,0.0])
                newBody.P[0]=data[2].astype(float)
            elif (data[0]=="XZP"):
                newBody.bType="PLA"
                newBody.bName=data[1]
                newBody.V=np.array([0.0,1.0,0.0])
                newBody.P[1]=data[2].astype(float)
            elif (data[0]=="XYP"):
                newBody.bType="PLA"
                newBody.bName=data[1]
                newBody.V=np.array([0.0,0.0,1.0])
                newBody.P[2]=data[2].astype(float)
            elif (data[0]=="SPH"):
                newBody.bType="SPH"
                newBody.bName=data[1]
                newBody.P=np.array(data[2:5]).astype(float)
                newBody.Rs[0]=data[5].astype(float)
            elif (tmpLine.startswith("*")):
                if (len(newBody.comment)==0):
                    newBody.comment=tmpLine
                else:
                    newBody.comment=newBody.comment+"\n"+tmpLine
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
                
    def rename(self,newName,lNotify=True):
        if (lNotify):
            self.tailMe("* NAME CHANGE: FROM %s TO %s"%(self.bName,newName))
        self.bName=newName

    def headMe(self,myString):
        '''
        simple method to head a string to the comment of the body declaration
        '''
        if (len(self.comment)>0):
            self.comment=myString+"\n"+self.comment
        else:
            self.comment=myString
            
    def tailMe(self,myString):
        '''
        simple method to tail a string to the comment of the body declaration
        '''
        if (len(self.comment)>0):
            self.comment=self.comment+"\n"+myString
        else:
            self.comment=myString
                
class Region():
    '''
    - no parsing/splitting of zones;
    - support only for preceding comments or commented lines between region
      definition lines, NO comments headed by !
    - a region definition always starts at column 1;
    - no check of length of lines for region definition;
    - ASSIGNMA cards: only 1:1 material:region assignment, NO material
      assignment for decay radiation simulation, NO magnetic/electric fields;
    '''

    def __init__(self):
        self.rName=""
        self.neigh=5
        self.definition=""
        self.material="BLACKHOLE"
        self.comment=""

    @staticmethod
    def fromBuf(myBuffer):
        newReg=Region()
        lHeadParsed=False
        for tmpLine in myBuffer.splitlines():
            if (not lHeadParsed):
                if (tmpLine.startswith("*")):
                    if (len(newReg.comment)==0):
                        newReg.comment=tmpLine
                    else:
                        newReg.comment=newReg.comment+"\n"+tmpLine
                else:
                    data=tmpLine.split()
                    newReg.rName=data[0]
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
        
    def rename(self,newName,lNotify=True):
        if (lNotify):
            self.tailMe("* NAME CHANGE: FROM %s TO %s"%(self.rName,newName))
        self.rName=newName

    def headMe(self,myString):
        '''
        simple method to head a string to the comment of the region declaration
        '''
        if (len(self.comment)>0):
            self.comment=myString+"\n"+self.comment
        else:
            self.comment=myString
            
    def tailMe(self,myString):
        '''
        simple method to tail a string to the comment of the region declaration
        '''
        if (len(self.comment)>0):
            self.comment=self.comment+"\n"+myString
        else:
            self.comment=myString
                
    def BodyNameReplaceInDef(self,oldNames,newNames):
        '''
        simple method to query-replace body names in region definition
        '''
        for oName,nName in zip(oldNames,newNames):
            self.definition=self.definition.replace(oName,nName)

    def assignMat(self,myMaterial):
        self.material=myMaterial
        
    def echo(self,lMat=False):
        # take into account comment lines
        myBuf=""
        if (len(self.comment)>0):
            myBuf=self.comment+"\n"
        if (lMat):
            # echo ASSIGNMA card
            return myBuf+"ASSIGNMA  %10s%10s" % ( self.material, self.rName )
        else:
            return myBuf+"%-8s   %4d %s" % \
                ( self.rName, self.neigh, self.definition )

    def merge(self,newReg,spacing=" "*17):
        # warn user in comment
        self.tailMe("* --> merged with region %s <--"%(newReg.rName))
        # merge definitions
        mergedDef=""
        newRegZones=newReg.definition.split("|")
        myRegZones=self.definition.split("|")
        lFirst=True
        for myRegZone in myRegZones:
            sMyRegZone=myRegZone.strip()
            if (len(sMyRegZone)>0):
                for newRegZone in newRegZones:
                    sNewRegZone=newRegZone.strip()
                    if (len(sNewRegZone)>0):
                        if (lFirst):
                            mergedDef="| %s %s"%(sMyRegZone,sNewRegZone)
                            lFirst=False
                        else:
                            mergedDef=mergedDef+"\n%s| %s %s"%(spacing,sMyRegZone,sNewRegZone)
        self.definition=mergedDef
        # merge comments
        if (len(new.comment)>0):
            self.taileMe(new.comment)
        # merge neighbours
        self.neigh=self.neigh+newReg.neigh

class Geometry():
    '''
    - name-based FLUKA geometry defition;
    - NO support of LATTICE cards;
    - NO support for #input cards or geo defitions in files other than that
       being parsed;
    - NO support for ROT-DEFI cards;
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
        self.title=""

    def addBod(self,tmpBod):
        self.bods.append(tmpBod)
 
    def addReg(self,tmpReg):
        self.regs.append(tmpReg)

    def setTitle(self,tmpTitle="My Geometry"):
        self.title=tmpTitle
        
    def assignma(self,tmpLine):
        data=tmpLine.split()
        myMat=data[1]
        myReg=data[2]
        print("...assigning material %s to region %s..."%(myMat,myReg))
        lFound=False
        for ii in range(len(self.regs)):
            if (self.regs[ii].rName==myReg):
                lFound=True
                self.regs[ii].assignMat(myMat)
        if (not lFound):
            print("...region %s NOT found in geometry!"%(myReg))
            exit(1)

    def headMe(self,myString):
        '''
        simple method to head a string to the geometry declaration (bodies,
           regions, assignma cards)
        '''
        self.bods[0].headMe(myString)
        self.regs[0].headMe(myString)
            
    def ret(self,myWhat,myName):
        lFound=False
        if (myWhat.upper()=="BODY"):
            if (myName.upper()=="ALL"):
                myEntry=[ body.bName for body in self.bos ]
                iEntry=[ ii for ii in range(len(self.bods)) ]
            else:
                for iEntry,myEntry in enumerate(self.bods):
                    if (myEntry.bName==myName):
                        lFound=True
                        break
        elif (myWhat.upper()=="REGION"):
            if (myName.upper()=="ALL"):
                myEntry=[ reg.rName for reg in self.regs ]
                iEntry=[ ii for ii in range(len(self.regs)) ]
            else:
                for iEntry,myEntry in enumerate(self.regs):
                    if (myEntry.rName==myName):
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
        tmpBuf=""
        regBuf=""
        for tmpLine in ff.readlines():
            if (lRead==0):
                # non-geometry input
                if (tmpLine.startswith("GEOBEGIN")):
                    lRead=1
                elif (tmpLine.startswith("ASSIGNMA")):
                    newGeom.assignma(tmpLine)
            elif (lRead==1):
                # title after GEOBEGIN
                newGeom.title=tmpLine[20:].strip()
                lRead=2
            elif (lRead==2):
                # definition of FLUKA bodies
                if (tmpLine.startswith("END")):
                    print("...acquired %d bodies;"%(len(newGeom.bods)))
                    lRead=3
                    tmpBuf="" # flush buffer
                else:
                    tmpBuf=tmpBuf+tmpLine
                    if (not tmpLine.startswith("*")):
                        # acquire body
                        newGeom.addBod(Body.fromBuf(tmpBuf.strip()))
                        tmpBuf="" # flush buffer
            elif (lRead==3):
                # definition of FLUKA regions
                if (tmpLine.startswith("END")):
                    if (len(regBuf)>0):
                        # acquire region
                        newGeom.addReg(Region.fromBuf(regBuf))
                        regBuf="" # flush region def buffer
                    print("...acquired %d regions;"%(len(newGeom.regs)))
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
                            newGeom.addReg(Region.fromBuf(regBuf))
                            regBuf="" # flush region def buffer
                        regBuf=tmpBuf+tmpLine
                        tmpBuf=""
                    
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
        if (myTitle is None):
            myTitle="appended geometries"
        new.title=myTitle
        return new

    def echo(self,oFileName,lSplit=False,what="all",dMode="w"):
        '''
        - what="all"/"bodies"/"regions"/"materials"
        '''
        import os
        
        if (not oFileName.endswith(".inp")):
            oFileName=oFileName+".inp"

        if (what.upper()=="BODIES"):
            print("saving bodies in file %s..."%(oFileName))
            ff=open(oFileName,dMode)
            ff.write("* \n")
            for tmpBody in self.bods:
                ff.write("%s\n"%(tmpBody.echo()))
            ff.write("* \n")
            ff.close()
            print("...done;")
        elif (what.upper()=="REGIONS"):
            print("saving regions in file %s..."%(oFileName))
            ff=open(oFileName,dMode)
            ff.write("* \n")
            for tmpReg in self.regs:
                ff.write("%s\n"%(tmpReg.echo()))
            ff.write("* \n")
            ff.close()
            print("...done;")
        elif (what.upper()=="MATERIALS"):
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
                self.echo(oFileName.replace(".inp","_bodies.inp",1),\
                          lSplit=False,what="bodies",dMode="w")
                self.echo(oFileName.replace(".inp","_regions.inp",1),\
                          lSplit=False,what="regions",dMode="w")
                self.echo(oFileName.replace(".inp","_assignmats.inp",1),\
                          lSplit=False,what="materials",dMode="w")
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
            oldBodyNames.append(self.bods[iBody].bName)
            newBodyNames.append(newNameFmt%(iBody+1))
            self.bods[iBody].rename(newBodyNames[-1],lNotify=lNotify)
        for iReg in range(len(self.regs)):
            self.regs[iReg].rename(newNameFmt%(iReg+1),lNotify=lNotify)
            self.regs[iReg].BodyNameReplaceInDef(oldBodyNames,newBodyNames)

class MergeGeo(Geometry):
    '''
    A dedicated class for geometries that should be merged with
      others, e.g. clones arranged on a grid or the hive containing
      them
    '''
    def __init__(self):
        Geometry.__init__(self)
        self.rCont=[]  # containment indicator (1 per region):
                       # -1: region to be contained into/sized by another one
                       #  0: regular region (neither contains nor it is contained)
                       #  1: cell region (it contains another region)
        self.rCent=[]  # central point of one or more region (3D arrays)

    def addReg(self,tmpReg,rCont=0,rCent=[0.0,0.0,0.0]):
        '''
        overriding Geometry.addReg(self,tmpReg)
        '''
        self.regs.append(tmpReg)
        self.rCont.append(rCont)
        self.rCent.append(np.array(rCent))

    def flagRegs(self,whichRegs,rCont,rCent):
        if (whichRegs is str):
            if (whichRegs.upper()=="ALL"):
                regs2mod,iRegs2mod=myGeo.ret("region","ALL")
            else:
                regs2mod=[whichRegs]
        else:
            if ( "ALL" in [tmpStr.upper() for tmpStr in whichRegs] ):
                regs2mod,iRegs2mod=myGeo.ret("region","ALL")
            else:
                regs2mod=whichRegs
        for whichReg in whichRegs:
            outReg,iOutReg=self.ret("region",whichReg)
            self.rCont[iOutReg]=rCont
            self.rCent[iOutReg]=rCent

    @staticmethod
    def appendGeometries(myGeos,myTitle=None):
        '''
        overriding Geometry.appendGeometries(myGeos,myTitle)
        myGeos is a list of MergeGeo instances
        '''
        new=MergeGeo()
        for myGeo in myGeos:
            new.bods=new.bods+myGeo.bods
            new.regs=new.regs+myGeo.regs
            new.rCont=new.rCont+myGeo.rCont
            new.rCent=new.rCent+myGeo.rCent
        if (myTitle is None):
            myTitle="appended geometries"
        new.title=myTitle
        return new

    @staticmethod
    def ImportFromGeometry(inGeo):
        '''
        converting an existing Geometry instance into a MergeGeo instance
        '''
        outMergeGeo=MergeGeo()
        outMergeGeo.bods=deepcopy(inGeo.bods)
        outMergeGeo.title=deepcopy(inGeo.title)
        for tmpReg in inGeo.regs:
            outMergeGeo.addReg(deepcopy(tmpReg))
        return outMergeGeo

    @staticmethod
    def DefineHive_SphericalShell(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,tmpTitle="Hive for a spherical shell",lDebug=True):
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
        newGeom=MergeGeo()
        
        print("   ...bodies...")
        spheres=[]
        for ii,RR in enumerate(RRs,1):
            tmpBD=Body()
            tmpBD.bName="HVRAD%03i"%(ii)
            tmpBD.bType="SPH"
            tmpBD.Rs[0]=RR
            tmpBD.comment="* Hive radial boundary at R[cm]=%g"%(RR)
            spheres.append(tmpBD)
        spheres[0].comment="* \n"+spheres[0].comment
        thetas=[]
        for ii,TT in enumerate(TTs,1):
            tmpBD=Body()
            tmpBD.bName="HVTHT%03i"%(ii)
            tmpBD.V=np.array([0.0,1.0,0.0])
            tmpBD.rotate(myMat=None,myTheta=-TT,myAxis=1,lDegs=True,lDebug=lDebug)
            tmpBD.comment="* Hive theta boundary at theta[deg]=%g"%(TT)
            thetas.append(tmpBD)
        thetas[0].comment="* \n"+thetas[0].comment
        phis=[]
        for ii,PP in enumerate(PPs,1):
            tmpBD=Body()
            tmpBD.bName="HVPHI%03i"%(ii)
            tmpBD.V=np.array([1.0,0.0,0.0])
            tmpBD.rotate(myMat=None,myTheta=PP,myAxis=2,lDegs=True,lDebug=lDebug)
            tmpBD.comment="* Hive phi boundary at phi[deg]=%g"%(PP)
            phis.append(tmpBD)
        phis[0].comment="* \n"+phis[0].comment
        newGeom.bods=spheres+thetas+phis

        print("   ...regions...")
        # - outside hive
        tmpReg=Region()
        tmpReg.rName="HV_OUTER"
        tmpReg.material="VACUUM"
        tmpReg.definition=''' | +%-8s | -%-8s
                 | +%-8s -%-8s -%-8s
                 | +%-8s -%-8s +%-8s
                 | +%-8s -%-8s +%-8s -%-8s +%-8s
                 | +%-8s -%-8s +%-8s -%-8s -%-8s
'''%(spheres[0].bName,spheres[-1].bName, \
     spheres[-1].bName,spheres[0].bName, thetas[-1].bName, \
     spheres[-1].bName,spheres[0].bName, thetas[ 0].bName, \
     spheres[-1].bName,spheres[0].bName, thetas[-1].bName, thetas[ 0].bName, phis[ 0].bName, \
     spheres[-1].bName,spheres[0].bName, thetas[-1].bName, thetas[ 0].bName, phis[-1].bName  )
        tmpReg.comment="* region outside hive"
        newGeom.addReg(tmpReg)
        # - inside hive
        iHive=0
        for iR in range(1,len(spheres)):
            for iT in range(1,len(thetas)):
                for iP in range(1,len(phis)):
                    tmpReg=Region()
                    tmpReg.rName="HVCL%04i"%(iHive)
                    tmpReg.material="VACUUM"
                    tmpReg.definition='+%-8s -%-8s +%-8s -%-8s +%-8s -%-8s'%\
                        (spheres[iR].bName,spheres[iR-1].bName,\
                         thetas [iT].bName,thetas [iT-1].bName,\
                         phis   [iP].bName,phis   [iP-1].bName)
                    tmpComment=             "* - hive region %4d: R[cm]=[%g:%g], theta[deg]=[%g:%g], phi[deg]=[%g:%g]"%(
                        iHive,RRs[iR-1],RRs[iR],TTs[iT-1],TTs[iT],PPs[iP-1],PPs[iP])
                    myCenter=cellGrid.ret(myWhat="POINT",iEl=iHive)
                    tmpComment=tmpComment+"\n*   center=[%g,%g,%g];"%(myCenter[0],myCenter[1],myCenter[2])
                    tmpReg.comment=tmpComment
                    newGeom.addReg(tmpReg,rCont=1,rCent=myCenter)
                    iHive=iHive+1

        newGeom.setTitle(tmpTitle=tmpTitle)
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
                print("MergeGeo.BuildGriddedGeo_SphericalShell(): unknown prototype %s!"%(\
                        myProtoList[iLoc]))
                exit(1)
            # - clone prototype
            myGeo=MergeGeo.ImportFromGeometry(myProtoGeos[myProtoList[iLoc]])
            # - move clone to requested location/orientation
            myGeo.solidTrasform(dd=myLoc.ret("POINT"),myMat=myLoc.ret("MATRIX"),lDebug=lDebug)
            # - flag the region(s) outside the prototypes or that should be sized
            #   by the hive cells
            myGeo.flagRegs(osRegNames,-1,myLoc.ret("POINT"))
            # - rename the clone
            baseName="GR%03d"%(iLoc)
            myGeo.rename(baseName)
            # - notify the user about original prototype and location
            myGeo.headMe("* \n"+ \
                         "* "+"="*108+"\n"+ \
                         "* GRID cell # %3d - family name: %s - prototype: %s\n"%(\
                            iLoc,baseName,myProtoList[iLoc])+ \
                         "* "+myLoc.echo(myFmt="% 13.6E",mySep="\n* ") + \
                         "-"*108 )
            # - append clone to list of geometries
            myGeos.append(myGeo)
        # return merged geometry
        return MergeGeo.appendGeometries(myGeos)

    @staticmethod
    def MergeGeos(hiveGeo,gridGeo):
        '''
        This method merges one FLUKA geometry onto another one.

        input parameters:
        - hiveGeo: MergeGeo instance of the hive;
        - gridGeo: MergeGeo instance of the grid of objects.

        The two geometries must not have common names - no check is performed
          for the time being.

        All the regions of gridGeo with rCont==-1 will be matched with regions
          of hiveGeo with rCont==1; a one-to-one mapping is established based
          on the rCent arrays. The merged regions will still belong to gridGeo
          and the respective region in hiveGeo will disappear.

        Similarly, all regions of gridGeo with rCont==-2 will be matched with regions
          of hiveGeo with rCont==1; a one-to-one mapping is established based
          on the rCent arrays. The merged regions will still belong to gridGeo
          and the respective region in hiveGeo will disappear.
        '''
        myGeos=[]
        # return merged geometry
        return MergeGeo.appendGeometries(myGeos)
    
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
    R=500
    dR=50
    NR=1
    Tmax=3    # theta [degs] --> range: -Tmax:Tmax
    NT=2      # number of steps (i.e. entities)
    Pmax=2    # phi [degs] --> range: -Pmax:Pmax
    NP=2      # number of steps (i.e. entities)

    # - hive geometry
    HiveGeo=MergeGeo.DefineHive_SphericalShell(R,R+dR,NR,-Tmax,Tmax,NT,-Pmax,Pmax,NP,lDebug=lDebug)
    HiveGeo.echo("hive.inp")

    # - gridded crystals
    #   acquire geometries
    fileNames=[ "caloCrys.inp" ] ; geoNames=fileNames
    myProtoGeos=acquireGeometries(fileNames,geoNames=geoNames);
    cellGrid=grid.Grid.SphericalShell(R,R+dR,NR,-Tmax,Tmax,NT,-Pmax,Pmax,NP,lDebug=lDebug)
    myProtoList=[ "caloCrys.inp" for ii in range(len(cellGrid)) ]
    GridGeo=MergeGeo.BuildGriddedGeo(cellGrid,myProtoList,myProtoGeos,osRegNames=["OUTER"],lDebug=lDebug)
    GridGeo.echo("grid.inp")
