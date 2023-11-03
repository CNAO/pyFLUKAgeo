# classes for managing FLUKA geometries
# A. Mereghetti, 2023/09/13
# python version: >= 3.8.10;

import numpy as np
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
                
    def rotate(self,myMat=None,myTheta=None,myAxis=3,lDegs=True,lDebug=True):
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
                
    def rename(self,newName):
        # notify user as a comment line
        if (len(self.comment)>0):
            # ...trailing to any existing comment
            self.comment=self.comment+"\n"
        self.comment=self.comment+"* NAME CHANGE: FROM %s TO %s"%(self.bName,newName)
        self.bName=newName
                
class Region():
    '''
    - no parsing/splitting of zones;
    - support only for preceding comments or commented lines between region
      definition lines, NO comments headed by !
    - a region definition always starts at column 1;
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
        
    def stringReplace(self,oldStrings,newStrings):
        for oString,nString in zip(oldStrings,newStrings):
            if (oString==self.rName):
                # notify user as a comment line
                if (len(self.comment)>0):
                    # ...trailing to any existing comment
                    self.comment=self.comment+"\n"
                self.comment=self.comment+\
                    "* NAME CHANGE: FROM %s TO %s"%(self.rName,nString)
                self.rName=nString
            self.definition.replace(oString,nString)

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
        if (len(self.bods[0].comment)>0):
            self.bods[0].comment=mySting+"\n"+self.bods[0].comment
        else:
            self.bods[0].comment=mySting
        if (len(self.regs[0].comment)>0):
            self.regs[0].comment=mySting+"\n"+self.regs[0].comment
        else:
            self.regs[0].comment=mySting
            
    def ret(self,myWhat,myName):
        lFound=False
        if (myWhat.upper()=="BODY"):
            for iEntry,myEntry in enumerate(self.bods):
                if (myEntry.bName==myName):
                    lFound=True
                    break
        elif (myWhat.upper()=="REGION"):
            for iEntry,myEntry in enumerate(self.bods):
                if (myEntry.bName==myName):
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
    def DefineHive_SphericalShell(RRs,TTs,PPs,tmpTitle="Hive for a spherical shell"):
        '''
        This method defines a the hive for a grid on a spherical shell.

        The grid is described in spherical coordinates:
        - r [cm];
        - theta [degs]: angle in yz-plane (positive when pointing towards y>0);
        - phi [degs]: angle in xz-plane (positive when pointing towards x>0).

        The grid is centred around the z-axis and symmetric in the two
          directions, up/down (y>0/y<0), left/right (x>0/x<0).

        The hive is delimited:
        - radially, by spheres;
        - on theta, by rotated XZPs;
        - on phi, by rotated YZPs;

        RRs, TTs, and PPs must be in increasing order
        '''
        print("Preparing the hive for a spherical shell...")
        
        print("...check of input info...")
        for tmpR in RRs:
            if (tmpR<0.0):
                print("Geometry.DefineHive_SphericalShell(): Error negative R: %g<0!"%(tmpR))
                exit(1)
        if ((TTs[-1]-TTs[0])>180):
            print("Geometry.DefineHive_SphericalShell(): Theta range too large: %g>180!"%(TTs[-1]-TTs[0]))
            exit(1)
        if ((PPs[-1]-PPs[0])>180):
            print("Geometry.DefineHive_SphericalShell(): Theta range too large: %g>180!"%(PPs[-1]-PPs[0]))
            exit(1)
                
        newGeom=Geometry()

        print("...generating bodies...")
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
            tmpBD.rotate(myMat=None,myTheta=-TT,myAxis=1,lDegs=True)
            tmpBD.comment="* Hive theta boundary at theta[deg]=%g"%(TT)
            thetas.append(tmpBD)
        thetas[0].comment="* \n"+thetas[0].comment
        phis=[]
        for ii,PP in enumerate(PPs,1):
            tmpBD=Body()
            tmpBD.bName="HVPHI%03i"%(ii)
            tmpBD.V=np.array([1.0,0.0,0.0])
            tmpBD.rotate(myMat=None,myTheta=PP,myAxis=2,lDegs=True)
            tmpBD.comment="* Hive phi boundary at phi[deg]=%g"%(PP)
            phis.append(tmpBD)
        phis[0].comment="* \n"+phis[0].comment
        newGeom.bods=spheres+thetas+phis

        print("...generating regions...")
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
                    iHive=iHive+1
                    tmpReg=Region()
                    tmpReg.rName="HVCL%04i"%(iHive)
                    tmpReg.material="VACUUM"
                    tmpReg.definition='+%-8s -%-8s +%-8s -%-8s +%-8s -%-8s'%\
                        (spheres[iR].bName,spheres[iR-1].bName,\
                         thetas [iT].bName,thetas [iT-1].bName,\
                         phis   [iP].bName,phis   [iP-1].bName)
                    tmpReg.comment="* - hive region: R[cm]=[%g:%g], theta[deg]=[%g:%g], phi[deg]=[%g:%g]"%(RRs[iR],RRs[iR-1],TTs[iR],TTs[iR-1],PPs[iR],PPs[iR-1])
                    newGeom.addReg(tmpReg)

        newGeom.setTitle(tmpTitle=tmpTitle)
        return newGeom

    @staticmethod
    def BuildGriddedGeo_SphericalShell(myGrid,myProtoList,myProtoGeos):
        '''
        This method defines a FLUKA geometry representing a grid of objects
           distributed on a spherical shell.

        input parameters:
        - grid: an instance of a GRID() class, i.e. a list of locations;
        - myProtoList: this list states which prototype should be used at each
                       location. NB: len(myProtoList)=len(myGrid);
        - myProtoGeos: dictionary of actual prototype geometries. The unique
                       entries of myProtoList are the full set or a subset of
                       the keys of this dictionary;
        '''
        myGeos=[]
        return Geometry.mergeGeometries(myGeos)
    
    @staticmethod
    def mergeGeometries(myGeos):
        '''
        Barely appending geometries one to another;
        '''
        new=Geometry()
        for myGeo in myGeos:
            new.bods=new.bods+myGeo.bods
            new.regs=new.bods+myGeo.regs
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

    def solidTrasform(self,dd=None,myMat=None,myTheta=None,myAxis=3,lDegs=True,lDebug=True):
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
            
def acquireGeometries(fileNames,geoNames=None):
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
    lDebug=True
    # # - manipulate a geometry
    # caloCrysGeo=Geometry.fromInp("caloCrys.inp")
    # myMat=RotMat(myAng=60,myAxis=3,lDegs=True,lDebug=lDebug)
    # caloCrysGeo.solidTrasform(dd=[0,10,-20],myMat=myMat)
    # caloCrysGeo.echo("pippo.inp")
    # - acquire geometries
    fileNames=[ "caloCrys.inp","caloCrys.inp" ]
    myGeometries=acquireGeometries(fileNames,geoNames=["pippo","pluto"]);
    
    # # - test generation of hive
    # R=500
    # dR=50
    # Tmax=30   # theta [degs] --> range: -Tmax:Tmax
    # NT=20     # number of steps (i.e. entities)
    # Pmax=20   # phi [degs] --> range: -Pmax:Pmax
    # NP=20     # number of steps (i.e. entities)
    # RRs,TTs,PPs=grid.DefHiveBoundaries_SphericalShell_OneLayer(R,dR,Tmax,NT,Pmax,NP)
    # HiveGeo=Geometry.DefineHive_SphericalShell(RRs,TTs,PPs,tmpTitle="Hive for a single-layer spherical shell")
    # HiveGeo.echo("pippo.inp")
