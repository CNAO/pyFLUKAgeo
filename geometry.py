# classes for managing FLUKA geometries
# A. Mereghetti, 2023/09/13
# python version: >= 3.8.10;

import numpy as np

class Body:
    '''
    - supported bodies: all planes and spheres, nothing else for the time being;
    - body definition on ONE line only, and always starts at col 1;
    - support only for preceding comments, NO comments headed by !
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
            return myBuf+"PLA %8s % 15.8E % 15.8E % 15.8E % 15.8E % 15.8E % 15.8E" \
             % ( self.bName, self.P[0], self.P[1], self.P[2], self.V[0], self.V[1], self.V[2] )
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
                newBody.V=np.array(data[2:5])
                newBody.P=np.array(data[5:8])
            elif (data[0]=="YZP"):
                newBody.bType="PLA"
                newBody.bName=data[1]
                newBody.V=np.array([1.0,0.0,0.0])
                newBody.P[0]=data[2]
            elif (data[0]=="XZP"):
                newBody.bType="PLA"
                newBody.bName=data[1]
                newBody.V=np.array([0.0,1.0,0.0])
                newBody.P[1]=data[2]
            elif (data[0]=="XYP"):
                newBody.bType="PLA"
                newBody.bName=data[1]
                newBody.V=np.array([0.0,0.0,1.0])
                newBody.P[2]=data[2]
            elif (data[0]=="SPH"):
                newBody.bType="SPH"
                newBody.bName=data[1]
                newBody.P=np.array(data[2:5])
                newBody.Rs[0]=data[5]
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
                
    def rotate(self,myMat=None,myTheta=None,myAxis=3):
        if (myMat is not None):
            print("rotation by matrix NOT supported yet!")
            exit(1)
        elif (myTheta is not None):
            if (myAxis is None):
                print("you must specify a rotation axis together with an angle!")
                exit(1)
            else:
                print("rotation by axis-angle NOT supported yet!")
                exit(1)
                
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
    - support only for preceding comments or commented lines between region definition lines,
      NO comments headed by !
    - a region definition always starts at column 1;
    - ASSIGNMA cards: only 1:1 material:region assignment, NO material assignment
      for decay radiation simulation, NO magnetic fields;
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
                    # remove region name and number of neighbour regions from definition
                    newReg.definition=newReg.definition.replace(data[0],"",1)
                    newReg.definition=newReg.definition.replace(data[1],"",1)
                    newReg.definition=newReg.definition.strip() # remove heading/trailing empty spaces
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
                self.comment=self.comment+"* NAME CHANGE: FROM %s TO %s"%(self.rName,nString)
                self.rName=nString
            self.definition.replace(oString,nString)

    def assignMat(self,myMaterial):
        self.material=myMaterial
        
    def echo(self):
        # take into account comment lines
        myBuf=""
        if (len(self.comment)>0):
            myBuf=self.comment+"\n"
        return myBuf+"%8s %4d %s" % ( self.rName, self.neigh, self.definition )

    def echoAssignMat(self):
        return "ASSIGNMA  %10s%10s" % ( self.material, self.rName )

class Geometry():
    '''
    - name-based FLUKA geometry defition;
    - NO support of LATTICE cards;
    - NO support for #input cards or geo defitions in files other than the one being parsed;
    - comments:
      . body: commented lines are considered always before the comment, and only
              commented lines before the body will be retained;
              --> trailing comments to body declaration section will disappear!
      . region: commented lines are kept where they are found, if they are found before
                or along the region declaration;
              --> trailing comments to region declaration section will disappear!
    '''
    def __init__(self):
        self.bods=[]
        self.regs=[]
        self.title=""

    def addBod(self,tmpBod):
        self.bods.append(tmpBod)
 
    def addReg(self,tmpReg):
        self.regs.append(tmpReg)

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
                        regBuf="" # flush buffer
                    print("...acquired %d regions;"%(len(newGeom.regs)))
                    lRead=0 # ignore LATTICE cards
                else:
                    if (tmpLine.startswith("*")):
                        # comment line: store in buffer
                        tmpBuf=tmpBuf+tmpLine
                        continue
                    if (tmpLine.startswith(" ")):
                        regBuf=regBuf+tmpBuf+tmpLine
                        tmpBuf="" # flush buffer
                    else:
                        if (len(regBuf)>0):
                            # acquire region
                            newGeom.addReg(Region.fromBuf(regBuf))
                        tmpBuf=tmpBuf+tmpLine
                        regBuf=tmpBuf
                    
        ff.close()
        print("...done;")
 
if (__name__=="__main__"):
    caloCrysGeo=Geometry.fromInp("caloCrys.inp")
