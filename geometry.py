# classes for managing FLUKA geometries
# A. Mereghetti, 2023/09/13
# python version: >= 3.8.10;
#
# - name-based FLUKA geometry defition;
#
# Class body:
# - only planes and spheres, for the time being;

import numpy as np

class body:

    def __init__(self):
        '''
        default body is a FLUKA PLA
        '''
        self.bName=""
        self.bType="PLA"
        self.P=np.zeros(3)
        self.V=np.array([0.0 0.0 1.0])
        self.Rs=np.zeros(2)

    def echo(self):
        if (self.bType=="PLA"):
            return "PLA %8s % 15.8E % 15.8E % 15.8E % 15.8E % 15.8E % 15.8E" \
             % ( self.bName, self.P[0], self.P[1], self.P[2], self.V[0], self.V[1], self.V[2] )
        elif (self.bType=="SPH"):
            return "SPH %8s % 15.8E % 15.8E % 15.8E % 15.8E" \
             % ( self.bName, self.P[0], self.P[1], self.P[2], self.Rs[0] )
        else:
            print("body %s NOT supported yet!"%(self.bType))
            exit(1)

    @staticmethod
    def fromLine(tmpLine):
        newBody=body()
        data=tmpLine.split()
        if (data[0]=="PLA"):
            newBody.bType="PLA"
            newBody.bName=data[1]
            newBody.V=np.array(data[2:5])
            newBody.P=np.array(data[5:8])
        elif (data[0]=="YZP"):
            newBody.bType="PLA"
            newBody.bName=data[1]
            newBody.V=np.array([1.0 0.0 0.0])
            newBody.P[0]=data[2]
        elif (data[0]=="XZP"):
            newBody.bType="PLA"
            newBody.bName=data[1]
            newBody.V=np.array([0.0 1.0 0.0])
            newBody.P[1]=data[2]
        elif (data[0]=="XYP"):
            newBody.bType="PLA"
            newBody.bName=data[1]
            newBody.V=np.array([0.0 0.0 1.0])
            newBody.P[2]=data[2]
        elif (data[0]=="SPH"):
            newBody.bType="SPH"
            newBody.bName=data[1]
            newBody.P=np.array(data[2:5])
            newBody.Rs[0]=data[5]
        else:
            print("body %s NOT supported yet!"%(self.bType))
            exit(1)
        return newBody

    def traslate(self,dd=None):
        if (dd is not None):
            if (self.bType=="PLA" or self.bType=="SPH" ):
                self.P=self.P+dd
        else:
            print("body %s NOT supported yet!"%(self.bType))
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
        self.bName=newName
                
class region():

    def __init__(self):
        self.rName=""
        self.neigh=5
        self.definition=""
        self.material="BLACKHOLE"

    @staticmethod
    def fromDefLines(myBuffer):
        newReg=region()
        data=myBuffer.split()
        newReg.rName=data[0]
        newReg.neigh=float(data[1])
        newReg.definition=myBuffer
        # remove region name and number of neighbour regions from definition
        newReg.definition=newReg.definition.replace(data[0],"",1)
        newReg.definition=newReg.definition.replace(data[1],"",1)
        return newReg
        
    def stringReplace(self,oldStrings,newStrings):
        for oString,nString in zip(oldStrings,newStrings):
            if (oString==self.rName):
                self.rName=oString
            self.definition.replace(oString,nString)

    def assignMat(self,myMaterial):
        self.material=myMaterial
