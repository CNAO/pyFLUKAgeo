# classes for managing FLUKA bodies
# A. Mereghetti, 2024/04/16
# python version: >= 3.8.10

import numpy as np
from myMath import RotMat
from FLUKA import GeoObject, echoFloats, assembleLine, MaxLineLength, LineHeader

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

    def echo(self,lFree=True,maxLen=MaxLineLength,header=LineHeader,lMultiLine=True):
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
        elif (self.bType=="RPP"):
            myFloats=[0.0 for ii in range(6)]
            myFloats[::2]=list(self.P)
            myFloats[1::2]=list(self.P+self.V)
        else:
            print("body %s NOT supported yet!"%(self.bType))
            exit(1)
        myStrings=[myStr]+echoFloats(myFloats,lFree=lFree)
        myStr=assembleLine(myStrings,maxLen=maxLen,header=header,lMultiLine=lMultiLine)
        if (self.isLinkedToTransform()):
            myStr="$Start_transform -%s\n%s\n$end_transform"%(self.TransformName,myStr)
        return GeoObject.echoComm(self)+myStr

    def isRotatable(self):
        return self.lIsRotatable

    def linkTransformName(self,Tname=None):
        if (Tname is not None or Tname!=""):
            self.TransformName=Tname
            
    def retTransformName(self):
        return self.TransformName

    def retCenter(self,myType=0):
        '''
        myType: specifies which center to return:
        actual (==0), upstream (<0), downstream (>0)
        '''
        if (self.bType=="PLA" or self.bType=="SPH" or \
            self.bType=="YZP" or self.bType=="XZP" or self.bType=="XYP"):
            return self.P
        elif (self.bType=="TRC" or self.bType=="RCC" or
              self.bType=="XCC" or self.bType=="YCC" or self.bType=="ZCC"):
            if (myType==0):
                return self.P+0.5*self.V
            elif (myType<0):
                return self.P
            elif (myType>0):
                return self.P+self.V
        elif (self.bType=="RPP"):
            myOrient=self.retOrient()
            if (myType==0):
                return self.P+0.5*self.V
            elif (myType<0):
                return self.P+np.multiply(0.5*self.V,1-myOrient)
            elif (myType>0):
                return self.P+np.multiply(0.5*self.V,1-myOrient)+\
                              np.multiply(    self.V,  myOrient)
        
    def retOrient(self):
        if (self.bType=="RPP"):
            if (self.V[0]==self.V[1] and self.V[2]==self.V[1]):
                # actually a cube: orientation as z-axis by default
                myV=np.array([0.0,0.0,1.0])
            elif (self.V[0]==self.V[1] and self.V[2]!=self.V[0]):
                # xy dims are identical: orientation as z-axis
                myV=np.array([0.0,0.0,1.0])
            elif (self.V[2]==self.V[0] and self.V[1]!=self.V[2]):
                # xz dims are identical: orientation as y-axis
                myV=np.array([0.0,1.0,0.0])
            elif (self.V[1]==self.V[2] and self.V[0]!=self.V[1]):
                # yz dims are identical: orientation as x-axis
                myV=np.array([1.0,0.0,0.0])
            else:
                # orientation is where the largest delta is found
                myV=np.zeros(3)
                myV[np.argmax(self.V)]=1.0
        else:
            myV=self.V/np.linalg.norm(self.V)
        return myV

    def retL(self):
        LL=0.0
        if (self.bType=="RCC" or self.bType=="REC"):
            LL=np.linalg.norm(self.V)
        elif (self.bType=="RPP"):
            LL=np.linalg.norm(np.multiply(self.retOrient(),self.V))
        return LL

    def isLinkedToTransform(self):
        return self.TransformName is not None

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
            elif (data[0]=="RPP"):
                newBody.bType=data[0]
                newBody.rename(data[1],lNotify=False)
                newBody.P=np.array(data[2:7:2]).astype(float)
                newBody.V=np.array(data[3:8:2]).astype(float)-newBody.P
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
               self.bType="PLA"
               self.V=np.array([1.0,0.0,0.0])
           elif (self.bType=="XZP"):
               self.bType="PLA"
               self.V=np.array([0.0,1.0,0.0])
           elif (self.bType=="XYP"):
               self.bType="PLA"
               self.V=np.array([0.0,0.0,1.0])
           elif (self.bType=="XCC"):
               self.bType="RCC"
               self.P[0]=-infL
               self.V=np.array([2*infL,0.0,0.0])
           elif (self.bType=="YCC"):
               self.bType="RCC"
               self.P[1]=-infL
               self.V=np.array([0.0,2*infL,0.0])
           elif (self.bType=="ZCC"):
               self.bType="RCC"
               self.P[2]=-infL
               self.V=np.array([0.0,0.0,2*infL])
           else:
               print("cannot make body %s rotatable!"%(myStr))
               exit(1)
           if (lDebug):
               print("...converted body %s into an %s"%(myStr,self.bType))
           self.lIsRotatable=True

    def makeUNrotatable(self,lDebug=False,infL=1000.0):
        if (self.isRotatable()):
            myStr=GeoObject.echoComm(self)+"%-3s %8s"%(self.bType,self.echoName())
            myV=self.retOrient()
            if (self.bType=="PLA"):
                if (np.allclose(myV,[1.0,0.0,0.0])):
                    self.bType="YZP"
                elif (np.allclose(myV,[0.0,1.0,0.0])):
                    self.bType="XZP"
                elif (np.allclose(myV,[0.0,0.0,1.0])):
                    self.bType="XYP"
            elif (self.bType=="RCC" and self.retL()>infL):
                if (np.allclose(myV,[1.0,0.0,0.0])):
                    self.bType="XCC"
                elif (np.allclose(myV,[0.0,1.0,0.0])):
                    self.bType="YCC"
                elif (np.allclose(myV,[0.0,0.0,1.0])):
                    self.bType="ZCC"
            else: return
            self.lIsRotatable=False
            if (lDebug):
                print("...converted body %s into an %s"%(myStr,self.bType))

    def traslate(self,dd=None,lDebug=False):
        if (dd is not None):
            if (lDebug):
                myStr="%-3s %8s"%(self.bType,self.echoName())
                print("applying traslation array [%g,%g,%g] to body %s..."%(dd[0],dd[1],dd[2],myStr))
            self.P=self.P+dd
        else:
            print("please provide me with a 3D traslation array!")
            exit(1)
                
    def rotate(self,myMat=None,myTheta=None,myAxis=3,lDegs=True,lDebug=False):
        if (myMat is not None):
            if (lDebug):
                myStr="%-3s %8s"%(self.bType,self.echoName())
                print("applying rotation matrix to body %s..."%(myStr))
            self.P=myMat.mulArr(self.P,lDebug=lDebug)
            self.V=myMat.mulArr(self.V,lDebug=lDebug)
            # DIRTY PATCH TO ROTATE RPPs by 90degs and multiples:
            if (self.bType=="RPP"):
                for ii in range(3):
                    if (self.V[ii]<0.0):
                        self.P[ii]=self.P[ii]+self.V[ii]
                        self.V[ii]=-self.V[ii]
        elif (myTheta is not None):
            if (myAxis is None):
                print("you must specify a rotation axis together with an angle!")
                exit(1)
            else:
                if (lDebug):
                    myStr="%-3s %8s"%(self.bType,self.echoName())
                    print("applying rotation by %g degs on axis %d to body %s..."%(myTheta,myAxis,myStr))
                myMat=RotMat(myAng=myTheta,myAxis=myAxis,lDegs=lDegs,lDebug=lDebug)
                self.rotate(myMat=myMat,myTheta=None,myAxis=myAxis,lDegs=lDegs,lDebug=lDebug)

    def resize(self,newL,lDebug=False):
        if (self.bType=="RCC" or self.bType=="RPP"):
           origCenter=self.retCenter()
           orient=self.retOrient()
           if (lDebug):
               myStr="%-3s %8s"%(self.bType,self.echoName())
               print("...resizing body %s to L=%g cm..."%(myStr,newL))
               oldP=self.P
               oldV=self.V
           if (self.bType=="RCC"):
               self.V=orient*newL
           elif (self.bType=="RPP"):
               self.V=self.V+(newL-np.dot(orient,self.V))*orient
           self.P=origCenter-self.V/2.
           if (lDebug):
               print("   ...center:",origCenter,"-->",self.retCenter())
               print("   ...orientation:",orient,"-->",self.retOrient())
               print("   ...self.P:",oldP,"-->",self.P)
               print("   ...self.V:",oldV,"-->",self.V)
                
