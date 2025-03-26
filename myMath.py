# classes for managing math aspects
# A. Mereghetti, 2023/09/29
# python version: >= 3.8.10;

import numpy as np
from copy import deepcopy

class Matrix:
    '''
    very simple class coding a 2D SQUARED matrix
    '''

    def __init__(self,nDim=3):
        self.vals=np.zeros((nDim,nDim))
        self.nDim=nDim

    def __getitem__(self,pos):
        ii,jj=pos
        return self.vals[ii][jj]
    
    def __setitem__(self,pos,value):
        ii,jj=pos
        self.vals[ii][jj]=value

    def echo(self,myFmt="% .12E"):
        buf=""
        for ii in range(self.nDim):
            buf=buf+"|"
            for jj in range(self.nDim):
                buf=buf+" "+myFmt%(self[ii,jj])
            buf=buf+" |\n"
        return buf

    def mulMat(self,RMat,lDebug=True):
        '''
        method to multiply squared matrices of the same number of rows:
                  R_i,j = SUM_k self_i,k x RMat_k,j
        '''
        if (self.nDim!=RMat.nDim):
            print("cannot multiply a %dx%d matrix times a %dx%d matrix"%\
                  (self.nDim,self.nDim,RMat.nDim,RMat.nDim))
            exit(1)
        newMat=deepcopy(self) # preserve original class in daughter classes
        for ii in range(self.nDim):
            for jj in range(self.nDim):
                newMat[ii,jj]=sum([ self[ii,kk]*RMat[kk,jj]\
                                        for kk in range(self.nDim)])
        if (lDebug):
            print("Matrix.mulMat():")
            print(newMat.echo())
        return newMat

    def mulArr(self,myArr,lDebug=True):
        'method to multiply a squared matrix by an array'
        if (self.nDim!=len(myArr)):
            print("cannot multiply a %dx%d matrix times a %d array"%\
                  (self.nDim,self.nDim,len(myArr)))
            exit(1)
        out=np.zeros((len(myArr)))
        for ii in range(self.nDim):
            out[ii]=sum([ self[ii,kk]*myArr[kk] \
                              for kk in range(self.nDim)])
        if (lDebug):
            print("Matrix.mulArr():",myArr,"-->",out)
        return out

    def mulSca(self,mySca):
        'method to multiply a squared matrix by an array'
        newMat=deepcopy(self) # preserve original class in daughter classes
        for ii in range(self.nDim):
            for jj in range(self.nDim):
                newMat[ii,jj]=self[ii,jj]*mySca
        return newMat

    def Minor(self,ll,mm):
        'calculate minor of self at pos ll,mm'
        newMat=Matrix(self.nDim-1)
        kk=0
        for ii in range(self.nDim):
            if (ii==ll):
                continue
            nn=0
            for jj in range(self.nDim):
                if (jj==mm):
                    continue
                newMat[kk,nn]=self[ii,jj]
                nn=nn+1
            kk=kk+1
        return newMat

    def Cofactor(self,ii,jj):
        if ((ii+jj)%2==1):
            return -1
        else:
            return  1

    def MinMat(self):
        'calculate matrix of minors'
        newMat=deepcopy(self) # preserve original class in daughter classes
        for ii in range(self.nDim):
            for jj in range(self.nDim):
                newMat[ii,jj]=self.Minor(ii,jj).det()
        return newMat

    def CofactMat(self):
        'calculate matrix of co-factors'
        newMat=deepcopy(self) # preserve original class in daughter classes
        for ii in range(self.nDim):
            for jj in range(self.nDim):
                if ((ii+jj)%2==1):
                    newMat[ii,jj]=newMat[ii,jj]*self.Cofactor(ii,jj)
        return newMat

    def AdjugateMat(self):
        'calculate adjugate matrix, i.e. the transposed'
        newMat=deepcopy(self) # preserve original class in daughter classes
        for ii in range(self.nDim):
            for jj in range(self.nDim):
                newMat[ii,jj]=self[jj,ii]
        return newMat

    def det(self):
        'calculate the determinant of the matrix'
        if (self.nDim==2):
            return self[0,0]*self[1,1]-self[1,0]*self[0,1]
        else:
            det=0.0
            for ii in range(self.nDim):
                det=det+self[0,ii]*self.Minor(0,ii).det()*self.Cofactor(0,ii)
            return det

    def inv(self):
        '''
        invert the matrix - based on:
        https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html
        '''
        newMat=deepcopy(self) # preserve original class in daughter classes
        newMat=newMat.MinMat()
        newMat=newMat.CofactMat()
        newMat=newMat.AdjugateMat()
        return newMat.mulSca(1/self.det())

class UnitMat(Matrix):
    '''
    very simple class coding a 2D unit matrix
    '''
    def __init__(self,nDim=3):
        Matrix.__init__(self,nDim=nDim)
        for ii in range(self.nDim):
            self[ii,ii]=1.0

class RotMat(UnitMat):
    '''
    very simple class coding 3D rotations

    Reminders:
    Rx = |  1  0  0 |    Ry = |  C  0  S |    Rz = |  C -S  0 |
         |  0  C -S |         |  0  1  0 |         |  S  C  0 |
         |  0  S  C |         | -S  0  C |         |  0  0  1 |
    where:
    . C=cos(theta)
    . S=sin(theta)
    . theta>0: anti-clockwise rotation when seen from the rotation axis
    '''
    def __init__(self,myAng=0.0,myAxis=3,lDegs=True,lDebug=False):
        # first, initialise matrix as unit matrix;
        UnitMat.__init__(self,nDim=3)
        # then, fill it with the actual trigonometric values
        if (myAxis<1 or myAxis>3):
            print("wrong indication of axis: %d (either 1=x,2=y,3=z)"%(myAxis))
            exit(1)
        tAng=myAng
        if (lDegs):
            tAng=np.radians(tAng)
        if (tAng!=0.0):
            if (myAxis==1):
                self[1,1]= np.cos(tAng)
                self[2,2]= np.cos(tAng)
                self[1,2]=-np.sin(tAng)
                self[2,1]= np.sin(tAng)
            elif (myAxis==2):
                self[0,0]= np.cos(tAng)
                self[2,2]= np.cos(tAng)
                self[2,0]=-np.sin(tAng)
                self[0,2]= np.sin(tAng)
            elif (myAxis==3):
                self[0,0]= np.cos(tAng)
                self[1,1]= np.cos(tAng)
                self[0,1]=-np.sin(tAng)
                self[1,0]= np.sin(tAng)
        if (lDebug):
            print("RotMat.__init__():")
            print(self.echo())

    @staticmethod
    def ConcatenatedRotMatrices(myAngs=[],myAxes=[],lDegs=True,lDebug=True):
        # go through the array on angles and axes to define the overall
        #   matrix transformation
        rotMatrices=[]
        if (lDebug):
            print("creating rotation matrices...")
        for myAng,myAx in zip(myAngs,myAxes):
            angDegs=myAng
            angRads=myAng
            if (lDegs):
                angRads=np.radians(angRads)
            else:
                angDegs=np.degrees(angDegs)
            if (lDebug):
                print("...creating rotation matrix: %g degs around %d axis..."%\
                      (angDegs,myAx))
            rotMatrices.append(RotMat(myAng=myAng,myAxis=myAx,\
                          lDegs=lDegs,lDebug=lDebug))
        if (lDebug):
            print("concatenating them...")
        newMat=RotMat(lDebug=lDebug)
        for rotMat in rotMatrices:
            newMat=rotMat.mulMat(newMat,lDebug=lDebug)
            
        return newMat

    def GetGimbalAngles(self,lDegs=True):
        '''
        The current matrix is analysed as:
           R=RxRyRz=
             |    CyCz        -CySz      Sy  |
             | CxSz+SxSyCz CxCz-SxSySz -SxCy |
             | SxSz-CxSyCz SxCz+CxSySz  CxCy |
        Ranges:
        - thetas[0]=ang_x: [-pi:pi];
        - thetas[1]=ang_y: [-pi/2:pi/2];
        - thetas[2]=ang_z: [-pi:pi];
        '''
        thetas=np.zeros(3)
        thetas[0]=np.arctan2(-self[1,2],self[2,2])
        thetas[1]=np.arcsin(self[0,2])
        thetas[2]=np.arctan2(-self[0,1],self[0,0])
        if (lDegs):
            thetas=np.degrees(thetas)
        return thetas
                
if (__name__=="__main__"):
    lDebug=False
    # myMat=RotMat(myAng=60,myAxis=3,lDegs=True,lDebug=lDebug)
    # myMat=RotMat.ConcatenatedRotMatrices(myAngs=[60,30,90],myAxes=[1,2,3],\
    #                               lDegs=True,lDebug=lDebug)
    # print(myMat.det())
    # print(myMat.echo())
    # print(myMat.inv().echo())
    myMat=RotMat.ConcatenatedRotMatrices(myAngs=[60,30,90],myAxes=[3,2,1],\
                                  lDegs=True,lDebug=lDebug)
    print(myMat.echo())
    print(myMat.GetGimbalAngles())
    print(myMat.inv().echo())
    print(myMat.inv().GetGimbalAngles())
    print(myMat.mulMat(myMat.inv()).echo())
    myMat=RotMat.ConcatenatedRotMatrices(myAngs=[60,30,90],myAxes=[1,2,3],\
                                  lDegs=True,lDebug=lDebug)
    print(myMat.echo())
    print(myMat.GetGimbalAngles())
    print(myMat.inv().echo())
    print(myMat.inv().GetGimbalAngles())
