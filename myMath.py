# classes for managing math aspects
# A. Mereghetti, 2023/09/29
# python version: >= 3.8.10;

import numpy as np

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
        'method to multiply squared matrices of the same number of rows'
        if (self.nDim!=RMat.nDim):
            print("cannot multiply a %dx%d matrix times a %dx%d matrix"%\
                  (self.nDim,self.nDim,RMat.nDim,RMat.nDim))
            exit(1)
        newMat=Matrix(self.nDim)
        for ii in range(self.nDim):
            for jj in range(self.nDim):
                newMat[ii,jj]=sum([ self[ii,kk]*RMat[kk,jj]\
                                        for kk in range(self.nDim)])
        if (lDebug):
            print("Matrix.mulMat():")
            print(self.echo())
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
            print("Matrix.mulArr():",out)
        return out

    def mulSca(self,mySca):
        'method to multiply a squared matrix by an array'
        newMat=Matrix(self.nDim)
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

    def MinMat(self):
        'calculate matrix of minors'
        newMat=Matrix(self.nDim)
        for ii in range(self.nDim):
            for jj in range(self.nDim):
                newMat[ii,jj]=self.Minor(ii,jj).det()
        return newMat

    def CofactMat(self):
        'calculate matrix of co-factors'
        newMat=Matrix(self.nDim)
        for ii in range(self.nDim):
            for jj in range(self.nDim):
                newMat[ii,jj]=self[ii,jj]
                if ((ii+jj)%2==1):
                    newMat[ii,jj]=-newMat[ii,jj]
        return newMat

    def AdjugateMat(self):
        'calculate adjugate matrix, i.e. the transposed'
        newMat=Matrix(self.nDim)
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
                det=det+self[0,ii]*self.Minor(0,ii).det()
            return det

    def inv(self):
        '''
        invert the matrix - based on:
        https://www.mathsisfun.com/algebra/matrix-inverse-minors-cofactors-adjugate.html
        '''
        newMat=self.MinMat()
        newMat=newMat.CofactMat()
        newMat=newMat.AdjugateMat()
        return newMat.mulSca(self.det())

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
    def __init__(self,myAng=0.0,myAxis=3,lDegs=True,lDebug=True):
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
        if (len(myAngs)!=len(myAxes)):
            print("number of angles (%d) and axes (%d) do not match!"%(\
                  len(myAngs),len(myAxes)))
            exit(1)
        # go through the array on angles and axes to define the overall
        #   matrix transformation
        rotMatrices=[]
        if (lDebug):
            print("creating rotation matrices...")
        for iTrasf in range(len(myAngs)):
            if (not lDegs):
                angDegs=np.degrees(myAngs[iTrasf])
                angRads=myAngs[iTrasf]
            else:
                angDegs=myAngs[iTrasf]
                angRads=np.radians(myAngs[iTrasf])
            if (lDebug):
                print("...creating rotation matrix: %g degs around %d axis..."%\
                      (angDegs,myAxes[iTrasf]))
            rotMatrices.append(RotMat(myAng=myAngs[iTrasf],myAxis=myAxes[iTrasf],\
                          lDegs=lDegs,lDebug=lDebug))
        if (lDebug):
            print("concatenating them...")
        newMat=UnitMat()
        for iTrasf in range(len(myAngs)-1,-1,-1):
            newMat=newMat.mulMat(rotMatrices[iTrasf],lDebug=lDebug)
            
        return newMat
                
if (__name__=="__main__"):
    lDebug=True
    myMat=RotMat(myAng=60,myAxis=3,lDegs=True,lDebug=lDebug)
    myMat=RotMat.ConcatenatedRotMatrices(myAngs=[90,90,90],myAxes=[1,2,3],\
                                  lDegs=True,lDebug=lDebug)
    print(myMat.det())
    print(myMat.echo())
    print(myMat.inv().echo())
