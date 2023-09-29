# classes for managing math aspects
# A. Mereghetti, 2023/09/29
# python version: >= 3.8.10;

import numpy as np

class Matrix:
    '''
    very simple class coding a 2D matrix
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
            for jj in range(self.nDim):
                buf=buf+" "+myFmt%(self[ii,jj])
            buf=buf+"\n"
        return buf

class RotMat(Matrix):
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
    def __init__(self,myAng=0.0,myAxis=3,lDegs=True):
        # first, initialise matrix as unit matrix;
        Matrix.__init__(self,nDim=3)
        for ii in range(self.nDim):
            self[ii,ii]=1.0
            
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
                
if (__name__=="__main__"):
    myMat=RotMat(myAng=30,myAxis=0)
    print(myMat.echo())
