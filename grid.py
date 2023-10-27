# Grid class
# A. Mereghetti, 2023/09/12
# python version: >= 3.8.10;
#
# This class describes the grid of points/orientations of crystals in the full
#    CALO geometry. Coordinates/orientations move the centre of each crystal,
#    defined at the origin of the FLUKA ref system of the single crystal, to the
#    actual position/orientation.
#
# Reference frame (FLUKA-like):
# - z-axis along beam;
# - y-axis anti-gravitational;
# - x-axis to get a right-handed ref sys.
# 
# The grid is described in spherical coordinates:
# - r [cm];
# - theta [degs]: angle in yz-plane (positive when pointing towards y>0);
# - phi [degs]: angle in xz-plane (positive when pointing towards x>0).
#
# the grid is centred around the z-axis and symmetric in the two directions
#   (up/down, left/right, i.e. x>0/x<0).
# 

import numpy as np

class GRID:
    def __init__(self,Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,lDebug=False):
        '''
        NR, NT, NP are number of points along grid;
        therefore, the number of steps are NR-1, NT-1, NP-1;
        '''
        # np.linspace: num is number of points
        self.RRs=np.linspace(Rmin,Rmax,num=NR)
        self.TTs=np.linspace(Tmin,Tmax,num=NT)
        self.PPs=np.linspace(Pmin,Pmax,num=NP)
        if (lDebug):
            print("GRID.__init__(): self.RRs:",self.RRs)
            print("GRID.__init__(): self.TTs:",self.TTs)
            print("GRID.__init__(): self.PPs:",self.PPs)

    @staticmethod
    def OneLayer(R,Tmax,NT,Pmax,NP,lDebug=True):
        return GRID(R,R,1,-Tmax,Tmax,NT,-Pmax,Pmax,NP,lDebug=lDebug)
    
    def DefineHiveBoundaries(self,lDebug=False,dR=None,dT=None,dP=None):
        # get dR
        if (len(self.RRs)==1):
            if (dR is None):
                print("dR NOT defined!")
                exit(1)
        else:
            if (dR is None):
                dR=self.RRs[1]-self.RRs[0]
        # get dT
        if (len(self.TTs)==1):
            if (dT is None):
                print("dT NOT defined!")
                exit(1)
        else:
            if (dT is None):
                dT=self.TTs[1]-self.TTs[0]
        # get dP
        if (len(self.PPs)==1):
            if (dP is None):
                print("dP NOT defined!")
                exit(1)
        else:
            if (dP is None):
                dP=self.PPs[1]-self.PPs[0]
        RRs=np.append(self.RRs-dR/2.,self.RRs[-1]+dR/2.)
        TTs=np.append(self.TTs-dT/2.,self.TTs[-1]+dT/2.)
        PPs=np.append(self.PPs-dP/2.,self.PPs[-1]+dP/2.)
        if (lDebug):
            print("GRID.DefineHiveBoundaries(): RRs,dR:",RRs,dR)
            print("GRID.DefineHiveBoundaries(): TTs,dT:",TTs,dT)
            print("GRID.DefineHiveBoundaries(): PPs,dP:",PPs,dP)
        return RRs, TTs, PPs

if ( __name__ == "__main__" ):
    # perform some tests
    lDebug=True
    R=500
    Tmax=27   # theta [degs] --> range: -Tmax:Tmax
    NT=10     # number of steps (i.e. entities)
    Pmax=20   # phi [degs] --> range: -Pmax:Pmax
    NP=21     # number of steps (i.e. entities)
    myGrid=GRID.OneLayer(R,Tmax,NT,Pmax,NP,lDebug=lDebug)
    RRs,TTs,PPs=myGrid.DefineHiveBoundaries(lDebug=lDebug,dR=100)
