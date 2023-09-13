# Grid class
# A. Mereghetti, 2023/09/12
# python version: >= 3.8.10;
#
# This class describes the grid of points/orientations of crystals in the full
#    CALO geometry. Coordinates/orientations move the centre of each crystal,
#    defined at the origin of the FLUKA ref system of the geometry of the single
#    crystal.
#
# Reference frame (FLUKA-like):
# - z-axis along beam;
# - y-axis anti-gravitational;
# - x-axis to get a right-hand ref sys.
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
    def __init__(self,Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP):
        self.RRs=np.linspace(Rmin,Rmax,num=NR)
        self.TTs=np.linspace(Tmin,Tmax,num=NT)
        self.PPs=np.linspace(Pmin,Pmax,num=NP)

    @staticmethod
    def SymOneLayer(R,Tmax,NT,Pmax,NP):
        return GRID(R,R,1,-Tmax,Tmax,NT,-Pmax,Pmax,NP)
    
    def CreateHive(self,dR=None,dT=None,dP=None):
        # get dR
        if (len(self.RRs)==1):
            if (dR is None):
                write("dR NOT defined!")
                exit(1)
        else:
            if (dR is None):
                dR=self.RRs[1]-self.RRs[0]
        # get dT
        if (len(self.TTs)==1):
            if (dT is None):
                write("dT NOT defined!")
                exit(1)
        else:
            if (dT is None):
                dR=self.TTs[1]-self.TTs[0]
        # get dP
        if (len(self.PPs)==1):
            if (dP is None):
                write("dP NOT defined!")
                exit(1)
        else:
            if (dP is None):
                dP=self.PPs[1]-self.PPs[0]
        RRs=np.append(self.RRs-dR/2.,self.RRs[-1]+dR/2.)
        TTs=np.append(self.TTs-dT/2.,self.TTs[-1]+dT/2.)
        PPs=np.append(self.PPs-dP/2.,self.PPs[-1]+dP/2.)
        return RRs, TTs, PPs

