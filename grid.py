# Grid class
# A. Mereghetti, 2023/10/28
# python version: >= 3.8.10;
#
# This class describes a grid of points/orientations onto which a prototype
#    geometry should be cloned. Coordinates/orientations are meant to move the
#    centre of the prototype (expected to be at the origin of the FLUKA ref
#    system) to the actual position/orientation.
#
# Reference frame (FLUKA-like):
# - z-axis along beam;
# - y-axis anti-gravitational;
# - x-axis to get a right-handed ref sys.

import numpy as np
import myMath

class LOCATION:
    '''
    a very simple class storing a position (3D array) and an orientation
       (3D rotation matrix)
    '''
    def __init__(self,myP=[0.,0.,0.],myW=myMath.UnitMat,myLab=""):
        self.P=myP
        self.W=myW
        self.label=myLab

    def echo(self,myFmt="% .6E"):
        buf="P=["
        for tmpX in self.P:
            buf=buf+" "+myFmt%(tmpX)
        buf=buf+" ]; W=|"
        for jj in range(self.W.nDim):
            for ii in range(self.W.nDim):
                buf=buf+" "+myFmt%(self.W[ii,jj])
            if (jj<self.W.nDim-1):
                buf=buf+" ;"
        buf=buf+" | label=\"%s\";"%(self.label)
        return buf
        
class GRID:
    def __init__(self):
        self.locs=[] # list of LOCATIONs
        
    def AddLoc(self,myP=[0.,0.,0.],myW=myMath.UnitMat,myLab=""):
        self.locs.append(LOCATION(myP=myP,myW=myW,myLab=myLab))

    def echo(self,myFmt="% .6E"):
        buf=""
        for myLoc in self.locs:
            buf=buf+myLoc.echo(myFmt=myFmt)+"\n"
        return buf

    @staticmethod
    def SphericalShell(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,\
                       lDebug=False):
        '''
        This method defines a grid describing a spherical shell.

        The grid is described in spherical coordinates:
        - r [cm];
        - theta [degs]: angle in yz-plane (positive when pointing towards y>0);
        - phi [degs]: angle in xz-plane (positive when pointing towards x>0).

        The grid is centred around the z-axis and symmetric in the two
          directions, up/down (y>0/y<0), left/right (x>0/x<0).

        NR, NT, NP are number of points along grid;
        therefore, the number of steps are NR-1, NT-1, NP-1;
        '''
        # define unique shell values of radius and angles
        # NB: np.linspace: num is number of points
        RRs=np.linspace(Rmin,Rmax,num=NR)
        TTs=np.linspace(Tmin,Tmax,num=NT)
        PPs=np.linspace(Pmin,Pmax,num=NP)
        
        # create set of locations (coordinates and orientations)
        newGrid=GRID()
        for tmpR in RRs:
            for tmpT in TTs:
                for tmpP in PPs:
                    myLab="R[cm]=%g,theta[deg]=%g,phi[deg]=%g"%(tmpR,tmpT,tmpP)
                    myW=myMath.RotMat.ConcatenatedRotMatrices( \
                        myAngs=[-tmpT,tmpP],myAxes=[1,2], \
                        lDegs=True,lDebug=lDebug)
                    myP=myW.mulArr([0.0,0.0,tmpR],lDebug=lDebug)
                    newGrid.AddLoc(myP=myP,myW=myW,myLab=myLab)
                    
        if (lDebug):
            print("GRID.SphericalShell(): RRs:",RRs)
            print("GRID.SphericalShell(): TTs:",TTs)
            print("GRID.SphericalShell(): PPs:",PPs)
            print(newGrid.echo())
            
        return newGrid
    
    @staticmethod
    def SphericalShell_OneLayer(R,Tmax,NT,Pmax,NP,lDebug=False):
        return GRID.SphericalShell(R,R,1,-Tmax,Tmax,NT,-Pmax,Pmax,NP,\
                                   lDebug=lDebug)

def DefHiveBoundaries_SphericalShell(Rmin,Rmax,dR,Tmin,Tmax,dT,Pmin,Pmax,dP):
    '''
    All input data in interface refer to the hive, not to the grid of objects
       therein contained!
    '''
    print("defining hive boundaries for a spherical shell:")
    print("* R[cm]=[%g:%g:%g];"%(Rmin,dR,Rmax))
    print("* theta[deg]=[%g:%g:%g];"%(Tmin,dT,Tmax))
    print("* phi[deg]=[%g:%g:%g];"%(Pmin,dP,Pmax))
    RRs=np.arange(Rmin,Rmax+dR,dR)
    TTs=np.arange(Tmin,Tmax+dT,dT)
    PPs=np.arange(Pmin,Pmax+dP,dP)
    return RRs, TTs, PPs
    
def DefHiveBoundaries_SphericalShell_OneLayer(R,dR,Tmax,NT,Pmax,NP):
    dT=2*Tmax/NT
    dP=2*Pmax/NP
    return DefHiveBoundaries_SphericalShell(R-dR/2,R+dR/2,dR,\
                                            -Tmax-dT/2,Tmax+dT/2,dT,\
                                            -Pmax-dP/2,Pmax+dP/2,dP)

if ( __name__ == "__main__" ):
    # perform some tests
    lDebug=True
    R=500
    Tmax=3    # theta [degs] --> range: -Tmax:Tmax
    NT=2      # number of steps (i.e. entities)
    Pmax=2    # phi [degs] --> range: -Pmax:Pmax
    NP=2      # number of steps (i.e. entities)
    GRID.SphericalShell_OneLayer(R,Tmax,NT,Pmax,NP,lDebug=lDebug)
    print(DefHiveBoundaries_SphericalShell_OneLayer(R,50,Tmax,NT,Pmax,NP))
