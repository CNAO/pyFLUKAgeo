# Grid class
# A. Mereghetti, 2023/10/28
# python version: >= 3.8.10;
#
# Reference frame (FLUKA-like):
# - z-axis along beam;
# - y-axis anti-gravitational;
# - x-axis to get a right-handed ref sys.

import numpy as np
import myMath

class Location:
    '''
    a very simple class storing a single position (3D array) and a single
       orientation (3D rotation matrix)
    '''
    def __init__(self,myP=[0.,0.,0.],myW=myMath.UnitMat,myLab=""):
        self.P=myP
        self.W=myW
        self.label=myLab

    def ret(self,what):
        if (what.upper()=="POINT"):
            return self.P
        elif (what.upper()=="MATRIX"):
            return self.W
        elif (what.upper()=="LABEL"):
            return self.label

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
        
class Grid:
    '''
    This class describes a grid of points/orientations onto which a prototype
       geometry should be cloned. Coordinates/orientations are meant to move the
       centre of the prototype (expected to be at the origin of the FLUKA ref
       system) to the actual position/orientation.
    '''
    def __init__(self):
        self.locs=[] # list of Locations

    def __len__(self):
        return len(self.locs)

    def __iter__(self):
        self.current_index=0
        return self
    
    def __next__(self):
        '''
        from https://blog.finxter.com/python-__iter__-magic-method/
        '''
        if self.current_index < len(self):
            x = self.locs[self.current_index]
            self.current_index += 1
            return x
        raise StopIteration

    def AddLoc(self,myP=[0.,0.,0.],myW=myMath.UnitMat,myLab=""):
        self.locs.append(Location(myP=myP,myW=myW,myLab=myLab))

    def echo(self,myFmt="% .6E"):
        buf=""
        for myLoc in self.locs:
            buf=buf+myLoc.echo(myFmt=myFmt)+"\n"
        return buf

    def ret(self,myWhat="point",iEl=-1):
        if (myWhat.upper()=="POINT"):
            return self.locs[iEl].P
        elif (myWhat.upper()=="REF"):
            return self.locs[iEl].W
        elif (myWhat.upper()=="LAB"):
            return self.locs[iEl].label
        else:
            print("Grid.ret(): unknown request %s!"%(myWhat))
            exit(1)

    @staticmethod
    def SphericalShell(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,\
                       lDebug=False):
        '''
        This method defines a grid describing a spherical shell.

        The grid is described in spherical coordinates:
        - r [cm];
        - theta [degs]: angle in yz-plane (positive when pointing towards y>0);
        - phi [degs]: angle in xz-plane (positive when pointing towards x>0).
        The grid starts around the z-axis.

        NR, NT, NP are number of points along grid (cell centers);
          therefore, the number of steps are NR-1, NT-1, NP-1;
        '''
        # define unique shell values of radius and angles
        # NB: np.linspace: num is number of points
        RRs=np.linspace(Rmin,Rmax,num=NR)
        TTs=np.linspace(Tmin,Tmax,num=NT)
        PPs=np.linspace(Pmin,Pmax,num=NP)
        
        # create set of locations (coordinates and orientations)
        newGrid=Grid()
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
            print("Grid.SphericalShell(): RRs:",RRs)
            print("Grid.SphericalShell(): TTs:",TTs)
            print("Grid.SphericalShell(): PPs:",PPs)
            print(newGrid.echo())
            
        return newGrid
    
    @staticmethod
    def SphericalShell_OneLayer(R,Tmax,NT,Pmax,NP,lDebug=False):
        '''
        Special case of a spherical shell instance, with:
        - only 1 radial layer;
        - symmetric angular ranges;
        '''
        return Grid.SphericalShell(R,R,1,-Tmax,Tmax,NT,-Pmax,Pmax,NP,\
                                   lDebug=lDebug)

class Hive:
    '''
    This class implements the coordinates of the hive for the grid.
    The actual structure of the class depends on the type.
    '''
    def __init__(self,hType=None):
        if (hType is not None):
            if (hType.upper()!="SPHERE"):
                print("Hive.__init__(): cannot create a hive of type %s!"%(\
                    hType))
                exit(1)
        self.hType=hType

    def ret(self,myWhat="all"):
        if (self.hType.upper()=="SPHERE"):
            if (myWhat.upper()=="RRS"):
                return self.RRs
            elif (myWhat.upper()=="TTS"):
                return self.TTs
            elif (myWhat.upper()=="PPS"):
                return self.PPs
            elif (myWhat.upper()=="ALL"):
                return self.ret(myWhat="RRs"), \
                       self.ret(myWhat="TTs"), \
                       self.ret(myWhat="PPs")
            else:
                print("Hive.ret(): what to return? %s"%(myWhat))

    @staticmethod
    def SphericalHive(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,lDebug=False):
        '''
        Method for creating a spherical hive - actually a spherical shell,
          for which a sphere is a special case.

        All the input info refer to the respective spherical grid of objects.
          Please see the docstring of the Grid.SphericalShell() method.
        The number of interfaces in the hive along r, theta and phy
          are NR+1, NT+1 and NP+1
        '''
        newHive=Hive("SPHERE")
        hdR=(Rmax-Rmin)
        if (NR==1):
            hRmax=Rmax-hdR/2.; hRmin=Rmin-hdR/2.
        else:
            hdR=hdR/(NR-1)
            hRmax=Rmax+hdR; hRmin=Rmin-hdR
        hdT=(Tmax-Tmin)
        if (NT==1):
            hTmax=Tmax-hdT/2.; hTmin=Tmin-hdT/2.
        else:
            hdT=hdT/(NT-1)
            hTmax=Tmax+hdT; hTmin=Tmin-hdT
        hdP=(Pmax-Pmin)
        if (NP==1):
            hPmax=Pmax-hdP/2.; hPmin=Pmin-hdP/2.
        else:
            hdP=hdP/(NP-1)
            hPmax=Pmax+hdP; hPmin=Pmin-hdP
        # define unique shell values of radius and angles
        if (lDebug):
            print("* R[cm]=[%g:%g:%g];"%(hRmin,hdR,hRmax))
            print("* theta[deg]=[%g:%g:%g];"%(hTmin,hdT,hTmax))
            print("* phi[deg]=[%g:%g:%g];"%(hPmin,hdP,hPmax))
        # NB: np.linspace: num is number of points
        newHive.RRs=np.linspace(hRmin,hRmax,num=NR+1)
        newHive.TTs=np.linspace(hTmin,hTmax,num=NT+1)
        newHive.PPs=np.linspace(hPmin,hPmax,num=NP+1)
        return newHive

    @staticmethod
    def SphericalHive_OneLayer(R,dR,Tmax,NT,Pmax,NP,lDebug=False):
    '''
    Special case of a spherical hive, with:
    - only 1 radial layer;
    - symmetric angular ranges;
    '''
    return Hive.SphericalHive(R,R+dR,1,-Tmax,Tmax,NT,-Pmax,Pmax,NP,lDebug=lDebug)

if ( __name__ == "__main__" ):
    # perform some tests
    lDebug=True
    R=500
    Tmax=3    # theta [degs] --> range: -Tmax:Tmax
    NT=2      # number of steps (i.e. entities)
    Pmax=2    # phi [degs] --> range: -Pmax:Pmax
    NP=2      # number of steps (i.e. entities)
    myGrid=Grid.SphericalShell_OneLayer(R,Tmax,NT,Pmax,NP,lDebug=lDebug)
    # print(DefHiveBoundaries_SphericalShell_OneLayer(R,50,Tmax,NT,Pmax,NP))
    for loc in myGrid:
        print(loc.ret("Point"))
