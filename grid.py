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
       orientation (3D rotation matrix).
    Additionally, the sequence of rotations giving rise to the final
       orientation is also stored.
    Rotation axes are labelled as: x=1, y=2, z=3;
    '''
    def __init__(self,myP=[0.,0.,0.],myW=myMath.UnitMat,myAngs=[0.,0.,0.],myAxes=[1,2,3],myLab=""):
        self.P=myP
        self.W=myW
        self.angs=myAngs # [degs]
        self.axes=myAxes
        self.label=myLab

    def ret(self,what):
        'a return method'
        if (what.upper()=="POINT"):
            return self.P
        elif (what.upper()=="MATRIX"):
            return self.W
        elif (what.upper().startswith("ANGLE")):
            return self.angs
        elif (what.upper().startswith("AX")):
            return self.axes
        elif (what.upper()=="LABEL"):
            return self.label
        else:
            print("Location.ret(): unknown request %s!"%(what))
            exit(1)

    def set(self,what,value):
        'a set method'
        if (what.upper()=="POINT"):
            self.P=value
        elif (what.upper()=="MATRIX"):
            self.W=value
        elif (what.upper().startswith("ANGLE")):
            self.angs=value
        elif (what.upper().startswith("AX")):
            self.axes=value
        elif (what.upper()=="LABEL"):
            self.label=value
        else:
            print("Location.set(): unknown request %s!"%(what))
            exit(1)

    def echo(self,myFmt="% 13.6E",mySep="; "):
        'echo method, to write out self.P and self.W'
        buf="P=["
        for tmpX in self.P:
            buf=buf+" "+myFmt%(tmpX)
        buf=buf+" ]"+mySep+"W=|"
        for jj in range(self.W.nDim):
            for ii in range(self.W.nDim):
                buf=buf+" "+myFmt%(self.W[ii,jj])
            if (jj<self.W.nDim-1):
                if ("\n" in mySep):
                    buf=buf+" | "+mySep+"  |"
                else:
                    buf=buf+" "+mySep
        buf=buf+" | "+mySep+"label=\"%s\""%(self.label)+mySep
        return buf

    def ComputeRotMat(self,lDebug=False):
        'compute rotation matrix based on angles and axes'
        if (len(self.angs)!=len(self.axes)):
            print("Location.ComputeRotMat(): number of angles (%d) and axes (%d) differ!"\
                  %(len(self.angs),len(self.axes)))
            exit(1)
        self.W=myMath.RotMat.ConcatenatedRotMatrices( \
               myAngs=self.angs,myAxes=self.axes, \
               lDegs=True,lDebug=lDebug)
        
class Grid:
    '''
    This class describes a grid of points/orientations onto which a prototype
       geometry should be cloned. Coordinates/orientations are meant to move the
       centre of the prototype (expected to be at the origin of the FLUKA ref
       system) to the actual position/orientation.
    The class is basically a list of locations.
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

    def AddLoc(self,myLoc):
        self.locs.append(myLoc)

    def echo(self,myFmt="% 13.6E"):
        buf=""
        for myLoc in self.locs:
            buf=buf+myLoc.echo(myFmt=myFmt)+"\n"
        return buf

    def ret(self,what="point",iEl=-1):
        return self.locs[iEl].ret(what=what)

    @staticmethod
    def SphericalShell(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,\
                       lDebug=False):
        '''
        This method defines a grid describing a spherical shell.

        The grid is described in spherical coordinates:
        - r [cm];
        - theta [degs]: angle in yz-plane (positive when pointing towards y>0);
                        range: [-90:90] degs;
        - phi [degs]: angle in xz-plane (positive when pointing towards x>0).
                        range: [-180:180] degs;
        The grid starts around the z-axis.

        NR, NT, NP are number of points along grid (cell centers);
          therefore, the number of steps are NR-1, NT-1, NP-1;
        '''
        # some checks
        oRmin,oRmax,oNR,oTmin,oTmax,oNT,oPmin,oPmax,oNP=\
            CheckSphericalRange(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP)
        
        # define mesh ranges
        if (oNT==1):
            oTmean=(oTmax+oTmin)/2.
            oTmin, oTmax = oTmean, oTmean
        if (oNP==1):
            oPmean=(oPmax+oPmin)/2.
            oPmin, oPmax = oPmean, oPmean
        # define unique shell values of radius and angles
        # NB: np.linspace: num is number of points, extremes are included
        RRs=np.linspace(oRmin,oRmax,num=oNR)
        TTs=np.linspace(oTmin,oTmax,num=oNT)
        PPs=np.linspace(oPmin,oPmax,num=oNP)
        
        # create set of locations (coordinates and orientations)
        newGrid=Grid()
        for tmpR in RRs:
            for tmpT in TTs:
                for tmpP in PPs:
                    myLab="R[cm]=%g,theta[deg]=%g,phi[deg]=%g"%(tmpR,tmpT,tmpP)
                    myLoc=Location(myAngs=[-tmpT,tmpP],myAxes=[1,2],myLab=myLab)
                    myLoc.ComputeRotMat(lDebug=lDebug)
                    myW=myLoc.ret("MATRIX")
                    myLoc.set("POINT",myW.mulArr([0.0,0.0,tmpR],lDebug=lDebug))
                    newGrid.AddLoc(myLoc)
                    
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
    This class implements the coordinates of the hive for a grid;
      this class is actually a parent class implementing general-purpose
      methods and fields.
    The actual structure of the class depends on the type - hence, the
      actual hives are implemented as child classes.
    '''
    def __init__(self,hType=None):
        if (hType is not None):
            if (hType.upper()!="SPHERE"):
                print("Hive.__init__(): cannot create a hive of type %s!"%(\
                    hType))
                exit(1)
        self.hType=hType

class SphericalHive(Hive):
    '''
    Class for creating a spherical hive - actually a spherical shell,
      for which a sphere is a special case.

    All the input info refer to the respective spherical grid of objects.
      Please see the docstring of the Grid.SphericalShell() method.
    The number of interfaces in the hive along r, theta and phy
      are NR+1, NT+1 and NP+1.

    The class is actually made of three arrays, which state radius and
       angles of the cutting spheres/planes building up the hive.
    '''

    def __init__(self,Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP,lDebug=False):
        # some checks
        oRmin,oRmax,oNR,oTmin,oTmax,oNT,oPmin,oPmax,oNP=\
            CheckSphericalRange(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP)
        
        # define mesh ranges
        Hive.__init__(self,"SPHERE")
        hdR=(oRmax-oRmin)
        if (NR==1):
            hRmax=oRmax-hdR/2.; hRmin=oRmin-hdR/2.
        else:
            hdR=hdR/(oNR-1)
            hRmax=oRmax+hdR/2.; hRmin=oRmin-hdR/2.
        hdT=(oTmax-oTmin)
        if (NT==1):
            tMean=(oTmax+oTmin)/2.
            hTmax=tMean+hdT/2.; hTmin=tMean-hdT/2.
        else:
            hdT=hdT/(oNT-1)
            hTmax=oTmax+hdT/2.; hTmin=oTmin-hdT/2.
        hdP=(oPmax-oPmin)
        if (oNP==1):
            pMean=(oPmax+oPmin)/2.
            hPmax=pMean+hdP/2.; hPmin=pMean-hdP/2.
        else:
            hdP=hdP/(oNP-1)
            hPmax=oPmax+hdP/2.; hPmin=oPmin-hdP/2.
            
        # define unique shell values of radius and angles
        # NB: np.linspace: num is number of points, extremes are included
        self.RRs=np.linspace(hRmin,hRmax,num=NR+1)
        self.TTs=np.linspace(hTmin,hTmax,num=NT+1)
        self.PPs=np.linspace(hPmin,hPmax,num=NP+1)

        if (lDebug):
            print("SphericalHive(): R[cm]=[%g:%g:%g] - N=%d;"%(hRmin,hdR,hRmax,oNR+1))
            print("SphericalHive(): theta[deg]=[%g:%g:%g] - N=%d;"%(hTmin,hdT,hTmax,oNT+1))
            print("SphericalHive(): phi[deg]=[%g:%g:%g] - N=%d;"%(hPmin,hdP,hPmax,oNP+1))
            print("")
            print("SphericalHive(): RRs:",self.ret("RRs"))
            print("SphericalHive(): TTs:",self.ret("TTs"))
            print("SphericalHive(): PPs:",self.ret("PPs"))
            
    @staticmethod
    def SphericalHive_OneLayer(R,dR,Tmax,NT,Pmax,NP,lDebug=False):
        '''
        Special case of a spherical hive, with:
        - only 1 radial layer;
        - symmetric angular ranges;
        '''
        return SphericalHive(R,R+dR,1,-Tmax,Tmax,NT,-Pmax,Pmax,NP,lDebug=lDebug)

    def ret(self,what="all"):
        if (what.upper()=="RRS"):
            return self.RRs
        elif (what.upper()=="TTS"):
            return self.TTs
        elif (what.upper()=="PPS"):
            return self.PPs
        elif (what.upper()=="ALL"):
            return self.ret(what="RRs"), \
                   self.ret(what="TTs"), \
                   self.ret(what="PPs")
        else:
            print("SphericalHive.ret(): what to return? %s"%(what))

def CheckSphericalRange(Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP):
    'See the docstring of the Grid.SphericalShell() method'
    oRmin,oRmax,oNR,oTmin,oTmax,oNT,oPmin,oPmax,oNP = \
        Rmin,Rmax,NR,Tmin,Tmax,NT,Pmin,Pmax,NP
    if (oRmin>oRmax):
        oRmin, oRmax = oRmax, oRmin
        if (lDebug):
            print("grid.CheckSphericalRange(): Rmin,Rmax sorted in increasing order;")
    if (oTmin>oTmax):
        oTmin, oTmax = oTmax, oTmin
        if (lDebug):
            print("grid.CheckSphericalRange(): Tmin,Tmax sorted in increasing order;")
    if (oPmin>oPmax):
        oPmin, oPmax = oPmax, oPmin
        if (lDebug):
            print("grid.CheckSphericalRange(): Pmin,Pmax sorted in increasing order;")
    if (oRmin<=0):
        print("grid.CheckSphericalRange(): Rmin<=0!")
        exit(1)
    if ( oTmin<-90 or oTmax>90 ):
        print("grid.CheckSphericalRange(): allowed range for theta: [-90:90] - given: [%g:%g]!"%(oTmin,oTmax))
        exit(1)
    if ( oPmin<-180 or oPmax>180 ):
        print("grid.CheckSphericalRange(): allowed range for phi: [-180:180] - given: [%g:%g]!"%(oPmin,oPmax))
        exit(1)
    if (oNR<=0):
        print("grid.CheckSphericalRange(): NR<=0!")
        exit(1)
    if (oNT<0):
        print("grid.CheckSphericalRange(): NT<=0!")
        exit(1)
    if (oNP<0):
        print("grid.CheckSphericalRange(): NP<=0!")
        exit(1)
    return oRmin,oRmax,oNR,oTmin,oTmax,oNT,oPmin,oPmax,oNP

if ( __name__ == "__main__" ):
    # perform some tests
    lDebug=True
    R=500
    dR=50
    NR=1
    Tmax=3    # theta [degs] --> range: -Tmax:Tmax
    NT=0      # number of steps (i.e. entities)
    Pmax=2    # phi [degs] --> range: -Pmax:Pmax
    NP=0      # number of steps (i.e. entities)
    myGrid=Grid.SphericalShell(R,R+dR,NR,-Tmax,Tmax,NT,-Pmax,Pmax,NP,lDebug=lDebug)
    # print(DefHiveBoundaries_SphericalShell_OneLayer(R,50,Tmax,NT,Pmax,NP))
    for loc in myGrid:
        print(loc.ret("Point"))
    myHive=SphericalHive(R,R+dR,NR,-Tmax,Tmax,NT,-Pmax,Pmax,NP,lDebug=lDebug)
    
