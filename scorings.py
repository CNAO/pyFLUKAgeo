# classes for managing FLUKA scoring cards and detectors
# A. Mereghetti, 2024/04/16
# python version: >= 3.8.10

from FLUKA import GeoObject, echoFloats
import numpy as np

class Scoring(GeoObject):
    '''
    A very basic class for handling scoring cards.
    For the time being, there is support only for:
    - the scoring is not really interpreted: the FLUKA definition is kept in
      memory, and manipulations are performed on-the-fly;
    - parsing/echoing, ALWAYS in FIXED format;
    - always 2 cards to fully define scoring;
    - NO comments between scoring cards defining the same scoring detector;
    - manipulation of unit;
    - support for AUXSCORE cards: only one AUXSCORE card per detector!
    '''
    def __init__(self,myName="",myComment="",scoType=""):
        GeoObject.__init__(self,myName=myName,myComment=myComment)
        self.definition=[] # array of strings storing the lines defining the scoring
        self.auxScoDef=""
        self.auxScoSDUM=""
        self.scoType=scoType

    def echo(self,what="all"):
        '''
        echo just the bare scoring cards and/or AUXSCORE card:
        * what="ALL": everything;
        * what="AUX": only AUXSCORE card;
        * what="COM": only comments;
        * what="SCO": only scoring cards;
        '''
        tmpBuf=""
        if (what.upper().startswith("ALL") or what.upper().startswith("COM")):
            tmpBuf=GeoObject.echoComm(self)
        if (what.upper().startswith("ALL") or what.upper().startswith("AUX")):
            if (self.hasAuxScoCard()):
                tmpBuf=tmpBuf+"%-10s%10s%10s%10s%10s%20s%-10s"%(\
                       "AUXSCORE",self.scoType,self.auxScoDef,"",self.echoName(),"",self.auxScoSDUM)
        if (what.upper().startswith("ALL") or what.upper().startswith("SCO")):
            if (what.upper().startswith("ALL") and self.hasAuxScoCard()):
                tmpBuf=tmpBuf+"\n"
            for myDef,mySdum,myEoL in zip(self.definition,[self.echoName(),"&"],["\n",""]):
                tmpBuf=tmpBuf+"%-10s%60s%-10s%-s"%(self.scoType,myDef,mySdum,myEoL)
        return tmpBuf

    @staticmethod
    def fromBuf(myBuffer,newScoDet=None):
        if (newScoDet is None): newScoDet=Scoring()
        for myLine in myBuffer.split("\n"):
            myLine=myLine.strip()
            if (myLine.startswith("*")):
                newScoDet.tailMe(myLine)
            elif (myLine.startswith("AUXSCORE")):
                newScoDet.auxScoDef=myLine[20:30]
                newScoDet.auxScoSDUM=myLine[70:].strip()
            else:
                newScoDet.definition.append(myLine[10:70])
                if (len(newScoDet.definition)==1):
                    # from first line defining USRBIN, get USRBIN name
                    newScoDet.rename(myLine[70:].strip(),lNotify=False)
        if (len(newScoDet.definition)!=2):
            print("Scoring.fromBuf(): something wrong when parsing USRBIN cards: got %d lines!"\
                  %(len(newScoDet.definition)))
            exit(1)
        if (len(newScoDet.name)==0):
            print("Scoring.fromBuf(): something wrong when parsing USRBIN cards: no name!")
            exit(1)
        return newScoDet

    def getUnit(self):
        return float(self.definition[0][20:30])
    def setUnit(self,myUnit,lFree=False):
        if (self.getUnit()<0):
            self.definition[0]=self.definition[0][ 0:20]+echoFloats(-abs(myUnit),lFree=lFree)[0]+ \
                               self.definition[0][30:]
        else:
            self.definition[0]=self.definition[0][ 0:20]+echoFloats( abs(myUnit),lFree=lFree)[0]+ \
                               self.definition[0][30:]

    def hasAuxScoCard(self):
        return len(self.auxScoDef)>0

class Usrbin(Scoring):
    '''
    A very basic class for handling USRBIN cards.
    In addition to Scoring class properties/methods, there is support only for:
    - a single ROTPRBIN card per transformation (1:1 mapping between ROTPRBIN and
      USRBIN cards); the mapping to the transformation is ALWAYS name-based;
      the ROTPRBIN card always PRECEEDs the respective USRBIN card;
    - manipulation of extremes and number of bins;
    '''
    def __init__(self,myName="",myComment=""):
        Scoring.__init__(self,myName=myName,myComment=myComment,scoType="USRBIN")
        self.TransfName=None

    def echo(self):
        '''take into account comment'''
        tmpBuf=Scoring.echo(self,what="com")
        auxBuf=Scoring.echo(self,what="aux")
        if (len(auxBuf)>0): tmpBuf=tmpBuf+auxBuf+"\n"
        if (self.isLinkedToTransform() and len(self.TransfName)>0):
            tmpBuf=tmpBuf+"%-10s%10s%10s%10s%10s\n"\
                    %("ROTPRBIN","",self.retTransformName(),"",self.echoName())
        return tmpBuf+Scoring.echo(self,what="sco")
                        
    @staticmethod
    def fromBuf(myBuffer):
        newBuffer=""
        RotPrBinBuff=""
        # fill in buffers
        for iLine,myLine in enumerate(myBuffer.split("\n")):
            if (myLine.strip().startswith("ROTPRBIN")):
                RotPrBinBuff=myLine.strip()
            else:
                if (iLine>0): newBuffer=newBuffer+"\n"
                newBuffer=newBuffer+myLine
        # parse buffers
        newUsrBin=Scoring.fromBuf(newBuffer,newScoDet=Usrbin())
        if (len(RotPrBinBuff)>0):
            newUsrBin.assignTransformName(RotPrBinBuff[20:30].strip())
            if (len(newUsrBin.TransfName)==0):
                print("Usrbin.fromBuf(): something wrong when parsing USRBIN cards: no ROT-DEFI card name!")
                exit(1)
        return newUsrBin

    def getBinType(self):
        return int(abs(float(self.definition[0][0:10]))+1E-4)
    def isSpecialBinning(self):
        binType=self.getBinType()
        return binType==2.0 or binType==12.0 or binType==8.0 or binType==18.0
    def isCartesianBinning(self):
        binType=self.getBinType()
        return binType==0.0 or binType==10.0
    def isCylindricalBinning(self):
        binType=self.getBinType()
        return binType==1.0 or binType==11.0

    def getNbins(self,axes="all"):
        'ONLY FIXED FORMAT DEFINITION'
        nBins=None
        if (not self.isSpecialBinning() ):
            if (isinstance(axes,list)):
                nBins=1
                for myAx in axes:
                    nBins=nBins*self.getNbins(axes=myAx)
            elif (isinstance(axes,str)):
                if (axes.upper()=="ALL"):
                    nBins=self.getNbins(axes=[1,2,3])
                elif (axes.upper()=="X"):
                    nBins=self.getNbins(axes=[1])
                elif (axes.upper()=="Y"):
                    nBins=self.getNbins(axes=[2])
                elif (axes.upper()=="Z"):
                    nBins=self.getNbins(axes=[3])
                else:
                    print("Usrbin.getNbins(): cannot get number of bins on ax specified as %s!"%(axes))
                    exit(1)
            elif (1<=axes and axes<=3):
                nBins=int(abs(float(self.definition[1][30+(axes-1)*10:40+(axes-1)*10]))+1E-4)
            else:
                print("Usrbin.getNbins(): cannot get number of bins on ax %d [1:3]!"%(axes))
                exit(1)
        return nBins
    def setNbins(self,nBins=[1],axes=[3],lFree=False):
        if (isinstance(axes,float) or isinstance(axes,int)): axes=[axes]
        elif (isinstance(axes,str)):
            if (axes.upper()=="ALL"):
                axes=[1,2,3]
            else:
                print("Usrbin.setNbins(): cannot set nBins to ax %s!"%(axes))
                exit(1)
        if (isinstance(nBins,float) or isinstance(nBins,int)): nBins=[nBins]
        if (len(axes)!=len(nBins)):
            print("Usrbin.setNbins(): len(axes)!=len(nBins)! axes=",axes,"; nBins=",nBins)
            exit(1)
        for myN,myAx in zip(nBins,axes):
            if (myAx<3):
                self.definition[1]=self.definition[1][0:30+10*(myAx-1)]+\
                   echoFloats(myN,lFree=lFree)[0]+self.definition[1][30+10*myAx:]
            elif (myAx>0):
                self.definition[1]=self.definition[1][0:30+10*(myAx-1)]+\
                   echoFloats(myN,lFree=lFree)[0]
            else:
                print("Usrbin.setNbins(): cannot set %d bins to axis %d!"%(myN,myAx))
                exit(1)

    def getExtremes(self,axes=3):
        'ONLY FIXED FORMAT DEFINITION'
        myMin=None; myMax=None
        if (isinstance(axes,list)):
            myMin=[]; myMax=[]
            for myAx in axes:
                tmpMin,tmpMax=self.getExtremes(axes=myAx)
                myMin.append(tmpMin); myMax.append(tmpMax)
        elif (isinstance(axes,float) or isinstance(axes,int)):
            if (axes==1):
                myMax=float(self.definition[0][30+(axes-1)*10:30+axes*10])
            elif (axes==2):
                if (self.isCylindricalBinning()):
                    myMax=np.pi
                else:
                    myMax=float(self.definition[0][30+(axes-1)*10:30+axes*10])
            elif (axes==3):
                myMax=float(self.definition[0][30+(axes-1)*10:          ])
            else:
                print("Usrbin.getExtremes(): cannot get extremes of bin on axis %d!"%(axes))
                exit(1)
            if (axes==2 and self.isCylindricalBinning()):
                myMin=-np.pi
            else:
                myMin=float(self.definition[1][ (axes-1)*10:axes*10])
        elif (isinstance(axes,str)):
            if (axes.upper()=="X"):
                myMin,myMax=self.getExtremes(axes=1)
            elif (axes.upper()=="Y"):
                myMin,myMax=self.getExtremes(axes=2)
            elif (axes.upper()=="Z"):
                myMin,myMax=self.getExtremes(axes=3)
            else:
                print("Usrbin.getExtremes(): cannot get extremes of bin on axis %s!"%(axes))
                exit(1)
        return myMin,myMax
    def setExtremes(self,myMin,myMax,axes=3,lFree=False):
        if (isinstance(myMin,list) and isinstance(myMax,list) and isinstance(axes,list)):
            if (len(myMin)!=len(myMax) or len(myMin)!=len(axes)):
                print("Usrbin.setExtremes(): cannot set min=",myMin,",max=",myMax,"axes",axes)
                exit(1)
            for tMin,tMax,tAx in zip(myMin,myMax,axes):
                self.setExtremes(tMin,tMax,axes=tAx,lFree=lFree)
        elif (not isinstance(myMin,list) and not isinstance(myMax,list) and not isinstance(axes,list)):
            myMinStr=echoFloats(myMin,lFree=lFree)[0]
            myMaxStr=echoFloats(myMax,lFree=lFree)[0]
            if (isinstance(axes,float) or isinstance(axes,int)):
                if (axes==1):
                    self.definition[0]=self.definition[0][:30+(axes-1)*10]+\
                                       myMaxStr+\
                                       self.definition[0][30+axes*10:]
                elif (axes==2):
                    if (self.isCylindricalBinning()):
                        print("...cylindrical binning: cannot set extremes on second axis!")
                        print("   skipping request;")
                    else:
                        self.definition[0]=self.definition[0][:30+(axes-1)*10]+\
                                           myMaxStr+\
                                           self.definition[0][30+axes*10:]
                elif (axes==3):
                    self.definition[0]=self.definition[0][:30+(axes-1)*10]+\
                                       myMaxStr
                else:
                    print("Usrbin.setExtremes(): cannot set extremes of bin on axis %d!"%(axes))
                    exit(1)
                self.definition[1]=self.definition[1][:(axes-1)*10]+\
                                   myMinStr+\
                                   self.definition[1][axes*10:]
            elif (isinstance(axes,str)):
                if (axes.upper()=="X"):
                    self.setExtremes(myMin,myMax,axes=1,lFree=lFree)
                elif (axes.upper()=="Y"):
                    self.setExtremes(myMin,myMax,axes=2,lFree=lFree)
                elif (axes.upper()=="Z"):
                    self.setExtremes(myMin,myMax,axes=3,lFree=lFree)
                else:
                    print("Usrbin.setExtremes(): cannot set extremes of bin on axis %s!"%(axes))
                    exit(1)
        else:
            print("Usrbin.setExtremes(): mixed arrays!")
            exit(1)

    def move(self,myCoord,axes=3,lAbs=True,lFree=False):
        if ((isinstance(myCoord,list) or isinstance(myCoord,np.ndarray)) and isinstance(axes,list)):
            if (len(myCoord)!=len(axes)):
                print("Usrbin.move(): cannot set absPos=",myCoord,"axes",axes)
                exit(1)
            for tCoord,tAx in zip(myCoord,axes):
                if (tCoord!=0.0): self.move(tCoord,axes=tAx,lAbs=lAbs,lFree=lFree)
        elif (not isinstance(myCoord,list) and not isinstance(myCoord,np.ndarray) and not isinstance(axes,list)):
            currMin,currMax=self.getExtremes(axes=axes)
            currDelta=currMax-currMin
            # a shift
            newMin=currMin+myCoord
            newMax=currMax+myCoord
            if (lAbs):
                # actually an absolute position
                currMean=(currMax+currMin)*0.5
                newMin=newMin-currMean; newMax=newMax-currMean
            self.setExtremes(newMin,newMax,axes=axes,lFree=lFree)
        else:
            print("Usrbin.move(): mixed arrays!")
            exit(1)

    def resize(self,newL,axis=3):
        if (axis!=3):
            print("Usrbin.resize(): For the time being, it is possible to re-size USRBINs only along the Z-axis!")
            exit(1)
        currMin,currMax=self.getExtremes(axes=axis); currNbins=self.getNbins(axes=axis)
        currDelta=currMax-currMin; currStep=currDelta*1.0/currNbins
        currMean=(currMax+currMin)*0.5
        # actually resize
        newNbins=newL*1.0/currDelta*currNbins
        if (newNbins%1<0.5):
            newNbins=int(newNbins)
        else:
            newNbins=int(newNbins)+1
        newMin=currMean-0.5*newNbins*currStep
        newMax=currMean+0.5*newNbins*currStep
        self.setExtremes(newMin,newMax,axes=axis)
        self.setNbins(nBins=newNbins,axes=axis)

    def assignTransformName(self,trasfName):
        if (len(trasfName)>0):
            self.TransfName=trasfName
    def retTransformName(self):
        return self.TransfName
    def isLinkedToTransform(self):
        return self.retTransformName() is not None
    
class TwoRegBasedScoring(Scoring):
    '''
    A very basic class for handling region-based scoring cards with two regions.
    '''
    def __init__(self,myName="",myComment="",scoType=""):
        Scoring.__init__(self,myName=myName,myComment=myComment,scoType=scoType)

    @staticmethod
    def fromBuf(myBuffer,newScoDet):
        return Scoring.fromBuf(myBuffer,newScoDet=newScoDet)
        
    def setRegName(self,whichReg,regName):
        if (whichReg==1 or whichReg==2):
            self.definition[0]=self.definition[0][:30+10*(whichReg-1)]+\
                               "%10s"%(regName.strip())+\
                               self.definition[0][40+10*(whichReg-1):]
        else:
            print("Usryield.setRegName(): which region do you choose?")
            exit(1)
    def retRegName(self,whichReg):
        if (whichReg==1 or whichReg==2):
            return self.definition[0][30+10*(whichReg-1):40+10*(whichReg-1)].strip()
        else:
            print("Usryield.retRegName(): which region do you choose?")
            exit(1)
    def regNameReplaceInDef(self,oldName,newName,nRegs=2):
        for iReg in range(nRegs):
            if (self.retRegName(iReg+1)==oldName):
                self.setRegName(iReg+1,newName)
                break

class Usryield(TwoRegBasedScoring):
    '''
    A very basic class for handling USRYIELD cards.
    '''
    def __init__(self,myName="",myComment=""):
        TwoRegBasedScoring.__init__(self,myName=myName,myComment=myComment,scoType="USRYIELD")

    @staticmethod
    def fromBuf(myBuffer):
        return TwoRegBasedScoring.fromBuf(myBuffer,newScoDet=Usryield())
        

class Usrbdx(TwoRegBasedScoring):
    '''
    A very basic class for handling USRBDX cards.
    '''
    def __init__(self,myName="",myComment=""):
        TwoRegBasedScoring.__init__(self,myName=myName,myComment=myComment,scoType="USRBDX")

    @staticmethod
    def fromBuf(myBuffer):
        return TwoRegBasedScoring.fromBuf(myBuffer,newScoDet=Usrbdx())
        
class Usrtrack(TwoRegBasedScoring):
    '''
    A very basic class for handling USRTRACK cards.
    '''
    def __init__(self,myName="",myComment=""):
        TwoRegBasedScoring.__init__(self,myName=myName,myComment=myComment,scoType="USRTRACK")

    @staticmethod
    def fromBuf(myBuffer):
        return TwoRegBasedScoring.fromBuf(myBuffer,newScoDet=Usrtrack())
        
    def setRegName(self,regName):
        TwoRegBasedScoring.setRegName(self,1,regName)
    def retRegName(self):
        return TwoRegBasedScoring.retRegName(self,1)
    def regNameReplaceInDef(self,oldName,newName):
        TwoRegBasedScoring.regNameReplaceInDef(self,oldName,newName,nRegs=1)

class Usrcoll(TwoRegBasedScoring):
    '''
    A very basic class for handling USRCOLL cards.
    '''
    def __init__(self,myName="",myComment=""):
        TwoRegBasedScoring.__init__(self,myName=myName,myComment=myComment,scoType="USRCOLL")

    @staticmethod
    def fromBuf(myBuffer):
        return TwoRegBasedScoring.fromBuf(myBuffer,newScoDet=Usrcoll())
        
    def setRegName(self,regName):
        TwoRegBasedScoring.setRegName(self,1,regName)
    def retRegName(self):
        return TwoRegBasedScoring.retRegName(self,1)
    def regNameReplaceInDef(self,oldName,newName):
        TwoRegBasedScoring.regNameReplaceInDef(self,oldName,newName,nRegs=1)

