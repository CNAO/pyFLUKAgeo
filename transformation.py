# classes for managing FLUKA transformations
# A. Mereghetti, 2024/04/16
# python version: >= 3.8.10

import numpy as np
from FLUKA import GeoObject, echoFloats, assembleLine, MaxLineLength, LineHeader

class RotDefi(GeoObject):
    '''
    Class handling ROT-DEFI cards.
    Please have a look also to the Transformation() class: the ROT-DEFI
       class simply stores infos of a specific rotation/traslation;
    - echo in FREE format uses an empty space as field separator;
    - parsing in FREE format expects an empty space as field separator;
    '''

    def __init__(self,myAx=3,myPhi=0.0,myTh=0.0,myDD=np.zeros(3),myComment=""):
        # mainly for comments, name NOT used!
        GeoObject.__init__(self,myComment=myComment)
        self.axis=myAx           # FLUKA coding: 1=x, 2=y, 3=z
        self.phi=myPhi           # polar angle [degrees] (FLUKA units)
        self.theta=myTh          # azimuth angle [degrees] (FLUKA units)
        self.DD=myDD             # translation [cm] (x,y,z components)

    def echo(self,lFree=True,myID=0,myName="",\
             maxLen=MaxLineLength,header=LineHeader,lMultiLine=False):
        '''echo in FREE format uses an empty space as field separator'''
        if (myID<=0):
            print("...cannot dump a ROT-DEFI card without an index!")
            exit(1)
        if (len(myName)==0):
            print("...cannot dump a ROT-DEFI card without a name!")
            exit(1)
        if (myID<100):
            myIDecho=myID
            if (self.phi!=0.0 or self.theta!=0.0):
                myIDecho=myIDecho+self.axis*100
        else:
            myIDecho=myID*1000
            if (self.phi!=0.0 or self.theta!=0.0):
                myIDecho=myIDecho+self.axis
        myFloats=[myIDecho,self.phi,self.theta]+list(self.DD)
        myStrings=["ROT-DEFI  "]+echoFloats(myFloats,lFree=lFree)+[" %-10s"%(myName)]
        return GeoObject.echoComm(self)+ \
            assembleLine(myStrings,maxLen=maxLen,header=header,lMultiLine=lMultiLine)

    @staticmethod
    def fromBuf(myBuffer,lFree=True):
        '''parsing in FREE format expects an empty space as field separator'''
        newRotDefi=RotDefi()
        for tmpLine in myBuffer.splitlines():
            if (tmpLine.startswith("*")):
                newRotDefi.tailMe(tmpLine.strip())
            else:
                # parse data
                if (lFree):
                    data=tmpLine.split()
                    pID=data[1]
                    phi=data[2]
                    tht=data[3]
                    DD=data[4:7]
                    myName=data[7]
                else:
                    pID=tmpLine[10:20].strip()
                    phi=tmpLine[20:30].strip()
                    tht=tmpLine[30:40].strip()
                    DD=[tmpLine[40:50].strip(),tmpLine[50:60].strip(),tmpLine[60:70].strip()]
                    myName=tmpLine[70:80].strip()
                # store data
                pID=float(pID)
                if (pID>1000):
                    myID=int(pID/1000)
                    myAxis=int(pID%1000)
                else:
                    myAxis=int(pID/100)
                    myID=int(pID%100)
                if (myAxis>0):
                    newRotDefi.axis=myAxis
                newRotDefi.phi=float(phi)
                newRotDefi.theta=float(tht)
                newRotDefi.DD=np.array([float(myDD) for myDD in DD])
        return newRotDefi, myID, myName

class Transformation(GeoObject):
    '''
    This is the actual class to be used when dealing with ROT-DEFI cards.
    The class implements an array of ROT-DEFI cards - hence, at least one card.
    '''

    def __init__(self,myID=0,myName=""):
        GeoObject.__init__(self,myName=myName) # mainly for names, comments NOT used!
        self.RotDefis=[] # list of ROT-DEFI cards
        self.ID=myID
        
    def __len__(self):
        return len(self.RotDefis)

    def __iter__(self):
        self.current_index=0
        return self
    
    def __next__(self):
        '''
        from https://blog.finxter.com/python-__iter__-magic-method/
        '''
        if self.current_index < len(self):
            x = self.RotDefis[self.current_index]
            self.current_index += 1
            return x
        raise StopIteration

    # override methods of GeoObject, since comments are of the single
    #   ROT-DEFI cards
    def headMe(self,myString):
        self.RotDefis[0].headMe(myString)
    def tailMe(self,myString):
        self.RotDefis[0].tailMe(myString)
    def echoComm(self):
        myComm=""
        if (len(self.RotDefis)>0):
            myComm=self.RotDefis[0].echoComm()
        return myComm

    def AddRotDefi(self,newRotDefi,iAdd=-1,lMoveComment=True):
        if (iAdd==-1 or iAdd==len(self.RotDefis)):
            # append ROT-DEFI card to list of ROT-DEFIs
            self.RotDefis.append(newRotDefi)
        elif (iAdd==0):
            # head ROT-DEFI card to list of ROT-DEFIs
            if (lMoveComment):
                myComm=self.echoComm().strip()
                if (len(myComm)>0):
                    # take care of main comment
                    newRotDefi.headMe(myComm)
                    self.RotDefis[0].comment=""
            self.RotDefis=[newRotDefi]+self.RotDefis
        else:
            # insert ROT-DEFI card at specific position in the list of ROT-DEFIs
            self.RotDefis=self.RotDefis[0:iAdd]+[newRotDefi]+self.RotDefis[iAdd:]
    def AddRotDefis(self,newRotDefis,iAdd=-1,lMoveComment=True):
        for myRotDefi in newRotDefis:
            self.AddRotDefi(myRotDefi,iAdd=iAdd,lMoveComment=lMoveComment)

    @staticmethod
    def fromBuf(self,myBuffer,lFree=True):
        newTransf=Transformation()
        myBuf=""
        for tmpLine in myBuffer.splitlines():
            if (tmpLine.startswith("*")):
                if (len(myBuf)==0):
                    myBuf=tmpLine
                else:
                    myBuf=myBuf+"\n"+tmpLine
            else:
                tmpRotDefi=RotDefi()
                myID,myName=tmpRotDefi.fromBuf(myBuf)
                if (myID!=0):
                    if (newTransf.ID==0):
                        # first ROT-DEFI card of transformation
                        newTransf.ID=myID
                    elif(newTransf.ID!=myID):
                        print("...cannot add ROT-DEFI with ID=%d to tranformation ID=%d!"%(myID,newTransf.ID))
                        exit(1)
                if (len(myName)>0):
                    if (len(newTransf.echoName())==0):
                        # first ROT-DEFI card of transformation
                        newTransf.rename(myName,lNotify=False)
                    elif (newTransf.echoName()!=myName):
                        print("...cannot add ROT-DEFI named %s to tranformation named %s !"%(myName,newTransf.echoName()))
                        exit(1)
                newTransf.AddRotDefi(tmpRotDefi)
                myBuf=""
        return newTransf

    def echo(self,lFree=True):
        buf=""
        for myRotDefi in self.RotDefis:
            if (len(buf)>0):
                buf=buf+"\n"
            buf=buf+myRotDefi.echo(myID=self.ID,myName=self.name)
        return buf

