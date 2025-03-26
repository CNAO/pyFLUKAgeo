# general python material for FLUKA stuff
# A. Mereghetti, 2024/04/16
# python version: >= 3.8.10

StringPrec=2.0E-13 # prec in determining integer values
MaxLineLength=132
LineHeader="%13s"%("")
begLineGLOB="* \n"+"* "+"="*(MaxLineLength-2)
endLineGLOB="* "+"-"*(MaxLineLength-2)+"\n* "
maxGeoNameLen=8
# FREE format:
cut0sFREE="0"*4
expDigitsFREE=22
numDigitsFREE=expDigitsFREE-7
# FIXED format:
cut0sFIXED="0"*2
expDigitsFIXED=10
numDigitsFIXED=expDigitsFIXED-7
# definition handling
remSpaceAfters=["+","-","(",")","|"]
addSpaceAfters=["(",")","|"]

class GeoObject():
    '''
    A very basic object, implementing simply a name and a comment.
    This object is supposed to collect all the stuff concerning
      handling of names and modification of comments of instances
      of the other classes;
    '''

    def __init__(self,myName="",myComment=""):
        self.name=myName
        self.comment=myComment

    def rename(self,newName,lNotify=True):
        'change instance name'
        if (lNotify):
            self.tailMe("* NAME CHANGE: FROM %s TO %s"%(self.name,newName))
        self.name=newName

    def headMe(self,myString):
        'head a string to the comment'
        if (len(self.comment)>0):
            self.comment=myString+"\n"+self.comment
        else:
            self.comment=myString
            
    def tailMe(self,myString):
        'tail a string to the comment'
        if (len(self.comment)>0):
            self.comment=self.comment+"\n"+myString
        else:
            self.comment=myString
                
    def echoComm(self):
        myBuf=""
        if (len(self.comment)>0):
            myBuf=self.comment+"\n"
        return myBuf

    def echoName(self):
        return self.name
            
def echoFloats(myFloats,expDigits=None,numDigits=None,prec=StringPrec,cut0s=None,lFree=True):
    'function for echoing floats in a human-readable way'
    if (lFree):
        if (expDigits is None): expDigits=expDigitsFREE
        if (numDigits is None): numDigits=numDigitsFREE
        if (cut0s is None): cut0s=cut0sFREE
    else:
        if (expDigits is None): expDigits=expDigitsFIXED
        if (numDigits is None): numDigits=numDigitsFIXED
        if (cut0s is None): cut0s=cut0sFIXED
    numFmt="%% %d.%dE"%(expDigits,numDigits)
    if (isinstance(myFloats,list)):
        floats2convert=myFloats
    else:
        floats2convert=[myFloats]
    myStrs=[]
    for float2convert in floats2convert:
        if (abs(float2convert)%1<prec): # integer
            myStr="% .1f"%(float2convert)
        else: # actual floating point number
            myStr=numFmt%(float2convert)
            data=myStr.split("E")
            mantissa=data[0]; exponent=data[1]
            if (cut0s is not None): # truncate numbers if a series of "0" is found
                if (mantissa.find(cut0s)>-1):
                    mantissa=mantissa[:mantissa.index(cut0s)]
            if (exponent=="+00"):
                exponent=None
            elif (exponent=="+01"):
                exponent=None
                mantissa=mantissa.replace(".","")
                mantissa=mantissa[0:3]+"."+mantissa[3:]
            elif (exponent=="+02"):
                exponent=None
                mantissa=mantissa.replace(".","")
                mantissa=mantissa[0:4]+"."+mantissa[4:]
            elif (exponent=="-01"):
                exponent=None
                mantissa=mantissa.replace(".","")
                mantissa=mantissa[0]+"0."+mantissa[1:]
            elif (exponent=="-02"):
                exponent=None
                mantissa=mantissa.replace(".","")
                mantissa=mantissa[0]+"0.0"+mantissa[1:]
            if (mantissa=="-0.0"): mantissa=" 0.0"
            if (exponent is None):
                myStr=mantissa
            else:
                myStr="%sE%s"%(mantissa,exponent)
        if (not lFree):
            if (len(myStr)<expDigits):
                myStr=" "*(expDigits-len(myStr))+myStr
        myStrs.append(myStr)
    return myStrs

def assembleLine(myStrings,maxLen=MaxLineLength,header=LineHeader,lMultiLine=True):
    'function for assembling a FLUKA line made of myStrings respecting all FLUKA constraits'
    if (len(myStrings)==0):
        finStr=""
    else:
        finStr=myStrings[0]; curLen=len(finStr)
        for myStr in myStrings[1:]:
            if (curLen+len(myStr)+1>=maxLen):
                if (lMultiLine):
                    finStr=finStr+"\n"+header+myStr
                    curLen=len(header+myStr)
                else:
                    print("Cannot prepare a line longer than %d chars!"%(maxLen))
                    print("%s %s"%(finStr,myStr))
                    exit(1)
            else:
                finStr=finStr+" "+myStr
                curLen=curLen+len(myStr)+1
    return finStr

def cleanRegLine(myLine,remSpaceAfters=remSpaceAfters,addSpaceAfters=addSpaceAfters):
    'function for cleaning away useless empty spaces from a region definition line'
    oLine=myLine
    # - remove useless empty spaces
    for remSpaceAfter in remSpaceAfters:
        myRep="%1s "%(remSpaceAfter)
        while (myRep in oLine):
            oLine=oLine.replace(myRep,remSpaceAfter)
    # - re-introduce the necessary ones
    for addSpaceAfter in addSpaceAfters:
        myRep="%1s "%(addSpaceAfter)
        oLine=oLine.replace(addSpaceAfter,myRep)
    return oLine

def HighLightComment(myString,begLine=None,endLine=None):
    if (begLine is None): begLine=begLineGLOB
    if (endLine is None): endLine=endLineGLOB
    return begLine+"\n* "+myString+" \n"+endLine

def TailNameInt(myString,maxLen=maxGeoNameLen,nDigits=2,addChar="_",lDebug=True):
    if (nDigits is None or nDigits<0): nDigits=0
    nNameFmt=None
    if (len(myString)>maxLen-nDigits):
        newName=myString[0:maxLen-nDigits]
        if (lDebug): print("TailNameInt(): chopping name %s to len(%s)==%d"%(myString,newName,len(newName)))
    elif (len(myString)<maxLen-nDigits):
        newName=myString+addChar*(maxLen-nDigits-len(myString))
        if (lDebug): print("TailNameInt(): extending name %s to len(%s)==%d"%(myString,newName,len(newName)))
    else:
        newName=myString
    if (nDigits>0):
        newNameFmt=newName+"%0"+"%d"%(maxLen-len(newName))+"d"
    return newName, newNameFmt

def AddNameStr(myString,maxLen=maxGeoNameLen,addStr="",addChar="_",lTail=True,lDebug=True):
    addStrLen=len(addStr)
    if (len(myString)>maxLen-addStrLen):
        suffix=myString[0:maxLen-addStrLen]; myAddStr=""
        if (lDebug): print("AddNameStr(): chopping name %s to len(%s)==%d"%(myString,suffix,len(suffix)))
    elif (len(myString)<maxLen-addStrLen):
        suffix=myString; myAddStr=addChar*(maxLen-addStrLen-len(suffix))
        if (lDebug): print("AddNameStr(): extending name %s with len(%s)==%d"%(myString,myAddStr,len(myAddStr)))
    else:
        suffix=myString
    if (lTail):
        return suffix+myAddStr+addStr
    else:
        return addStr+myAddStr+suffix
