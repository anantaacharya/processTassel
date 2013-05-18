#written by ananta acharya,
#ananta@uga.edu
#has some common utilities like file I/O and other Sequence related
#commontools.py
#may change time to time

import os
def OKtoRead(inputfile):
    """Checks whether the input file is readable"""
    if not os.path.exists(inputfile):
        print "File does not Exist"
        return False
    elif not os.path.isfile(inputfile):
        print "Is not a file"
        return False
   
    return True
    
def OKtoWrite(outputfile):
    """Checks whether the input file is writable"""
    if not os.path.exists(os.path.dirname(outputfile)):
        print "Path Does Not Exist"
        return False
    elif os.path.isdir(outputfile):
        print "Is not a file but folder"
        return False
    elif os.path.exists(outputfile):
        print "File already exists, Can not Overwite or Merge"
        return False
    else:
        return True

def isNumber(s):
    try:
        float(s)
        return True
    except ValueError:
        return False

def checkRatio(s):
    try:
        s=float(s)
        if 0<=s<=1:
            return True
    except:
        return False    

def IUPAC(n1,n2):

    listt=[n1.upper(),n2.upper()]

    if "A" in listt and "T" in listt:
        return "W"
    elif "A" in listt and "C" in listt:
        return "M"
    elif "A" in listt and "G" in listt:
        return "R"
    elif "C" in listt and "T" in listt:
        return "Y"
    elif "C" in listt and "G" in listt:
        return "S"
    elif "G" in listt and "T" in listt:
        return "K"
    else:
        return "N"
    

     
    
def hamdist(str1, str2):
    assert len(str1) == len(str2)
    
    ne = operator.ne
    return sum(imap(ne, str1, str2))

def compare(str1, str2):
    mismatch=0
    if len(str1)<>len(str2):
        return False
    for i in range(len(str1)):
        
            if str1[i]<>str2[i]:
                mismatch+=1
                if mismatch==2:
                    return False
    else:
        return True
