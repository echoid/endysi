# Utility functions

import os
import time

def join(p1, p2):
    p3 = ''
    # '.../' and '/...'
    if p1.endswith('/') and p2.startswith('/'):
        p3 = p1 + p2[1:len(p2)]

    # '...' and '/...'
    elif not p1.endswith('/') and p2.startswith('/'):
        p3 = p1 + p2

    elif p1.endswith('/') and not p2.startswith('/'):
        p3 = p1 + p2

    # '...' and '...'
    elif not p1.endswith('/') and not p2.startswith('/'):
        p3 = p1 + '/' + p2

    return p3

def getNameFrom(path):
    (d, n) = os.path.split(path)
    (f, e) = os.path.splitext(n)
    return f

def getDirFrom(path):
    (d, n) = os.path.split(path)
    return d

def genDirNameFrom(caLevel):
    caStr = str(caLevel)
    return caStr.replace('.', 'p')

def genCaLevelFrom(dirName):
    caStr = dirName.replace('p', '.')
    return float(caStr)

def genModelFileName(dirName):
    return '_' + dirName + '.bngl'

def makeDirs(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return

def getEndDirFrom(path):
    dirs = path.split('/')
    while '' in dirs:
        dirs.remove('')

    return dirs[-1]


def genTimeString():
    return time.strftime('%m-%d-%Y_%H-%M-%S')

def genDateString():
    t1 = time.localtime()
    dateString = '%2d-%2d-%4d' % (t1.tm_mon, t1.tm_mday, t1.tm_year)
    return dateString.replace(' ', '0')

def removeEquilFileFrom(filelist):
    for f in filelist:
        if 'equil' in f:
            filelist.remove(f)

    return

def genStimLevelFrom(dirName):
    return float(dirName.replace('p', '.'))

def getConditionFrom(filename):
    name = getNameFrom(filename)
    s = name.split('_')
    return float(s[1])

def getStrConditionFrom(filename):
    name = getNameFrom(filename)
    s = name.split('_')
    return s[1]

def getIntConditionFrom(filename):
    name = getNameFrom(filename)
    s = name.split('_')
    return int(s[1])





