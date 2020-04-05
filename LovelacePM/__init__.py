import os
import sys
import datetime
import pickle
import LoveUpdate
ordir=os.getcwd()
os.chdir(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.getcwd())
from utils import *
from wing import *
from body import *
from control import *
from aircraft import *
from xfoil_visc import *
from aerodynamic_output import *
from toolkit import *
tdy=LoveUpdate.getdate()
if not os.path.exists('updatestat.lup'):
    lupfile=open('updatestat.lup', 'wb')
    statusdict={'lastupdate':tdy, 'updateme':True}
    pickle.dump(statusdict, lupfile)
    lupfile.close()
lupfile=open('updatestat.lup', 'rb')
statusdict=pickle.load(lupfile)
lupfile.close()
if statusdict['updateme']:
    retval=LoveUpdate.update('LovelacePM', statusdict['lastupdate'])
    if type(retval)==tuple:
        code=retval[0]
        tdy=retval[1]
        if not code:
            print('LovelacePM has autoupdated. If you wish to deactivate the autoupdate feature, import LovelacePM in another python instance use LovelacePM.updatecancel()')
    elif type(retval)==bool:
        code=retval
    else:
        code=False
    if not code:
        statusdict['lastupdate']=tdy
        lupfile=open('updatestat.lup', 'wb')
        pickle.dump(statusdict, lupfile)
        lupfile.close()
    else:
        print('WARNING: error in autoupdate. Check internet connection if you wish to update LovelacePM. Moving on')
else:
    retval=True
    print('WARNING: LovelacePM autoupdate is deactivated. if you wish to activate autoupdate, import LovelacePM in another python instance use LovelacePM.updateset()')
def updateset():
    pdir=os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    lupfile=open('updatestat.lup', 'rb')
    statusdict=pickle.load(lupfile)
    lupfile.close()
    statusdict['updateme']=True
    lupfile=open('updatestat.lup', 'wb')
    pickle.dump(statusdict, lupfile)
    lupfile.close()
    os.chdir(pdir)
def updatecancel():
    pdir=os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    lupfile=open('updatestat.lup', 'rb')
    statusdict=pickle.load(lupfile)
    lupfile.close()
    statusdict['updateme']=False
    lupfile=open('updatestat.lup', 'wb')
    pickle.dump(statusdict, lupfile)
    lupfile.close()
    os.chdir(pdir)
del retval, tdy, lupfile
os.chdir(ordir)
del ordir
