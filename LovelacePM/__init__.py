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
from pytoolkit import *

#Define update module control functions
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
def update_message_enable():
    pdir=os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    lupfile=open('updatestat.lup', 'rb')
    statusdict=pickle.load(lupfile)
    lupfile.close()
    statusdict['message']=True
    lupfile=open('updatestat.lup', 'wb')
    pickle.dump(statusdict, lupfile)
    lupfile.close()
    os.chdir(pdir)
def update_message_suppress():
    pdir=os.getcwd()
    os.chdir(os.path.dirname(os.path.abspath(__file__)))
    lupfile=open('updatestat.lup', 'rb')
    statusdict=pickle.load(lupfile)
    lupfile.close()
    statusdict['message']=False
    lupfile=open('updatestat.lup', 'wb')
    pickle.dump(statusdict, lupfile)
    lupfile.close()
    os.chdir(pdir)

#Generate update registry file if none is found
tdy=LoveUpdate.getdate()
if not os.path.exists('updatestat.lup'):
    print('NOTICE: update registry file not found. Creating anew')
    lupfile=open('updatestat.lup', 'wb')
    statusdict={'lastupdate':tdy, 'updateme':True, 'message':True}
    _=LoveUpdate.update('LoveUpdate', '')
    pickle.dump(statusdict, lupfile)
    lupfile.close()

#Read update registry file
lupfile=open('updatestat.lup', 'rb')
statusdict=pickle.load(lupfile)
lupfile.close()

#Declare variables
hasupdated=False
hasupdated_lup=False
code_lup=False
code=False

#If an older version update file doesn't include message suppression/enabling variable, create it
if not 'message' in statusdict:
    update_message_enable() #activate message issuing when autoupdate is not performed
    statusdict['message']=True

#Update if activated
if statusdict['updateme']:
    code_lup, hasupdated_lup, tdy=LoveUpdate.update('LoveUpdate', statusdict['lastupdate'])
    if not code_lup:
        if hasupdated_lup:
            print('LoveUpdate has autoupdated. If you wish to deactivate the autoupdate feature, import LovelacePM in another python instance and use LovelacePM.updatecancel()')
        code, hasupdated, tdy=LoveUpdate.update('LovelacePM', statusdict['lastupdate'])
        code=code and code_lup
        if not code:
            if hasupdated:
                print('LovelacePM has autoupdated. If you wish to deactivate the autoupdate feature, import LovelacePM in another python instance and use LovelacePM.updatecancel()')
            statusdict['lastupdate']=tdy
            lupfile=open('updatestat.lup', 'wb')
            pickle.dump(statusdict, lupfile)
            lupfile.close()
        else:
            print('WARNING: error in autoupdate. Check internet connection if you wish to update LovelacePM. Moving on')
    else:
        print('WARNING: error in autoupdate. Check internet connection if you wish to update LoveUpdate. Moving on')
elif statusdict['message']: #warn if message flag is activated
    print('WARNING: LovelacePM autoupdate is deactivated. if you wish to activate autoupdate, import LovelacePM in another python instance and use LovelacePM.updateset()')
    print('To suppress this message, use LovelacePM.update_message_suppress(). To undo this action, use LovelacePM.update_message_enable()')


del code, code_lup, lupfile, tdy, hasupdated_lup, hasupdated
os.chdir(ordir)
del ordir
