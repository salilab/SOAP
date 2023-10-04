"""
SOAP env module, defining the computing environments.

Dependence (in both local host and sgeserver): numpy, scipy, bottleneck, modeller, mdt, h5py, pymc, matplotlib.

If the installation paths for these packages are not in the python search path (in both local host and sgeserver), please add the paths to $PYTHONPATH or sys.path.

You might need to run "module load modeller" to load modeller and MDT.

**Rename this file to env.py after proper configuration**

"""
from __future__ import print_function
import sys
import os
hostname=os.getenv('HOSTNAME')
print(hostname)
from modeller import *
if sys.version_info[0] >= 3:
    import queue
else:
    import Queue as queue
import threading


class SOAPenv(object):
    """This class defines the directory and other environment variables for all the other classes.
    These variables need to be defined::

        self.localHost=''#local server for storing results and tables, and running SOAP
        self.sgeServer=''#for submitting jobs on sge cluster, ssh key authencation needed for passwordless access both way (local host<->server)
        self.userName="" # username on local host and sge server, must use RSA Key for passwordless SSH Authentication
        self.localStatPotPath=''# the root path of statistical potential files in local host (could take TB of spaces)
        self.localDecoysPath='' # the directory for storing decoys on local host
        self.localProcessedPDBPaths=['',''] # paths to preprocessed pdb files on local host
        self.serverProcessedPDBPaths=['',''] # paths to processed pdb files on sge server
        self.localRawPDBPath='' # path to original pdb files on local host
        self.serverRawPDBPath=''# paths to processed pdb files on sge server
        self.serverLogPath='' # for storing logs when path name contains "\#"
        self.serverUserPath='' # the working directory on sge server
        self.serverHDF5Path='' # for storing hdf5 tables on sge server, temporary use only
        self.localInstallPath='' # path to SOAP module on localHost
        self.serverInstallPath='' # path to SOAP module on sge server
        self.serverScrathPath=''# for storing temporary results on sge cluster
        self.localRunCachePath='' # temporary run path on local host
        self.localResultDB='x.sql3'  # database file location on local host

    .. note::
        the statpot directory must contain:

        1. file rsn with a 5-6 digit number in it (run serial number), eg. 10000

        2. folder results for storing optimization results
        3. folder Decoys for storing decoys metadata used during training
        4. folder pdbpir for storing the filtered pir sequence files for PDB structures

    .. warning::
        The paths in lib/universal_script.sh also need to be configured correctly for SOAP to run on SGE cluster.

    .. warning::
        All paths defined in this file as well as the database file must exist before SOAP can be run.

    """


    def __init__(self):
        #define hostnames
        self.localHost=''#local server for storing results and tables, and running SOAP
        self.sgeServer=''#for submitting jobs on sge cluster, ssh key authencation needed for passwordless access both way (local host<->server)
        self.userName="" # username on local host and sge server, must use RSA Key for passwordless SSH Authentication
        self.localStatPotPath=''# the root path of statistical potential files in local host (could take TB of spaces)
        self.localDecoysPath='' # the directory for storing decoys on local host
        self.localProcessedPDBPaths=['',''] # paths to preprocessed pdb files on local host
        self.serverProcessedPDBPaths=['',''] # paths to processed pdb files on sge server
        self.localRawPDBPath='' # path to original pdb files on local host
        self.serverRawPDBPath=''# paths to processed pdb files on sge server
        self.serverLogPath='' # for storing logs when path name contains "\#"
        self.serverUserPath='' # the working directory on sge server
        self.serverHDF5Path='' # for storing hdf5 tables on sge server, temporary use only
        self.localInstallPath='' # path to SOAP module on localHost
        self.serverInstallPath='' # path to SOAP module on sge server
        self.serverScrathPath=''# for storing temporary results on sge cluster
        self.localRunCachePath='' # temporary run path on local host
        self.localResultDB='x.sql3'  # database file location on local host
        self.hostname=hostname
        if self.hostname==self.localHost:
            self.hostn=0#local
            self.basedir=self.localStatPotPath
            self.loopStructurePath=os.path.join(self.basedir, 'loop')
            self.ddbdir=self.localDecoysPath
            self.pdbdir=self.localProcessedPDBPaths
            self.opdbdir=self.localRawPDBPath
            self.decoysbasedir=os.path.join(self.basedir, 'Decoys')
            self.pdbpirdir=os.path.join(self.basedir, 'pdbpir')
            self.jobserver=self.userName+'@'+self.sgeServer
            self.scriptdir=self.localInstallPath
            self.logpath=os.path.join(self.basedir, 'scorelog')
            self.runs=self.localRunCachePath
            self.otherrun=True #True, there are other instances doing the same run or there aren't
            self.currentrundirlist=[]
            self.jobserverhdf5path=self.serverHDF5Path
            self.resultdbpath=self.localResultDB
            self.libdir=os.path.join(self.localInstallPath, 'lib')
            self.debug=True
        else:
            self.hostn=1
            self.loopStructurePath=os.path.join(self.serverUserPath, 'loop')
            self.basedir=os.getenv("RUNPATH")
            self.ddbdir=os.path.join(self.serverUserPath, 'decoys')
            self.pdbdir=self.serverProcessedPDBPaths #'/netapp/database/pdb/remediated/pdb/'#
            self.opdbdir=self.serverRawPDBPath
            self.decoysbasedir=None
            self.libdir=os.path.join(self.serverInstallPath, 'lib')
        self.hostn=0
        self.runpriority=0
        self.env=Environ()
        log.none()
        self.mcmc_temperature_ratios1=[3,3,2.5,2.5,2,2,1.5,1.5,1,1]
        self.acc_rate_ratio=1.0
        self.stepblockind=-1
        self.mcmc_pt=(1,10)
        self.queues=''
        self.numofthread=0
        self.queue=queue.Queue() #queue item format [func, [paras]]
        self.spdsscoreTasks={}

    def log(self,sdict):#obsolete
        logstrfile=open(self.basedir+'score.summary','a')
        logstrfile.write(sdict['string']+'\n')
        logstrfile.close()
        if os.path.isfile(self.basedir+'score.summary.pickle'):
            logcPicklefile=open(self.basedir+'score.summary.pickle','rb')
            loglist=cPickle.load(logcPicklefile)
            loglist.append(sdict)
            logcPicklefile.close()
        else:
            loglist=[]+[sdict]
        logcPicklefile=open(self.basedir+'score.summary.pickle','wb')
        cPickle.dump(loglist,logcPicklefile)
        logcPicklefile.close()

    def setup_threadpool(self):
        print("###setup thread pool")
        for i in range(self.numofthread):
            t = ThreadClass()
            t.setDaemon(True)
            t.start()
        print('###setup thread pool finished')

runenv=SOAPenv()

class ThreadClass(threading.Thread):
    # will take any jobs from queue, while the queue
    def __init__(self):
        threading.Thread.__init__(self)
        self.queue = runenv.queue

    def run(self):
        while True:
            inputlist=self.queue.get()
            #st=time.time()
            inputlist[0](*inputlist[1])
            #print("calculation time "+str(time.time()-st))
            self.queue.task_done()
        print("quit")

sys.path.extend([])
scriptname='universal_script.sh'
runtimetest=False
localtest=False
runtime='336:0:0' #336:0:0
scoreminnumber=500.0 #(try not to increase this number becasue random error on mdt when the number if large)
smoothingjobratio=2*10000000.0
optimizerjobtime=5.0
#alllabqueue=0
accessibilitytypes=['l']
mcverbose=1

from error import *
import re
import time
from modeller import optimizers
from modeller import automodel
from modeller.optimizers import actions
from modeller.automodel import *
from modeller.scripts import complete_pdb
from modeller.parallel import *
import mdt
import mdt.features
import pickle
if sys.version_info[0] >= 3:
    cPickle = pickle
else:
    import cPickle
import subprocess
from operator import itemgetter, attrgetter
import random
import gzip
import traceback
import gc
import h5py
import numpy as np
import numpy.linalg
import scipy.cluster
import scipy.interpolate
import scipy.optimize
import scipy.stats
from pymc import *
import copy
import bottleneck


#set up scratch dir and host specific modules
if runenv.hostn==0:
    import pdb
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pylab as plt
    import signal
    import shelve
    #import pypchip
    scratchdir='/scratch'
else:
    if not os.getenv('JOB_ID'):
        scratchdir='/scratch/'+runenv.userName+'/temp/'
        if not os.path.isdir(scratchdir):
            os.makedirs(scratchdir)
    else:
        import statvfs
        print(os.statvfs("/scratch"))
        dfs=os.statvfs("/scratch")
        if (dfs[statvfs.F_BAVAIL]*dfs[statvfs.F_BSIZE])>100000000:
            scratchdir=runenv.serverScrathPath+os.getenv('JOB_ID')+'/'+os.getenv('SGE_TASK_ID')+'/'
        else:
            scratchdir=runenv.serverScrathPath+os.getenv('JOB_ID')+'/'+os.getenv('SGE_TASK_ID')+'/'
        if not os.path.isdir(scratchdir):
            os.makedirs(scratchdir)

from utility import *
from jobSchedule import *
