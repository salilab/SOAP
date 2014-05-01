"""
   SOAP job control(local parallel or SGE parallel).

"""

from env import *


class task(object):
    def __init__(self,dir='',rdirname='',afterprocessing=None, preparing=None, targetfile=''):
        self.dir=dir
        self.rdirname=rdirname
        if '#' in self.rdirname:# '#' is not recoginized by sge as valid filename
            self.logpath=os.path.join(runenv.runlogpath, self.rdirname.replace('#',''))
        else:
            self.logpath=os.path.join(runenv.serverbasedir,rdirname,'output')
        self.runname='runme.sh' #should not be initilized with other values
        self.rtrys=0
        self.afterprocessing=afterprocessing
        self.k=0
        self.runstatusdict={}
        self.baton2diskfull=False
        self.unnoticederror=False
        self.started=False
        self.preparing=preparing
        self.justwait=0
        self.targetfile=targetfile
        self.queues=runenv.queues
        #if isinstance(afterprocessing,k2cvcluster) and alllabqueue==0:
        #    self.queues='-q lab.q'
        self.crashjobdict={}
        self.rsn=-9999
        self.runduration=-9999 # the time it takes for the run to finish after it is submitted
        self.finished=False
        self.waitinglist=[]
       
    def get_runf(self):
        rfh=open(self.dir+self.runname)
        rfs=rfh.read()
        rfh.close()
        rer=re.search('1-([0-9]+)',rfs)
        self.rfs=rfs
        self.numofruns=int(rer.group(1))
 
    def get_inputlist(self):
        ilfh=open(self.dir+'inputlist')
        il=ilfh.read()
        self.inputlist=il.split(',')
    
    def get_sn(self):
        if self.rsn==-9999:
            while os.path.isfile(runenv.basedir+'rsnlock'):
                time.sleep(1)
            print os.system('touch '+runenv.basedir+'rsnlock')
            rsnf=open(runenv.basedir+'rsn','r')
            rsn=int(rsnf.read())
            rsnf.close()
            rsnf=open(runenv.basedir+'rsn','w')
            rsn=rsn+1
            rsnf.write(str(rsn))
            rsnf.close()
            print os.system('rm '+runenv.basedir+'rsnlock')
            self.rsn=rsn
        
    def get_tasksn(self):
        #get a unique number for this task to be used as the unique identifier...
        if self.rsn<0:
            self.get_sn()
        rsn=self.rsn
        runname='r'+str(rsn)+'.sh'
        self.runname=runname
        os.system('mv '+self.dir+'runme.sh '+self.dir+self.runname)
        
    def run_task_cluster(self):
        self.get_tasksn()
        self.get_runf()
        self.get_inputlist()
        self.submit_task()
        self.started=True
        self.monitor_task()

    def reload_existingruns(self):
        self.started=True
        self.justwait=0
        self.waittinglist=[]
        os.chdir(self.dir)
        afl=os.listdir('./')
        runname=[f for f in afl if f.endswith('.sh.gz') and f.startswith('r')][0]
        self.runname=runname[:-3]
        print os.system('gunzip '+runname)
        self.rsn=int(self.runname[1:-3])
        self.get_runf()
        print os.system('gzip '+self.runname)
        self.recheck_runstatus()
        self.afterprocessing.nors=self.numofruns
     
    def start_task_cluster(self):
        self.started=True
        if os.path.isfile(self.targetfile):
            print 'Already executed by some other instances :)'
            return 1
        elif os.path.isdir(self.dir):
            if runenv.otherrun or (self.dir in runenv.currentrundirlist):
                print self.dir+' exist, meaning some other guys are generating what you want, just be patient and wait:)'
                self.waitinglist=[self.targetfile]
                self.justwait=1
                return 0
            else:
                print "reloading existing runs"
                self.reload_existingruns()
                runenv.currentrundirlist.append(self.dir)
                return 0
        self.waitinglist=self.preparing.prepare_task_input()
        runenv.currentrundirlist.append(self.dir)
        if isinstance(self.waitinglist,list) and len(self.waitinglist)>0:
            print 'waiting for others to calc this:'
            print self.waitinglist
            self.justwait=2
            return 0
        self.get_tasksn()
        self.get_runf()
        self.get_inputlist()
        self.submit_task()
        if '#' in self.rdirname:# '#' is not recoginized by sge as valid filename
            self.logpath=os.path.join(runenv.runlogpath, self.rdirname.replace('#',''))
        else:
            self.logpath=os.path.join(runenv.serverbasedir,self.rdirname,'output')         
        return 0

    def run_task_local(self,command,inputlist=None):
        self.get_runf()
        self.get_inputlist()
        if inputlist:
            self.inputlist=inputlist
        for item in self.inputlist:
            command(item)
   
    def submit_task(self):
        try:
            print os.system('ssh '+runenv.jobserver+' \' rm -r '+self.rdirname+' ;mkdir '+self.rdirname+'; mkdir '+self.rdirname+'/output''\'') 
            os.chdir(self.dir)
            print os.system('gzip -r -1 *')
            rc1=1
            while rc1:
                fl=os.listdir('./')
                hfl=[f for f in fl if f.endswith('hdf5')]
                if len(hfl)>0 and os.path.getsize(hfl[0])>1000000000:
                    hf=hfl[0]
                    ssh = subprocess.Popen(["ssh", "%s" % runenv.jobserver, 'ls '+runenv.jobserverhdf5path],
                       shell=False,
                       stdout=subprocess.PIPE)
                    result = ssh.communicate()#ssh.stdout.readlines()
                    if not hf in result[0]:
                        print os.system('scp '+hf+' '+runenv.jobserver+':'+runenv.jobserverhdf5path) 
                    print os.system('mv '+hf+' ../temp')
                    rc1=os.system('scp -r * '+runenv.jobserver+':~/'+self.rdirname+'/')
                    print os.system('mv ../temp'+' ./'+hf)
                    print os.system('ssh '+runenv.jobserver+' \' ln -s '+runenv.jobserverhdf5path+hf+' ~/'+self.rdirname+'/'+hf+'\'')
                else:
                    rc1=os.system('scp -r * '+runenv.jobserver+':~/'+self.rdirname+'/')
                print os.system('scp '+runenv.scriptpath+' '+runenv.jobserver+':~/'+self.rdirname+'/')
                 
                time.sleep(0.1)
            rc2=os.system('ssh '+runenv.jobserver+' \'cd '+self.rdirname+';gunzip *gz;qsub '+self.queues+' ./'+self.runname+'\'')#
            self.runstatus=np.ones(self.numofruns,dtype=np.int8) # 0 finished, 1 running,  2  failed, 3 marked as finished
            self.runtime=np.zeros(self.numofruns)
            self.runduration=time.time()
            for jn in range(0,self.numofruns):
                self.runstatusdict[jn]={'errorcodelist':[],'starttime':0}
            if runtimetest:
                self.runtime=np.ones(self.numofruns,dtype=[('runtime','i4'),('cn','f4'),('r2n','f4'),('host','a10')])
        except Exception,e:
            traceback.print_exc()
            time.sleep(60)
        
    def start_sr(self,jnl):
        for jn in jnl:
            #pdb.set_trace()
            self.runstatus[jn]=1
            jn=str(jn+1)
            runmdt2=open(self.dir+self.runname[:-3]+jn+'.sh','w')
            runmdt2s=self.rfs.replace(str(1)+'-'+str(self.numofruns),jn)
            runmdt2.write(runmdt2s)
            runmdt2.flush()
            os.chdir(self.dir)
            rc1=os.system('scp '+self.runname[:-3]+jn+'.sh'+' '+runenv.jobserver+':~/'+self.rdirname+'/')
            print os.system('rm '+self.runname[:-3]+jn+'.sh')
            #rc2=os.system('scp -r '+self.inputlist[int(jn)-1]+' '+runenv.jobserver+':~/'+self.rdirname+'/')
            rc3=os.system('ssh '+runenv.jobserver+' '+'\'qsub -p 0 '+self.queues+' ./'+self.rdirname+'/'+self.runname[:-3]+jn+'.sh\'')
            if rc3:
                print rc3
                raise NetworkError('Can not restart runs, possibly due to network problems')
            self.runstatusdict[int(jn)-1]['starttime']=0
            self.rtrys=self.rtrys+1
    
    def check_start_failed_runs(self,jnl):
        self.finishedrn=(self.runstatus==0).sum()
        self.activern=(self.runstatus==1).sum()
        self.waitingrn=(self.runstatus==4).sum()
        if len(jnl)==0:
            return 0
        if self.finishedrn==0 and self.activern==0 and len(jnl)>min(20,0.1*self.numofruns):
            raise Bugs('Bugs in code, all runs failed.')
        else:
            for jn in jnl:
                self.runstatusdict[jn]['errorcodelist'].append(self.check_runlog(jn))
            self.check_failed_runs(jnl)
            if len(jnl)>0.5*(self.numofruns-self.waitingrn):
                print "50% of the runs failed with error codes (): "+str(jnl)
                raise Exception('50% of the runs failed with error codes')
                self.errorjnl=jnl
            print "starting runs"
            self.start_sr(jnl)
            #repeated error for a single run should be noticed...
            
    def check_failed_runs(self,jnl):
        fcl=[]
        fkl=[]
        ct=True
        cn=0
        for key in self.runstatusdict:
            if len(self.runstatusdict[key]['errorcodelist'])>0 and self.runstatus[key]>0:
                fcl.append(self.runstatusdict[key]['errorcodelist'])
                fkl.append(key)
                if ct:
                    cn=cn+1
            else:
                ct=False
        if self.baton2diskfull:
            print "Disk quota exceeded, please clean disk before continue!!!!!!!!!!!!!!"
            pdb.set_trace()
            return 0
        if cn>max(20,0.1*self.numofruns):
            raise Bugs('First 20 or 10% has all failed, probably bugs in code')
        if len(squeeze_list(fcl))==1 and len(fcl[0])>1 and (len(fcl[0])**2*len(fcl)>50):
            raise Bugs('Bugs for some of the runs '+str(jnl))
            
    def check_runlog(self,jn):
        jn=str(jn+1)
        print 'checking run number '+jn
        proc1=subprocess.Popen('ssh '+runenv.jobserver+' '+'\'cd '+self.logpath+'; ls -c\'',shell=True,stdout=subprocess.PIPE)
        fc=proc1.communicate()[0]
        fl=fc.split('\n')
        fl2=[item for item in fl if item.endswith('.'+jn)]
        if len(fl2)==0:
            print 'no output log'
            return 'No output log'
        else:
            f=fl2[0]
        proc1=subprocess.Popen('ssh '+runenv.jobserver+' '+'\'cd '+self.logpath+'; cat '+f+' \'',shell=True,stdout=subprocess.PIPE)
        fc=proc1.communicate()[0]
        self.baton2diskfull=False
        if re.search('fork failure',fc) or re.search('Dynamic memory allocation failed',fc):
            print "memory problem with the run"
            ecode="memory problem with the run"
        elif re.search('Bugs',fc):
            print "Bugs"
            ecode='Bugs'
            pdb.set_trace()
        elif re.search('EOFError',fc):
            #rer=re.search('\\n.*\\nEOFError.*\\n')
            print "EOFError when reading "+fc
            ecode='EOFError'
        elif re.search('No space left on device',fc):
            print "The node or baton2 has no space left"
            ecode='No space'
        elif re.search('Killed',fc):
            print "Killed?????"
            ecode='Killed'
        elif re.search('Segmentation fault',fc):
            print "Segmentation fault - code problem"
            ecode='Segmentation fault'
        elif re.search('Error',fc):
            print "general Error"
            ecode='Errors'
        elif re.search('Aborted',fc):
            print "Aborted"
            ecode='Aborted'
        elif re.search('MemoryError',fc):
            print "MemoryError"
            ecode='MemoryError'
        elif re.search('Disk quota exceeded',fc):
            print "Disk quota exceeded, please clean disk before continue!!!!!!!!!!!!!!"
            ecode='Disk quota exceeded'
            self.baton2diskfull=True
        else:
            print "Unknown problem"
            ecode='Unknown'
        return ecode
        #_mdt.MDTError: unable to open file - no hdf5 file or file corrupt
            
    def copy_result(self,djl):
        #the limit of command length is 131072, so big lists need to be separated into small ones
        dl=len(djl)
        ep=0
        while ep<dl:
            if ep+400<dl:
                self.copy_result_400(djl[ep:ep+400])
                ep=ep+400
            else:
                self.copy_result_400(djl[ep:dl])
                ep=dl

    def copy_result_400(self,djl):
        #might worth check whether the file sizes agress after copy... in the future
        if len(djl)==0:
            return 0
        filelist=''
        rfl=[]
        rd='~/'+self.rdirname+'/'
        for (fjn, fjf) in djl:
            if self.runstatus[int(fjn)-1]>0:
                filelist=filelist+rd+fjf+' '
                rfl.append(fjf)
        if filelist=='':
            return 0
        mt=10
        if os.system('scp -r '+runenv.jobserver+':\"'+filelist+'\" '+self.dir):
            print 'file copy failure'
            efl,nefl=check_remote_file_exist(runenv.jobserver,rd,rfl)
            if len(nefl)>0:
                ndjl=[]
                for item in djl:
                    if item[1] in efl:
                        ndjl.append(item)
                    else:
                        self.runstatus[int(item[0])-1]=2
                djl=ndjl
                filelist=''
                for (fjn, fjf) in djl:
                    if self.runstatus[int(fjn)-1]>0:
                        filelist=filelist+rd+fjf+' '
            if len(efl)>0:
                ntry=10
                while os.system('scp -r '+runenv.jobserver+':\"'+filelist+'\" '+self.dir) and ntry>0:
                    ntry=ntry-1
                if ntry==0:
                    raise Exception('Can not copy files, disk full?'+filelist)
        tr=os.system('for i in *.tar.gz; do tar -xvzf $i; done')
        if tr:
            tr=os.system('for i in *.tar.gz; do tar -xvzf $i; done')
            if not tr:
                print os.system('rm *.tar.gz')
            else:
                raise Exception('can not untar file ')
        else:
            print os.system('rm *.tar.gz')
        for (fjn, fjf) in djl:
            if self.runstatus[int(fjn)-1]>0:
                self.runstatus[int(fjn)-1]=0
                print os.system('rm '+os.path.join(self.logpath,'*.'+fjn))
        print os.system('ssh '+runenv.jobserver+' '+'\'cd '+self.rdirname+'; rm '+filelist+'\'')
                    
    def get_runtimestats(self):
        #if not (runtimetest and os.path.isfile(self.dir+'pdbs1')):
        #    return 0
        proc1=subprocess.Popen('ssh '+runenv.jobserver+' '+'\'cd '+self.logpath+'; ls \'',shell=True,stdout=subprocess.PIPE)
        runstats=proc1.communicate()[0]
        rsl=runstats.split('\n')
        for runstatus in rsl:
            if runstatus.startswith('r'):
                continue
            rl=runstatus.split('.')
            print rl
            if len(rl)<3:
                continue
            fjn=rl[0]
            self.runtime[int(fjn)-1]['runtime']=int(rl[1])
            self.runtime[int(fjn)-1]['cn']=pir(path=self.dir+'pdbs'+fjn).count_pir()
            self.runtime[int(fjn)-1]['r2n']=pir(path=self.dir+'pdbs'+fjn).count_residue2()
            self.runtime[int(fjn)-1]['host']=rl[2]
     
    def get_wait_time(self):
        #not implemented yet
        self.wait_time=40
        return 40

    def process_result(self):
        self.runduration=time.time()-self.runduration
        if self.justwait==1:
            self.finished=True
            return 1
        print 'run finished:)'
        #try:
        #    self.get_runtimestats()
        #except Exception,e:
        #    traceback.print_exc()
        print os.system('ssh '+runenv.jobserver+' \'rm -r '+self.rdirname+' \'')
        if not self.logpath in ['/netapp/sali/gqdong/output','/netapp/sali/gqdong/output/']:
            print os.system('ssh '+runenv.jobserver+' \'rm -r '+self.logpath+' \'')
        #if runtimetest:
        #    np.save("../runtime"+str(self.rsn),self.runtime)
        if self.afterprocessing:
            #try:
            rv=self.afterprocessing.afterprocessing()
            #except Exception, e:
            #    print "error while processing the result of this run "
            #    traceback.print_exc()
            #    print "Please use the break point to look for whatever it wrong and continue without affecting other runs"
            #    pdb.set_trace()
            #    print "saving the failed object..."
            #    fh=open('run_failed_bugs.pickle','wb')
            #    cPickle.dump(self.afterprocessing,fh)
            #    fh.close()
            #    return 1
        else:
            rv=1
        self.finished=True
        return rv
    
    def monitor_task(self):
        time.sleep(60)
        k=0
        while True:
            if self.monitor():
                break
            time.sleep(60)
        
    def monitor(self):
        try:
            if self.finished:
                return 1
            if not self.started:
                print "##########monitoring:########## "+self.rdirname                
                return self.start_task_cluster()
            print "##########monitoring:########## "+self.rdirname
            if self.justwait:
                #print '############ Waiting for others to finish '
                nwl=[]
                for f in self.waitinglist:
                    if not os.path.isfile(f):
                        nwl.append(f)
                        #print 'WAITING FOR:  '+f
                self.waitinglist=nwl
                if len(nwl)==0:
                    return self.process_result()
                else:
                    print 'waiting for others'
                    return 0            
            if self.unnoticederror:
                self.process_keyboard_interrupt(self.error)
                print 'continue looping'
                return 0
            os.chdir(self.dir)
            proc1=subprocess.Popen('ssh '+runenv.jobserver+' '+'\' qstat | grep '+self.runname[:-3]+'\'',shell=True,stdout=subprocess.PIPE)
            runstatus=proc1.communicate()[0]
            if proc1.returncode==255:
                raise NetworkError('Network problem')
            proc2=subprocess.Popen('ssh '+runenv.jobserver+' '+'\' ls ~/'+self.rdirname+'/job*; rm ~/'+self.rdirname+'/jobfailed* \'',shell=True,stdout=subprocess.PIPE)
            jobstatlist=proc2.communicate()[0]
            if proc1.returncode==255:
                raise NetworkError('Network problem')
            print "obtaining runstats from server finished"
            print "***Don't interrupt***"
            self.runstatusstr=runstatus
            self.runstatus[self.runstatus>0]=2 #set status to failed
            if runstatus:
                self.analyze_runstatus(runstatus)
            fjl=re.findall('jobfailed\.([0-9]+)',jobstatlist)
            self.runstatus[[int(item)-1 for item in fjl if self.runstatus[int(item)-1]>0]]=2
            djl=re.findall('jobfinished\.([0-9]+)\.(.+)',jobstatlist)
            self.copy_result(djl)
            self.check_start_failed_runs(np.nonzero(self.runstatus==2)[0])
            self.k=self.k+1
            print "*********************"
            print 'Num of finished runs: '+str(self.finishedrn)+'/'+str(self.numofruns)+' active:'+str(self.activern) +' @iter#'+str(self.k)
            if (self.runstatus==0).sum()==self.numofruns:
                return self.process_result()
            elif self.rtrys >= 2*self.numofruns:
                raise Exception('jobfailed because of too many trys, check input files or scripts for errors')
        except NetworkError, e:
            traceback.print_exc()
            time.sleep(60)
        except Bugs, e:
            traceback.print_exc()
            #pdb.set_trace()
            #self.qdel_jobs()
            self.process_keyboard_interrupt(e)
        except KeyboardInterrupt,e:
            self.process_keyboard_interrupt(e)
        except Exception, e:
            print e
            traceback.print_exc()
            print self.dir
            self.process_keyboard_interrupt(e)
        return 0

    def recheck_k2cv_runstatus(self):
        os.chdir(self.dir)
        fl=os.listdir('./')
        fnl=[int(f.split('.')[0])-1 for f in fl if f.endswith('optimizer.pickle')]
        self.runstatus=np.ones(self.numofruns)        
        self.runstatus[fnl]=0
        
    def recheck_runstatus(self):
        os.chdir(self.dir)
        fl=os.listdir('./')
        fl=[f for f in fl if not (f.endswith('gz') or f.endswith('sh'))]
        #every input are gzipped ...
        fnl=[int(f.split('.')[0])-1 for f in fl ]
        self.runstatus=np.ones(self.numofruns)        
        self.runstatus[fnl]=0
        for jn in range(0,self.numofruns):
            self.runstatusdict[jn]={'errorcodelist':[],'starttime':0}   
        
    def reload_k2cvruns(self,rsn,numofruns=200):
        self.rsn=rsn
        self.dir=runenv.runs+str(rsn)+'/'
        self.rdirname=str(rsn)
        self.numofruns=numofruns
        self.runname='r'+str(rsn)+'.sh'
        rfh=open(self.dir+self.runname)
        rfs=rfh.read()
        rfh.close()
        self.rfs=rfs
        self.started=True
        self.unnoticederror=False
        for jn in range(0,self.numofruns):
            self.runstatusdict[jn]={'errorcodelist':[],'starttime':0}
        self.inputlist=['runme.py' for i in range(200)]
        self.recheck_k2cv_runstatus()
        self.monitor()

    def check_long_runs(self,runstatus): # will not be used
        runtime=self.runtime['runtime']
        runtime=runtime[runtime>0]
        runtimemean=runtime.mean()
        runtimestd=runtime.std()
        runtimelimit=min(runtimemean+5*runtimestd,runtimemean*2)
        rsl=self.decode_runstatus(runstatus)
        currenttime=time.time()
        nodechecklist=[]
        for item in rsl:
            if item[4] in ['r','Rr','t']:
                if currenttime-self.runstatusdict[int(item[-1])-1]['starttime']>rumtimelimit:
                    nodechecklist.append(item[-3])
        self.delete_runs_on_crashed_nodes(nodechecklist)
        
    def delete_runs_on_crashed_nodes(self):
        try:
            proc1=subprocess.Popen('ssh '+runenv.jobserver+' '+'\' qstat -qs u \'',shell=True,stdout=subprocess.PIPE)
            runstatus=proc1.communicate()[0]
            rsl=self.decode_runstatus(runstatus)
            autl=[]
            for item in rsl:
                if item[4] in ['r','Rr','t','dr']:
                    crashtask=item[0]+' -t '+item[-1]
                    if crashtask in self.crashjobdict:
                        self.crashjobdict[crashtask]+=1
                    else:
                        self.crashjobdict[crashtask]=1
            for key in self.crashjobdict.keys():
                if self.crashjobdict[key]>5:
                    autl.append(key)
                    del self.crashjobdict[key]
            print 'The following node crashed >5: '+str(autl)
            if len(autl)>0:
                print os.system('ssh '+runenv.jobserver+' '+'\' qdel -f '+', '.join(autl)+' \'')
        except:
            return 
        
    def get_single_node_status(self,node):
        proc1=subprocess.Popen('ssh '+runenv.jobserver+' '+'\' qstat -F -q '+node+'\'',shell=True,stdout=subprocess.PIPE)
        nodestatus=proc1.communicate()[0]
        
    def decode_runstatus(self,runstatus):
        runstatuslist=runstatus.split('\n')[:-1]
        rsl=[]
        for item in runstatuslist:
            if re.search('@', item):
                rsl.append(item.split())
        return rsl

    def analyze_runstatus(self,runstatus):
        runstatuslist=runstatus.split('\n')[:-1]
        for i in range(0,len(runstatuslist)):
            rstatus=runstatuslist[i].split()[4].strip()
            rl=runstatuslist[i].split()[-1].split(',')
            rl=[ri.split(':')[0] for ri in rl]
            rnl=[]
            for sr in rl:
                if re.search('([0-9]+)-([0-9]+)',sr):
                    rer=re.search('([0-9]+)-([0-9]+)',sr)
                    rnl=rnl+range(int(rer.group(1))-1,int(rer.group(2)))
                else:
                    rnl.append(int(sr)-1)
            for trn in rnl:
                if self.runstatus[trn]>0:
                    if rstatus in ['r','Rr','t','Rt','hr']:
                        self.runstatus[trn]=1
                        if not self.runstatusdict[trn]['starttime']:
                            self.runstatusdict[trn]['starttime']=time.time()
                    elif rstatus in ['qw','hqw']:
                        self.runstatus[trn]=4# waiting in the queue
                    elif rstatus.startswith('d'):
                        self.runstatus[trn]=2
                    else:
                        print 'can not decide the status of the run'
                        pdb.set_trace()

    def qdel_jobs(self,delete=False, qw=False):
        proc1=subprocess.Popen('ssh '+runenv.jobserver+' '+'qstat | grep '+self.runname[:-3]+' | awk \'{split($0,c,\" \");print  c[1]}\' ',shell=True,stdout=subprocess.PIPE)
        runnums=proc1.communicate()[0]
        runnums=set(runnums.split('\n'))
        for runnum in runnums:
            if qw and len(runnum)>5:
                print os.system('ssh '+runenv.jobserver+' '+'qdel -f '+runnum+' &')
            elif len(runnum)>5:
                print os.system('ssh '+runenv.jobserver+' '+'qdel -f '+runnum+' &')
        if delete:
            print os.system('ssh '+runenv.jobserver+' \'rm -r '+self.rdirname+' \'')
            print os.system('rm -r '+self.dir)

    def process_keyboard_interrupt(self,e, pos='task-monitor'):
        print '@@@@@@@@ '+pos
        print ""
        self.error=e
        TIMEOUT = 120 # number of seconds your want for timeout
        def interrupted(signum, frame):
            print "waiting more than 120s, continue looping......"
            raise Exception('timeout')
        signal.signal(signal.SIGALRM, interrupted)
        #def input():
        signal.alarm(TIMEOUT)
        try:
            print 'You have 120s to interrupt the code to handle the error, otherwise code will continue looping \n 0-quit and delete;\n 1-quit this job;\n 2-re-raise;\n 3-quit whole script and delete;\n 4-enter debug mode;\n 5-continue\n 6-restart runs with new sp.py\n 7-re-parepare and restart job \n 8-ignore erros \n '
            s = raw_input()
        except:
            return 0
        signal.alarm(0)
        self.unnoticederror=False
        if s=='0':
            self.qdel_jobs(delete=True)
            raise Bugs('Keyboard interruption')
        elif s=='1':
            self.qdel_jobs()
            raise Bugs('Keyboard interruption')
        elif s=='2':
            raise e
        elif s=='3':
            if pos=='tasklist-monitorend':
                self.qdel_jobs(delete=True)
                sys.exit(str(e))
            else:
                raise FatalError(str(e))
        elif s=='4':
            pdb.set_trace()
            self.started=True
            return 0
        elif s=='5':
            pass
        elif s=='6':
            self.submit_task()
        elif s=='7':
            self.start_task_cluster()            
        elif s=='8':
            self.error=None
            self.start_sr(self.errorjnl)
        else:
            self.unnoticederror=True
        
    def parallel_local(self,command,inputlist=[],nors=8):
        listoflist=[]
        for i in range(0,nors):
            list=[]
            listoflist.append(list)
            k=0
        for item in inputlist:
            listoflist[k%nors].append(item)
            k=k+1
        pids=[]
        for i in range(0,len(listoflist)):
            time.sleep(0.1)
            pid=os.fork()
            if pid == 0:
                print 'running child#'+str(i)
                for item in listoflist[i]:
                    try:
                        print 'child#'+str(i)+'start running'
                        command(item)
                    except Exception,e:
                        print e
                sys.exit()
                break
            elif pid > 0:
                pids.append(pid)
                if i==(len(listoflist)-1):
                    for j in range(0,len(pids)):
                        print 'waiting '+str(j)
                        os.wait()
                else:
                    continue

    def parallel_local_withreturn(self,command,inputlist=[],nors=8):
        if localtest or len(inputlist)==1:
            return self.serial_local(command,inputlist)
        listoflist=[]
        pipelistoflist=[]
        for i in range(0,nors):
            listoflist.append([])
            pipelistoflist.append([])
            k=0
        pipelist=[]
        for item in inputlist:
            r,w=os.pipe()
            pipelist.append([r,w])
            pipelistoflist[k%nors].append([r,w])
            listoflist[k%nors].append(item)
            k=k+1
        rl=[]
        pids=[]
        for i in range(0,len(listoflist)):
            time.sleep(0.1)
            pid=os.fork()
            if pid == 0:
                print 'running child#'+str(i)
                k=-1
                for item in listoflist[i]:
                    #try:
                    if 1:
                        k=k+1
                        p=pipelistoflist[i][k]
                        os.close(p[0])
                        w=os.fdopen(p[1],'w')
                        print 'child#'+str(i)+'start running'
                        result=command(item)
                        #print 'resultresultresultresult'
                        #print result
                        w.write(cPickle.dumps(result))
                        w.close()
                    #except Exception,e:
                    #    print e
                sys.exit()
                break
            elif pid > 0:
                pids.append(pid)
                if i==(len(listoflist)-1):
                    for j in range(0,len(pids)):
                        print 'waiting '+str(j)
                        os.wait()
                    for p in pipelist:
                        os.close(p[1])
                        r=os.fdopen(p[0])
                        rc=r.read()
                        print p
                        print 'result'+rc
                        rl.append(cPickle.loads(rc))
                        r.close()
                    return rl
                else:
                    continue

    def serial_local(self,command,inputlist=[]):
        outputlist=[]
        for input in inputlist:
            outputlist.append(command(input))
        return outputlist

class localtask(task):
    def __init__(self,func=None, inputlist=[]):
        self.func=func
        self.inputlist=inputlist
        self.finished=False
        
    def monitor(self):
        if self.finished:
            return 1
        tr=self.func(*self.inputlist)
        self.finished=True 
        return tr
        
    def qdel_jobs(self,delete=False):
        pass

class taskchain(task):
    #this calss will chain a list of tasks, and excute the chained functions one by one, including local runs
    #the inididual function should either excute some code, return some result as the input to the next function
        #or the function should result a task object, which can be monitored, and the monitor() return 0 when not success
    def __init__(self, chains=[]):#, sharedinput=[], initialinput=[]
        
        self.taskchains=chains
        #self.currenttask=[]
        self.currentpos=0
        #self.currentresult=initialinput
        #self.sharedinput=sharedinput
        self.started=False
    
    def monitor(self):
        self.currentresult=self.taskchains[self.currentpos].monitor()
        while self.currentresult!=0 and self.currentpos<(len(self.taskchains)-1):
            self.taskchains[self.currentpos]=0
            self.currentpos+=1
            self.currentresult=self.taskchains[self.currentpos].monitor()
        return self.currentresult #all finished

    def qdel_jobs(self,delete=False):
        if self.currentpos<len(self.taskchains):
            self.taskchains[self.currentpos].qdel_jobs(delete)

class tasklist(task): 
    def __init__(self,tasklist=[],afterprocessing=None,other=None,reload_rundir='',reload_dirlist=[]):
        if reload_rundir:
            fl=os.listdir(reload_rundir)
            fl=[os.path.join(reload_rundir,f) for f in fl if os.path.isdir(os.path.join(reload_rundir,f))]
            reload_dirlist=fl
        if len(reload_dirlist)>0:
            tasklist=self.reload_tasks(reload_dirlist)
        self.tasklist=tasklist
        self.afterprocessing=afterprocessing
        self.other=other
        self.resultlist=[0 for item in tasklist]
        for i in range(0,len(tasklist)):
            if tasklist[i]==0:
                self.resultlist[i]=123456789
        self.bugs=False
        self.buginfo=''
        self.unnoticederror=False
        self.crashjobdict={}
        
    def reload_tasks(self,dl):
        #only works for tasks with pickle saved...
        import gzip
        tl=[]
        for d in dl:
            ro=pickle.load(gzip.open(os.path.join(d,'input.pickle.gz')))
            ro.task.reload_existingruns()
            tl.append(ro.task)
        return tl
        
    def monitor(self):
        tli=range(0,len(self.tasklist))
        tli.reverse()
        for k in tli:
            if self.resultlist[k]:
                continue
            try:
                self.resultlist[k]=self.tasklist[k].monitor()
            except Bugs,e:
                self.resultlist[k]=987654321
                self.bugs=True
                self.buginfo=e
            except FatalError,e:
                raise e
            except:
                try:
                    traceback.print_exc()
                except:
                    pass
                pdb.set_trace()
            if self.resultlist[k]!=0:
                del self.tasklist[k]
                del self.resultlist[k]
                #self.tasklist[k]==0
        if 0 in self.resultlist:
            return 0
        else:
            print "All task finished or failed"
            if self.bugs:
                print "bugs during the run "+str(self.buginfo)
                pdb.set_trace()
            if self.afterprocessing:
                print "processing tasklist result"
                return self.afterprocessing(self.resultlist, self.tasklist,self.other)
            else:
                return self.resultlist
        
    def monitor2end(self):
        try:
            while 0 in self.resultlist:
                try:
                    for k in range(0,len(self.tasklist)):
                        if self.resultlist[k]:
                            continue
                        try:
                            self.resultlist[k]=self.tasklist[k].monitor()
                        except Bugs,e:
                            self.bugs=True
                            self.buginfo=e
                            self.resultlist[k]=987654321
                        #except Exception,e:
                        #    pdb.set_trace()
                        if self.resultlist[k]!=0:
                            #del self.tasklist[k]
                            self.tasklist[k]==0
                    if 0 in self.resultlist:
                        print 'waiting'
                        time.sleep(60)
                    self.delete_runs_on_crashed_nodes()
                except KeyboardInterrupt,e:
                    task.process_keyboard_interrupt(self,e,pos='tasklist-monitorend')
            print "All all task finished or failed"
            if self.bugs:
                print "bugs during the run "+str(self.buginfo)
                raw_input('Input anything to enter debug mode(others finished): ')
                pdb.set_trace()
            if self.afterprocessing:
                print "processing this group of tasks"
                return self.afterprocessing(self.resultlist,self.tasklist,self.other)
            else:
                return self.resultlist
        except FatalError,e:
            self.qdel_jobs(delete=True)
            sys.exit(str(e))
        except KeyboardInterrupt,e:
            task.process_keyboard_interrupt(self,e,pos='tasklist-monitorend')

    def qdel_jobs(self,delete=False):
        for k in range(0,len(self.tasklist)):
            if self.resultlist[k]:
                continue
            self.tasklist[k].qdel_jobs(delete)
               
            
def report_job_runstatus(runpath, runsuccess, runnumber, outputname,inputname='runme.py',temppath=''):
    if runenv.hostn:
        if temppath:
            os.chdir(temppath)
        else:
            temppath=runpath
        if runsuccess:
            fjf=str(runnumber)+outputname
            fl=os.listdir('./')
            nfl=[f for f in fl if f.startswith(fjf)]
            if len(nfl)==0:
                print 'Bugs in code, output file does not exist '
                print fl
                runsuccess=False
                #if temppath!=runpath:
                #    print os.system('rm -rf '+temppath)
            else:
                tr=os.system('tar cvzf '+fjf+'.tar.gz '+fjf+'* --remove-files')
                if tr:
                    print os.system('rm -f '+fjf+'.tar.gz ')
                    tr=os.system('tar cvzf '+fjf+'.tar.gz '+fjf+'*  --remove-files')
                if tr:
                    print os.system('rm -f '+fjf+'*')
                    print 'FatalError: can not tar the result file, disk  full?'
                    runsuccess=False
                else:
                    if temppath!=runpath:
                        cr=os.system('mv '+fjf+'.tar.gz '+runpath)
                        if cr:
                            cr=os.system('mv '+fjf+'.tar.gz '+runpath)
                        if cr:
                            cr=os.system('mv '+fjf+'.tar.gz'+runpath)
                        #print os.system('rm -rf '+temppath)
                        if cr:
                            print os.system('rm -f '+runpath+fjf+'.tar.gz')
                            print 'FatalError: can not copy the result file, netapp disk quota full?'
                            print "FatalError: exiting..."
                            runsuccess=False
        if runsuccess:
            print os.system('touch '+runpath+'jobfinished.'+str(runnumber)+'.'+fjf+'.tar.gz')
        else:
            print os.system('touch '+runpath+'jobfailed.'+str(runnumber)+'.runme.py')
    else:
        if not runsuccess:
            raise Exception('Run failed '+str(runnumber))
            
def generate_job_submit_script(freememory,spname,runtime,nors,parallel=0,mem_total=0,additional_resources=''):
    runmdt=open('runme.sh','w')
    submdt=file('/bell1/home/gqdong/Dropbox/Code/'+scriptname,'r')
    submdtf=submdt.read()
    if parallel:
        submdtf=submdtf.replace('memundefined','#$ -l mem_free='+str(freememory)+'G\n#$ -pe smp '+str(parallel)+'\n'+additional_resources)
    else:
        if mem_total:
            submdtf=submdtf.replace('memundefined','#$ -l mem_free='+str(freememory)+'G\n'+'#$ -l mem_total='+str(mem_total)+'G\n'+additional_resources)
        else:
            submdtf=submdtf.replace('memundefined','#$ -l mem_free='+str(freememory)+'G\n'+additional_resources)
    #submdtf=submdtf.replace('memundefined','#$ -l mem_free='+str(freememory)+'G\n')
    if "#" in spname:
        nspname=spname.replace('#','')
        submdtf=submdtf.replace('undefinedrundir/output','output/'+nspname)
        print os.system('ssh '+runenv.jobserver+' mkdir '+'/netapp/sali/gqdong/output/'+nspname)
    submdtf=submdtf.replace('undefinedrundir',spname)
    submdtf=submdtf.replace('timeundefined',runtime)
    submdtf=submdtf.replace('#$ -p 0',"#$ -p "+str(runenv.runpriority))
    submdtf=submdtf.replace('tundefined',str(1)+'-'+str(nors))
    runmdt.write(submdtf)
    runmdt.flush()