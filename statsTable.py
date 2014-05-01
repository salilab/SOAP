"""
   SOAP statistical table

"""

from feature import *
from env import *
from sequences import *
                    
class rawsp(object):
    def __init__(self,pdbset='',features='',genmethod='',tlib=None,decoy=False,routine='calc_sp_all', model={}):
        if model:
            features=model['features']
            genmethod=model['genmethod']
            pdbset=model['pdbset']
        if len(pdbset)>0 and pdbset[0]=='o': # o indicate that we want to use unprocessed pdb.
            self.opdb=True
            pdbset=pdbset[1:]
        else:
            self.opdb=False
        self.features=features
        self.genmethod=genmethod
        self.pdbset=pdbset
        self.decoy=decoy
        self.routine=routine
        self.tlib=tlib
        if self.opdb:
            pdbset='o'+pdbset
        if runenv.hostn==0:
            try:
                self.dir=runenv.basedir+pdbset+'/'+features+'-'+genmethod+'/'
            except:
                pdb.set_trace()
            self.spname=pdbset+'.'+features+'.'+genmethod
            self.sppath=self.dir+self.spname
        elif runenv.hostn==1:
            self.dir=runenv.basedir
            self.spname=pdbset+'.'+features+'.'+genmethod
            self.sppath=self.dir+self.spname
        self.nors=200
        self.flagdict={}
        self.currentdsname='dummy'
        self.genmethodlist=decode_genmethod(self.genmethod)
        if self.features:
            self.fo=feature(self.features)
            self.bins=self.fo.get_all_bins_centers()
            for i in range(0,len(self.bins)):
                if len(self.bins[i])>0:
                    self.bl=len(self.bins[i])
                    self.bp=i
            rl1,rl2=self.fo.get_featuretypepos()
            self.ifl=len(rl1)
            self.rs1=rl1+rl2
            self.rs2=copy.copy(self.rs1)
            for i in range(0,len(self.rs1)):
                self.rs2[self.rs1[i]]=i
            if self.rs1==range(0,len(self.rs1)):
                self.permute=False
            else:
                self.permute=True
            glib=self.fo.get_lib()
            self.tlib=glib
        self.bndsepcalc=False
        self.runpath='./'
        if self.pdbset=='Otherpot':
            self.mkdir()
        self.justwait=False
        self.waittarget=[]
        if routine=='calc_sp_i':
            self.prd='no-backupi'
            self.ftype='.npy'
        else:
            self.prd='no-backup'
            self.ftype='.hdf5'
        self.mkdir() #make data dir: IT IS NOT THE RUNDIR-THE TEMP DIR
        self.isatompairclustering()
          
    def mkdir(self):
        if runenv.hostn:
            return 0
        basedir=runenv.basedir
        pdbset=self.pdbset
        if self.opdb:
            pdbset='o'+pdbset
        genmethod=self.genmethod
        features=self.features
        if not (os.access(basedir+pdbset,os.F_OK)):
            os.mkdir(basedir+pdbset)
        if not (os.access(basedir+pdbset+'/'+features+'-'+genmethod,os.F_OK)):
            os.mkdir(basedir+pdbset+'/'+features+'-'+genmethod)
        os.chdir(basedir+pdbset+'/'+features+'-'+genmethod)
        if self.tlib:
            print os.system('cp '+runenv.libdir+self.tlib+' ./')

    def isatompairclustering(self):
        genmethod=self.genmethod
        self.atompairclustering=False;
        self.apca=0
        if re.search('(ap[0-9]{1,2})',genmethod):
            rer=re.search('(ap[0-9]{1,2})',genmethod)
            apname=rer.group(1)
            self.atompairclustering=True
            self.apca=np.load(runenv.libdir+apname+'.npy')

    def subset_of_existing(self):
        if not os.path.isdir(runenv.basedir+self.pdbset):
            return []
        rangelist=self.get_allsp_range()
        nal=[]
        rl=[]
        revelentsp=[]
        for sprange in rangelist:
            if sprange[2][0:4]==self.genmethodlist[0:4] and sprange[2][-2:]==self.genmethodlist[-2:]:
                subsetindex=self.fo.issubset(feature(sprange[1]))
                if len(subsetindex)==0:
                    continue
                bsrange=sprange[2][4:6]
                if not bsrange in rl:
                    rl.append(bsrange)
                    revelentsp.append([sprange[0],bsrange,subsetindex])
                    na=np.zeros([20,1])
                    if bsrange[0]==-1 and bsrange[1]==-1:
                        na[...]=1
                    else:
                        na[bsrange[0]:min(bsrange[1]+1,20)]=1
                    nal.append(na)
        if len(nal)==0:
            return []
        existingsparray=np.hstack(nal)
        na=np.zeros([20,1])
        if self.genmethodlist[4]==-1 and self.genmethodlist[5]==-1:
            na[...]=1
        else:   
            na[self.genmethodlist[4]:min(self.genmethodlist[5]+1,20)]=1
        mixweight=numpy.linalg.lstsq(existingsparray,na)[0]
        residue=sum((np.dot(existingsparray,mixweight)-na)**2)
        logstr=''
        if np.allclose(residue,0):
            spl=[]
            for r1,r2 in zip(list(mixweight),revelentsp):
                if np.abs(r1[0])>0.0001:
                    spl.append([r1[0],r2[0],r2[2]]) #target
                    logstr=str(r1[0])+'*'+str(r2[1])+' + '
            print logstr
            return spl
        else:
            return []
                    
    def get_allsp_range(self):
        fl=os.listdir(runenv.basedir+self.pdbset)
        fl=[f for f in fl if os.path.isdir(runenv.basedir+self.pdbset+'/'+f)]
        rangelist=[]
        for f in fl:
            fs=f.split('-')
            fs=['-'.join(fs[:-1]),fs[-1]]
            spdir=runenv.basedir+self.pdbset+'/'+f+'/'+'.'.join([self.pdbset]+fs)
            if os.path.isfile(spdir+self.ftype) or os.path.isdir(runenv.basedir+self.pdbset+'/'+f+'/'+self.prd):
                rangelist.append([spdir, fs[0],decode_genmethod(fs[1])])
        return rangelist
        
    def calc_sp_all(self,runnumber):
        """Methods to generate the statistical potential from PDB files  in the currnt directory."""
        self.runpath=get_path_to_thisfile()
        features=self.features
        pdbfile='pdbs'+str(runnumber)
        print 'running pdb: '+pdbfile
        env = environ()
        #log.minimal()
        io=env.io
        if self.opdb:
            pdbfdir=runenv.opdbdir
            io.atom_files_directory=[scratchdir]
        else:
            pdbfdir=runenv.pdbdir
            io.atom_files_directory=pdbfdir+[scratchdir]
        mlib=mdt.Library(env)
        mlib.bond_classes.read('${LIB}/bndgrp.lib')
        featurelist=feature(self.features,mlib).get_featurelist()
        #print featurelist
        m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
        fl=len(featurelist)
        a=alignment(env)
        f=modfile.File(self.runpath+pdbfile,'r')
        csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(self.genmethod)
        refy=0
        rmdt=None
        if re.search('ref(.*)$',self.genmethod):
            env2=environ()
            mlib2=mdt.Library(env2)
            featurelist=feature(self.features,mlib2).get_featurelist()
            refrer=re.search('ref(.*)$',self.genmethod)
            refname=refrer.group(1)
            refy=1
            rmdt=mdt.Table(mlib2,features=featurelist,shape=[-1]*len(featurelist))
            rmdt.read_hdf5(refname+'.h5.lib')
        k=0
        lstartpos=0
        lendpos=0
        lscale=0
        if re.search('s([0-9\.]*)e([0-9\.]*)s([0-9\.]*)',self.genmethod):
            rer=re.search('s([0-9\.]*)e([0-9\.]*)s([0-9\.]*)',self.genmethod)
            lstartpos=float(rer.group(1))
            lendpos=float(rer.group(2))
            lscale=float(rer.group(3))
        while a.read_one(f, allow_alternates=True):
            k=k+1
            print k
            fd,naln=self.get_pdbfile(pdbfdir,[],a[0].code,env,a,scratchdir)
            try:
                if self.opdb:
                    aln=alignment(env)
                    mdl=model(env)
                    if len(a[0].code)==5:
                        cc=a[0].code[-1]
                        m2=model(env)
                        m2.read(a[0].code[0:4])
                        cl=[c.name for c in m2.chains]
                        ci=ord(cc)-65
                        cc=cl[ci]
                        mdl.read(a[0].code[0:4],model_segment=('FIRST:'+cc,'LAST:'+cc))
                    elif len(a[0].code)==4:
                        mdl.read(a[0].code[0:4])
                    else:
                        raise Bugs('Do not know how to calculate mdt table using opdb for '+a[0].code)
                    aln.append_model(mdl,align_codes=a[0].code)
                    toaln=aln
                else:
                    toaln=a
                print toaln[0].code
                m.add_alignment_witherr(toaln,chain_span_range=(-csmax,-csmin,csmin,csmax),residue_span_range=(-kcut1,-kcut,kcut,kcut1),bond_span_range=(bs,bsm),errorscale=errsv,disulfide=ssp)
            except Exception,e:
                traceback.print_exc()
                if re.search('fork failure',str(e)):
                    raise Exception('fork failure, quiting')
        pdbfilel=re.match('.*(pdb.*)',pdbfile)
        pdbfilename=pdbfilel.group(1)
        if pdbfilename[-4:]=='.pir':
            pdbfilename=pdbfilename[:-4]
        if m.sum() != 0:
            try:
                temppath=scratchdir
                runsuccess=True
                m.write_hdf5(temppath+str(runnumber)+'.result.hdf5')
            except:
                try:
                    temppath='/scrapp2/gqdong/'+os.getenv('JOB_ID')+'/'
                    os.makedirs(temppath)
                    print os.system('mkdir '+temppath)
                    m.write_hdf5(temppath+str(runnumber)+'.result.hdf5')
                    runsuccess=True
                except:
                    runsuccess=False
        else:
            runsuccess=False
            print m.sum()
            print 'mdt sum to 0, code problem'
        report_job_runstatus(self.runpath, runsuccess, runnumber, '.result',inputname='runme.py',temppath=temppath)
            
    def calc_sp_cluster(self):
        sj=self.get_task()
        tasklist(tasklist=[sj]).monitor2end()

    def get_task(self):
        if re.search('(ap[0-9]+n[0-9]+)$',self.genmethod):
            rer=re.search('(ap[0-9]+)(n[0-9]+)$',self.genmethod)
            ngm=self.genmethod[:-len(rer.group(2))]
            nsp=rawsp(pdbset=self.pdbset,features=self.features,genmethod=ngm,tlib=self.tlib,decoy=self.decoy,routine=self.routine)
            return nsp.get_task()
        print 'getting task '+self.spname
        if self.routine=='calc_sp_all':
            self.targetfile=self.sppath+'.hdf5'
            self.filetype='.hdf5'
            tc=os.path.isfile(self.sppath+'.hdf5')
            self.nors=200        
        elif self.routine=='calc_sp_i':
            pdbsets=self.pdbset.split('_')
            if len(pdbsets)==1:
                self.targetfile=self.sppath+'.npy'
                self.filetype='.npy'
                tc=os.path.isfile(self.sppath+'.npy')
                self.nors=1000
            else:
                tl=[]
                for pset in pdbsets:
                    nspo=rawsp(pset,self.features,self.genmethod,self.tlib,self.decoy,self.routine)
                    tl.append(nspo.get_task())
                print "running break down sps"
                return tasklist(tasklist=tl,afterprocessing=dummy)
            self.nors=1000
        if not tc:
            sj=task(self.dir+self.prd+'/',self.spname,afterprocessing=self,preparing=self,targetfile=self.sppath+self.filetype)
            return sj
        else:
            return 0
        
    def afterprocessing(self):
        print "run finished -> processing result"
        if self.justwait and len(self.waittarget)==0:
            pass
        elif self.justwait and len(self.waittarget)>0: # the sp need to be combined with others
            os.chdir(self.dir+self.prd+'/')
            self.get_targetsp_fromexistingsps()
        else:
            os.chdir(self.dir+self.prd+'/')
            if self.routine=='calc_sp_all':
                self.sp_sum()
                print os.system('cp sum.hdf5 '+self.sppath+'.hdf5')
                if 0:
                    self.sp_stderr2(nolog=True)
                    print os.system('mv ostderr2.hdf5 '+self.sppath+'.ostderr2.hdf5')
                    print os.system('mv omean.hdf5 '+self.sppath+'.omean.hdf5')                   
                    self.sp_stderr2(nolog=False)
                    print os.system('mv stderr2.hdf5 '+self.sppath+'.stderr2.hdf5')
                    print os.system('mv mean.hdf5 '+self.sppath+'.mean.hdf5')     
            elif self.routine=='calc_sp_i':
                allarray=self.cat_spi()
                if self.atompairclustering:
                    #process atom pair features
                    for i in range(allarray.shape[1]):
                        sdir=runenv.basedir+self.pdbset+'/'+self.features+'-'+self.genmethod+'n'+str(i)+'/'
                        scorepath=sdir+self.pdbset+'.'+self.features+'.'+self.genmethod+'n'+str(i)
                        if not os.path.isdir(sdir):
                            os.mkdir(sdir)
                        np.save(scorepath,allarray[:,i,...].squeeze())
                print os.system('mv all.npy ../'+self.spname+'.npy')#
            print os.system('rm -r *')
            os.chdir('..')
        print os.system('rm -r '+self.dir+self.prd)
        return 1

    def get_targetsp_fromexistingsps(self):
        sparray=[]
        for target in self.waittarget:
            sa=self.get_sp_subset(target[1],target[2])
            try:
                if len(sparray)==0:
                    if self.routine=='calc_sp_i':
                        nas=[item for item in [sa.shape[0]]+self.fo.fdlist if item!=1]
                        sparray=np.zeros(nas)
                    else:
                        sparray=np.zeros(self.fo.fdlist)
                sparray=sparray+target[0]*sa
            except:
                print traceback.print_exc()
                pdb.set_trace()
        self.save_npy(sparray)

    def get_empty_table(self,write=True):
        env = environ()
        mlib=mdt.Library(env)
        featurelist=feature(self.features,mlib).get_featurelist()
        m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
        if write:
            m.write_hdf5(self.sppath+'.hdf5')
        return m

    def save_npy(self,sparray):
        if self.routine=='calc_sp_all':
            self.get_empty_table()
            f=h5py.File(self.sppath+'.hdf5','r+')
            mdtb=f['mdt']
            mdtb[...]=sparray[...]
            f.close()            
        elif self.routine=='calc_sp_i':
            np.save(self.sppath+'.npy',sparray)
            
    def get_sp_subset(self,f,subsetindex):
        if self.routine=='calc_sp_i':
            subsetindex[0]=[[0]]+subsetindex[0]
        mdta=self.load2npy(f)
        return self.get_array_subset(mdta,subsetindex)

    def get_array_subset(self,mdta,subsetlist):
        subsetindex=subsetlist[0]
        subsetindex.reverse()
        nsubsetindex=copy.deepcopy(subsetindex)
        rr=range(len(subsetindex))
        rr.reverse()
        for si,i in zip(subsetindex,rr): #eliminate the dimensions summed over first
            if len(si)==0:
                mdta=mdta.sum(axis=i)
                nsubsetindex.pop(i)
        rr=range(len(nsubsetindex))
        rr.reverse()
        if len(mdta.shape)!=len(nsubsetindex):
            for si,i in zip(subsetindex,rr): #delete dimension with size 1
                if si==[0,1]:
                    nsubsetindex.pop(i)
            rr=range(len(nsubsetindex))
            rr.reverse()
        try:
            for si,i in zip(nsubsetindex,rr):
                nas=list(mdta.shape)
                if len(si)>1 and nas[i]>(len(si)-1):
                    nas[i]=len(si)-1
                    na=np.zeros(nas)
                    for j,ej in zip(range(si[0],si[1]),range(si[-2]+1,si[-1]+1)):
                        nai=[':' for item in rr]
                        nai[i]=str(j)+':'+str(ej)+':'+str(si[1]-si[0])
                        try:
                            na=na+eval('mdta['+','.join(nai)+']')#there must be a simple way that does not need permuation.
                        except:
                            print "line 1071 in sp"
                            pdb.set_trace()
                    mdta=na
        except:
            print traceback.print_exc()
            pdb.set_trace()
        if subsetlist[1]!=range(len(subsetlist[1])):
            na=np.transpose(na,subsetlist[1])
        return na
                    
    def load2npy(self,f):
        print "loading "+f
        if self.routine=='calc_sp_all':
            f=h5py.File(f,'r')
            mdtb=f['mdt']
            ms=mdtb.shape
            print ms
            mdta=np.zeros(np.array(list(ms)))
            mdta[...]=mdtb[...]
            f.close()
        elif self.routine=='calc_sp_i':
            mdta=np.load(f)
        return mdta
            
    def calc_sp_local(self,nors,routine='calc_sp_all'):
        self.routine=routine
        self.nors=nors
        self.filetype='.npy'
        self.targetfile=self.spname+'.npy'
        self.prepare_task_input()
        sj=task(self.dir+'no-backup/')
        if routine=='calc_sp_all':
            sj.run_task_local(self.calc_sp_all)
            self.sp_mean()
            #self.sp_stderr()
            print os.system('mv sum.hdf5 ../'+self.spname+'.hdf5')
            #print os.system('mv stderr.hdf5 ../'+self.spname+'.stderr.hdf5')
        elif routine=='calc_sp_i':
            sj.run_task_local(self.calc_sp_i)
            self.cat_spi()
            print os.system('mv all.npy ../'+self.spname+'.npy')
        
    def sp_sum(self,wdir=None):
        if wdir:
            #wdir=self.sppath
            os.chdir(wdir)
        #combine mdts in the current directory
        env = environ()
        mlib=mdt.Library(env)
        featurelist=feature(self.features,mlib).get_featurelist()
        m1=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
        m2=mdt.Table(mlib, features=featurelist,shape=[-1]*len(featurelist))
        k=0
        files=os.listdir('./')
        for name in files:
            try:
                if name.endswith('.result.hdf5'):
                    fname=name
                else:
                    continue
                print 'Adding '+name
                k=k+1          
                if k==1:
                    m1.read_hdf5(fname)
                else:
                    m2.read_hdf5(fname)
                    m1+=m2
                print 'Adding '+name+'Finished'
            except Exception, e:
                pdb.set_trace()
                traceback.print_exc()
                return 1
        m1.write_hdf5('sum.hdf5')
        m1=[]
        m2=[]    
        print "Combining finished"
        return 0

    def should_wait(self):
        if self.atompairclustering:
            return []
        self.waittarget=self.subset_of_existing()
        if len(self.waittarget)>0: # the sp can be calculated by combining other sps
            self.justwait=True
            wl=[item[1]+self.filetype for item in self.waittarget]
            for item in self.waittarget:
                item[1]=item[1]+self.filetype
            print "waiting for other runs to finish to combine to our sp"
            print wl
            self.mkdir()
            print os.system('rm -r '+self.dir+self.prd)
            os.mkdir(self.dir+self.prd)
            return wl
        else:
            return []
    
    def prepare_task_input(self):
        wl=self.should_wait()
        if len(wl)>0:
            return wl
        nors=self.nors
        routine=self.routine
        print "prepare task input for generating"+self.spname
        os.chdir(self.dir)
        #os.system('rm -r '+self.prd)
        freememory=feature(self.features).get_runmem()
        self.runpath='/netapp/sali/gqdong/'+self.spname+'/'
        if freememory>100: #put the temporary time somewhere else when the resulting file is large.
            temp='/bell4/gqdong/'+self.spname
            os.mkdir(temp)
            os.system('ln -s '+temp+' '+'./'+self.prd)
        else:
            if not (os.access(self.prd,os.F_OK)):
                os.mkdir(self.prd)
        freememory=freememory+0.2
        if freememory>2:
            self.nors=min(int(1800/freememory),50)
        else:
            self.nors=nors
        nors=self.nors
        os.chdir(self.prd)
        if self.tlib:
            os.system('cp '+runenv.libdir+self.tlib+' ./')
        if 'l' in self.fo.fclist:
            residuepower=1
        else:
            residuepower=2
        if self.decoy:
            pir(decoysets( self.pdbset.split('_')).pirpath).sep_pir('./',nors,permin=scoreminnumber,residuepower=residuepower)
        else:
            print os.system('cp '+runenv.basedir+'pdbpir/pdb_'+self.pdbset+'.pir ./')
            pir('pdb_'+self.pdbset+'.pir').sep_pir('./',nors,residuepower=residuepower)
        nors=int(os.popen('ls | grep -c pdbs').read())
        self.nors=nors
        os.system('cp '+runenv.scriptpath+' ./')
        inputlist=open('inputlist','w')
        inputlist.write(','.join(['pdbs'+str(i) for i in range(1,nors+1)]))
        inputlist.close()
        makemdt=open('runme.py','w')
        makemdt.write('keeptry=5\nwhile keeptry:\n    try:\n        from SOAP.statsTable import *\n        keeptry=0\n    except:\n        keeptry-=1\n'
                      +'import sys \n \nrsp=rawsp(\''+self.pdbset+'\',\''
                      +self.features+'\',\''+self.genmethod+'\')\n'+'rsp.opdb='+str(self.opdb)+'\n'
                      +'rsp.'+routine+'(sys.argv[1])')
        makemdt.flush()
        generate_job_submit_script(freememory,self.spname,runtime,nors)
        return []
    
    def calc_sp_i(self,runnumber,reportstatus=True):
        self.runpath=get_path_to_thisfile()
        pirfile='pdbs'+str(runnumber)
        dnlist=pir(pirfile).get_codelist()
        if self.atompairclustering:
            idist=self.calc_individual_sp(pirfile,dnlist,'apidist')
        else:
            idist=self.calc_sp_individual1D_frompir(pirfile,dnlist)           
        if reportstatus:
            try:
                temppath=scratchdir
                runsuccess=True
                np.save(temppath+str(runnumber)+'.result.npy',idist)
            except:
                try:
                    temppath='/scrapp/gqdong/'+os.getenv('JOB_ID')+'/'
                    print os.system('mkdir '+temppath)
                    np.save(temppath+str(runnumber)+'.result.npy',idist)
                    runsuccess=True
                except:
                    runsuccess=False
            report_job_runstatus(self.runpath, True, runnumber, '.result',inputname='runme.py',temppath=temppath)
        else:
            return idist
        
    def calc_sp_individual1D_frompir(self, pirfile,dnlist):
     #generate distributions for calculate benchmark score for reference state
        features=self.features
        bm=self.genmethod
        tmpdir=scratchdir
        env = environ()
        log.minimal()
        mlib=mdt.Library(env)
        mlib.bond_classes.read('${LIB}/bndgrp.lib')
        env.io.atom_files_directory=[scratchdir,'./']+runenv.pdbdir
        io=env.io
        featurelist=feature(self.features,mlib).get_featurelist()
        csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(self.genmethod)
        m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
        fl=len(featurelist)
        #m=m.reshape(featurelist,[0 for i in range(0,fl)], [-1 for i in range(0,fl)])
        tal=1
        ms=m.shape
        for si in ms:
            tal=tal*si
        refdista=np.zeros([len(dnlist)]+[tal])
        f=modfile.File(pirfile,'r')
        aln=alignment(env)
        k=0
        while aln.read_one(f):
            sr=0;
            alncode=aln[0].code
            print alncode
            codelist=alncode.split('.')
            if len(codelist)>1:#decoys
                print runenv.ddbdir+codelist[0]+'/'+codelist[1]+'/'
                fd,naln=self.get_pdbfile(runenv.ddbdir+codelist[0]+'/'+codelist[1]+'/',codelist,alncode,env,aln,tmpdir)
                io.atom_files_directory=[fd]+io.atom_files_directory
            else:#pdbs
                io.atom_files_directory=runenv.pdbdir+io.atom_files_directory
                #fd,naln=self.get_pdbfile(runenv.pdbdir,[],alncode,env,aln,tmpdir)
            print k
            m.add_alignment_witherr(naln,chain_span_range=(-csmax,-csmin,csmin,csmax),residue_span_range=(-kcut1,-kcut,kcut,kcut1),bond_span_range=(bs,bsm),errorscale=errsv,io=io,disulfide=ssp)#,startpos=0,endpos=0,scale=0,refyes=False,refmdt='',io=io)
            if len(m.shape)>1 and (np.array(m.shape)>1).sum()>1:
                raise Bugs('Currently the code does not support 2-d rawsp-i')
            else:
                for i in range(0,tal):
                    if len(m.shape)==2:
                        refdista[k,i]=m[0,i]
                    #np.unravel_index(i,ms)
                    else:
                        refdista[k,i]=m[i]
            k=k+1
            m=[]
            aln=[]
            m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
            fl=len(featurelist)
            #m=m.reshape(featurelist,[0 for i in range(0,fl)], [-1 for i in range(0,fl)])
            aln=alignment(env)
        #print os.system('rm -rf '+tmpdir)
        return refdista
                
    def calc_sp_individual1D_fromlist(self, dnlist):
     #generate distributions for calculate benchmark score for reference state
        features=self.features
        bm=self.genmethod
        env = environ()
        log.minimal()
        mlib=mdt.Library(env)
        mlib.bond_classes.read('${LIB}/bndgrp.lib')
        env.io.atom_files_directory=[scratchdir,'./'] +runenv.pdbdir
        featurelist=feature(self.features,mlib).get_featurelist()
        csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(self.genmethod)
        k=0
        m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
        fl=len(featurelist)
        #m=m.reshape(featurelist,[0 for i in range(0,fl)], [-1 for i in range(0,fl)])
        refdista=np.zeros([len(dnlist)]+ m.shape)
        tal=1
        ms=m.shape
        for si in ms:
            tal=tal*si
        for k in range(0,len(dnlist)):
            mdl=model(env)
            mdl.read(file=dnlist[k])
            aln=alignment(env)
            codelist=dnlist[k].split('.')
            if len(codelist)>1:
                print runenv.ddbdir+codelist[0]+'/'+codelist[1]+'/'
                fd,naln=self.get_pdbfile(runenv.ddbdir+codelist[0]+'/'+codelist[1]+'/',codelist,dnlist[k],env,aln)
                env.io.atom_files_directory=[fd]+env.io.atom_files_directory
            else:
                fd,naln=self.get_pdbfile(runenv.ddbdir+codelist[0]+'/'+codelist[1]+'/',[],alncode,env,aln,tmpdir)
            aln.append_model(mdl,align_codes=dnlist[k],atom_files=dnlist[k])
            print k
            #if k> 100:
            #    return refdista
            m.add_alignment_witherr(aln,chain_span_range=(-csmax,-csmin,csmin,csmax),residue_span_range=(-kcut1,-kcut,kcut,kcut1),bond_span_range=(bs,bsm),errorscale=errsv,disulfide=ssp)
            for i in range(0,tal):
                id=np.unravel_index(i,ms)
                rid=tuple([k]+list(id))
                refdista[rid]=m[id]
            m=[]
            aln=[]
            mdl=[]
            m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
            fl=len(featurelist)
        return refdista
    
    def calc_individual_sp(self, pirfile,dnlist,ptype):
        #get the mdt tables only for ddb
        features=self.features
        env = environ()
        tmpdir=scratchdir     
        log.minimal()
        mlib=mdt.Library(env)
        mlib.bond_classes.read('${LIB}/bndgrp.lib')
        env.io.atom_files_directory=[scratchdir,'./']+runenv.pdbdir
        io=env.io
        featurelist=feature(self.features,mlib).get_featurelist()
        m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
        ma=m.get_array_view()
        mshape=[]
        mshape=[i for i in m.shape if i>1]
        m1s=m.shape[0]
        m2s=1
        for item in m.shape[1:]:
            m2s*=item
        #decode generation method
        bm=self.genmethod
        #initialize the array for storing result
        if ptype=='idist':
            ra=np.zeros([len(dnlist)]+ mshape)
            csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(self.genmethod)
        elif ptype=='apidist':
            apc=self.apca
            apcl=[]
            apcn=apc.max()+1
            pairdim=len(apc.shape)
            sds=list(mshape[:-pairdim])
            sdl=[]
            for i in range(len(sds)):
                indexshape=[1 for j in range(len(sds)+1)]
                indexshape[i]=sds[i]
                sdl.append(np.arange(0,sds[i]).reshape(indexshape))
            for i in range(apcn):
                cdl=sdl+[]
                apci=np.nonzero(apc==i)
                indexshape=[1 for j in range(len(sds)+1)]
                indexshape[-1]=len(apci[0])
                for item in apci:
                    cdl.append(item.reshape(indexshape))
                apcl.append(tuple(cdl))
            ra=np.zeros([len(dnlist),apcn]+ sds)
            csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(self.genmethod)
        elif ptype=='apdecompscore':
            apc=self.apca
            apcl=[]
            apcn=apc.max()+1
            pairdim=len(apc.shape)
            sds=list(mshape[:-pairdim])
            sumaxis=range(1,len(mshape)-pairdim)
            sumaxis.reverse()
            for i in range(apcn):
                cdl=[]
                apci=np.nonzero(apc==i)
                for item in apci:
                    cdl.append(item.reshape([len(item)]))
                apcl.append(tuple(cdl))
            ms=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
            ms.read_hdf5(self.ssppath+'.hdf5')
            msa=ms.get_array_view()
            self.svdu,self.svdv=self.get_svd()
            ra=np.zeros([len(dnlist),apcn,self.svdu.shape[0]])
            csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(self.bm)            
        elif ptype=='apscore':
            apc=self.apca
            apcl=[]
            apcn=apc.max()+1
            pairdim=len(apc.shape)
            sds=list(mshape[:-pairdim])
            sumaxis=range(0,len(mshape)-pairdim)
            sumaxis.reverse()
            for i in range(apcn):
                cdl=[]
                apci=np.nonzero(apc==i)
                for item in apci:
                    cdl.append(item.reshape([len(item)]))
                apcl.append(tuple(cdl))
            ms=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
            ms.read_hdf5(self.ssppath+'.hdf5')
            msa=ms.get_array_view()
            ra=np.zeros([len(dnlist),apcn])
            csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(self.bm)
        elif ptype=='decompscore':
            self.svdu,self.svdv=self.get_svd()
            ra=np.zeros([len(dnlist),self.svdu.shape[0]])
            csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(self.bm)
        #loop through all the structures to generate the table for individual structures
        f=modfile.File(pirfile,'r')
        aln=alignment(env)
        k=0
        while aln.read_one(f):
            sr=0;
            alncode=aln[0].code
            print alncode
            codelist=alncode.split('.')
            if len(codelist)>1:
                print runenv.ddbdir+codelist[0]+'/'+codelist[1]+'/'
                fd,naln=self.get_pdbfile(runenv.ddbdir+codelist[0]+'/'+codelist[1]+'/',codelist,alncode,env,aln,tmpdir)
                io.atom_files_directory=[fd]+io.atom_files_directory
            else:
                naln=aln
                io.atom_files_directory=runenv.pdbdir+io.atom_files_directory
            print k
            m.add_alignment_witherr(naln,chain_span_range=(-csmax,-csmin,csmin,csmax),residue_span_range=(-kcut1,-kcut,kcut,kcut1),bond_span_range=(bs,bsm),errorscale=errsv,io=io,disulfide=ssp)#,startpos=0,endpos=0,scale=0,refyes=False,refmdt='',io=io)
            #process the result from this single table
            #sma=np.squeeze(np.copy(ma))
            sma=ma
            if ptype=='idist':
                try:
                    ra[k,...]=sma[i]
                except:
                    pdb.set_trace()
            elif ptype=='apdecompscore':
                rso=(np.dot(self.svdu.T,sma.reshape([m1s,m2s]))*self.svdv).reshape(m.shape)
                for ax in sumaxis:
                    rso=rso.sum(axis=ax)
                ra[k,:]=rso[:,apcl[i]].sum(axis=1)
            elif ptype=='apidist':
                for i in range(apcn):
                    ra[k,i,...]=sma[apcl[i]].sum(axis=-1)
            elif ptype=='apscore':
                for i in range(apcn):
                    nma=(sma*msa)
                    for ax in sumaxis:
                        nma=nma.sum(axis=ax)
                    ra[k,i]=nma[apcl[i]].sum()
            elif ptype=='decompscore':
                ra[k,:]=(np.dot(self.svdu.T,sma.reshape([m1s,m2s]))*self.svdv).sum(axis=1)
            k=k+1
            m.clear()
        #print os.system('rm -rf '+tmpdir)
        return ra
                    
    def cat_spi(self,name='result'):
        filelist=[str(item)+'.'+name+'.npy' for item in range(1, self.nors+1)]
        al=[]
        alen=0
        aw=0
        alp=[]
        for name in filelist:
            temp=np.load(name)
            alen=alen+temp.shape[0]
            aw=temp.shape[1:]
            alp.append(temp.shape[0])
            al.append(name)
        ddist=np.zeros([alen]+list(aw),dtype=np.float32)
        sp=0
        for i in range(0,len(al)):
            ddist[sp:sp+alp[i]]=np.load(al[i])
            sp=sp+alp[i]
        np.save('all.npy',ddist)
        return ddist
    
    def sp_stderr2(self,wdir=None,nolog=False):
        ##!!!!!!!!!!!!!!!will calculate the stderr here 
        if wdir:
            os.chdir(wdir)
        print "Calculating Standard deviation"
        filelist=[]
        files=os.listdir('./')
        for filename in files:
            if filename!='sum.hdf5' and filename.endswith('.result.hdf5'):
                filelist.append(filename)
        if nolog:
            prefix='o'
        else:
            prefix=''
        os.system('cp sum.hdf5 '+prefix+'mean.hdf5')
        os.system('cp sum.hdf5 '+prefix+'stderr2.hdf5')
        f=[]
        f2=[]
        f=h5py.File(prefix+'mean.hdf5','r+')
        f2=h5py.File(prefix+'stderr2.hdf5','r+')
        mdtmean=f['mdt']
        mdtmean[...]=0
        if nolog:
            ssp=scaledsp('nbsumnolog',self)
        else:
            ssp=scaledsp('nbsum',self)
        for filename in filelist:
            print 'calcting error for '+filename
            ssp.scalesp(ssppath=filename[:-5])
        for filename in filelist:
            #try:
            ft=h5py.File(filename,'r')
            mdtmean[...]=mdtmean[...]+ft['mdt'][...]
            ft.close()
            ft=[]
            #except:
            #    print 'error in sp_stderr2 1'
            #    pdb.set_trace()
        mdtmean[...]=mdtmean[...]/len(filelist)
        mdtstd=f2['mdt']
        mdtstd[...]=0
        for filename in filelist:
            #try:
            ft=h5py.File(filename,'r')
            mdtstd[...]=mdtstd[...]+(ft['mdt'][...]-mdtmean[...])**2
            ft.close()
            ft=[]
            #except:
            #    print 'error in sp_stderr2 2'
            #    pdb.set_trace()
        mdtstd[...]=mdtstd[...]/len(filelist)**2
        f.flush()
        f.close()
        f2.flush()
        f2.close()
        print "Calculating Standard error square Finished"

    def load_i(self):
        idist=np.load(self.sppath+'.npy')
        ms=idist.shape
        if 0 and  self.bl<ms[self.bp+1]: #outdated
            print "reshaping the original stats array"
            if len(ms)==2:
                idist=np.copy(idist[:,:-1])
            elif len(ms)==3:
                idist=np.copy(idist[:,:-1,:-1])
            elif len(ms)==4:
                idist=np.copy(idist[:,:-1,:-1,:-1])
            elif len(ms)==5:
                idist=np.copy(idist[:,:-1,:-1,:-1,:-1])
            np.save(self.sppath+'.npy',idist)
        if idist.dtype!='float32':
            idist=idist.astype(np.float32)
            np.save(self.sppath+'.npy',idist)
        return idist

    def get_i(self,nors=600):
        self.nors=nors
        self.routine='calc_sp_i'
        pdbsets=self.pdbset.split('_')
        nalen=0
        if len(pdbsets)==1:
            if not os.path.isfile(self.sppath+'.npy'):
                self.calc_sp_cluster()
            return self.load_i()
        else:
            if os.path.isfile(self.sppath+'.npy'):
                return self.load_i()
            else:
                al=[]
                self.calc_sp_cluster()
                for pset in pdbsets:
                    nspa=rawsp(pset,self.features,self.genmethod,self.tlib,self.decoy,self.routine).load_i()
                    al.append(nspa)
                    nalen+=nspa.shape[0]
                if len(nspa.shape)>1:
                    newshape=tuple([nalen]+list(nspa.shape[1:]))
                else:
                    newshape=(nalen)
                print newshape
                na=np.empty(newshape,dtype=np.float32)
                k=0
                for nspa in al:
                    na[k:k+nspa.shape[0]]=nspa[...]
                    k+=nspa.shape[0]
                #os.makedirs(os.path.split(self.sppath)[0])
                np.save(self.sppath+'.npy',na)
                return na
    
    def combine_flag(self,pdbdir):
        pass
    
    def get_pdbfile(self,pdbdir,codelist,alncode,env,aln,tmpdir='./'):
        code=alncode
        if isinstance(pdbdir, list) and len(pdbdir)==1:
            pdbdir=pdbdir[0]
        print "unziping "+code
        if isinstance(pdbdir, list):
            return [pdbdir,aln]
        if os.path.isfile(pdbdir+'needattachtobase'):
            basefile=pdbdir+codelist[1]+'.base.pdb'
            if not os.path.isfile(basefile):
                basefile=pdbdir+codelist[0]+'.'+codelist[1]+'.base.pdb'
            fh=open(basefile,'r')
            fc=fh.read()
            fh.close()
            fh=gzip.open(pdbdir+alncode+'.pdb.gz','r')
            dfc=fh.read()
            fh.close()
            fc=fc+'\n'+dfc
            fh=open(os.path.join(tmpdir,alncode+'.pdb'),'w')
            fh.write(fc)
            fh.close()
            return [tmpdir,aln]
        elif os.path.isfile(pdbdir+'needcombinewithbase'):
            basefile=pdbdir+codelist[1]+'.base.pdb'
            if not os.path.isfile(basefile):
                basefile=pdbdir+codelist[1]+'.base'
            if not os.path.isfile(basefile):
                basefile=pdbdir+codelist[0]+'.'+codelist[1]+'.base.pdb'
            if not os.path.isfile(basefile):
                basefile=pdbdir+codelist[1]+'.base.gz'
                fh=gzip.open(basefile)
            else:
                fh=open(basefile,'r')
            fc=fh.read()
            fh.close()
            fh=gzip.open(pdbdir+alncode+'.pdb.gz','r')
            dfc=fh.read()
            fh.close()
            fc=fc.replace('replacethis',dfc)
            fh=open(os.path.join(tmpdir,alncode+'.pdb'),'w')
            fh.write(fc)
            fh.close()
            return [tmpdir,aln]
        elif os.path.isfile(pdbdir+'needtransformation'):
            tranformationmatrix=aln[0].name
            tm=tranformationmatrix.split(" ")
            tm=[float(item) for item in tm]
            if not aln[0].code.startswith(self.currentdsname):
                env.libs.topology.read(file='$(LIB)/top_heav.lib')
                env.libs.parameters.read(file='$(LIB)/par.lib')
                self.currentdsname=codelist[0]+'.'+codelist[1]
                #print self.currentdsname
                self.basemodel=model(env,file=pdbdir+self.currentdsname+'.base.pdb')
                taln=alignment(env)
                taln.append_model(self.basemodel,align_codes='original')
                taln.append_model(self.basemodel,align_codes='tbt')
                self.basealn=taln
                self.currentmodel=model(env)
                self.currentmodel.generate_topology(taln[0])
                if len(self.basemodel.chains)>2:
                    raise Excpetion('more than 2 chains in base pdb')
            self.currentmodel.transfer_xyz(self.basealn)
            sel=selection(self.currentmodel.chains[1])
            sel.rotate_origin((1,0,0),tm[0]*57.29577951308232)
            sel.rotate_origin((0,1,0),-tm[1]*57.29577951308232)
            sel.rotate_origin((0,0,1),tm[2]*57.29577951308232)
            sel.translate(tm[3:6])
            naln=alignment(env)
            naln.append_model(self.currentmodel,align_codes='formdt')
            return [pdbdir,naln]
        elif 0 and isinstance(pdbdir, list):
            for spdbdir in pdbdir:
                if os.path.isfile(os.path.join(tmpdir,code+'.pdb')) or os.path.isfile(os.path.join(tmpdir,'pdb'+code+'.pdb')) or os.path.isfile(os.path.join(tmpdir,'pdb'+code[:4]+'.pdb')):
                    return [tmpdir,aln]
                elif os.path.isfile(os.path.join(spdbdir,code+'.pdb.gz')):
                    print os.system('gzip -d -c  '+os.path.join(spdbdir,code+'.pdb.gz')+' > '+os.path.join(tmpdir,code+'.pdb'))
                    return [tmpdir,aln]
                elif os.path.isfile(os.path.join(spdbdir,'pdb'+code+'.pdb.gz')):
                    print os.system('gzip -d -c  '+os.path.join(spdbdir,'pdb'+code+'.pdb.gz')+' > '+os.path.join(tmpdir,'pdb'+code+'.pdb'))
                    return [tmpdir,aln]            
                elif os.path.isfile(os.path.join(spdbdir,'pdb'+code[:4]+'.pdb.gz')):
                    print os.system('gzip -d -c  '+os.path.join(spdbdir,'pdb'+code[:4]+'.pdb.gz')+' > '+os.path.join(tmpdir,'pdb'+code[:4]+'.pdb'))
                    return [tmpdir,aln]   
        #elif os.path.isfile(os.path.join(tmpdir,code+'.pdb')) or os.path.isfile(os.path.join(tmpdir,'pdb'+code+'.pdb')) or os.path.isfile(os.path.join(tmpdir,'pdb'+code[:4]+'.pdb')):
        #    return [tmpdir,aln]
        #elif os.path.isfile(os.path.join(pdbdir,code+'.pdb.gz')):
        #    print os.system('gzip -d -c  '+os.path.join(pdbdir,code+'.pdb.gz')+' > '+os.path.join(tmpdir,code+'.pdb'))
        #    return [tmpdir,aln]
        #elif os.path.isfile(os.path.join(pdbdir,'pdb'+code+'.pdb.gz')):
        #    print os.system('gzip -d -c  '+os.path.join(pdbdir,'pdb'+code+'.pdb.gz')+' > '+os.path.join(tmpdir,'pdb'+code+'.pdb'))
        #    return [tmpdir,aln]            
        #elif os.path.isfile(os.path.join(pdbdir,'pdb'+code[:4]+'.pdb.gz')):
        #    print os.system('gzip -d -c  '+os.path.join(pdbdir,'pdb'+code[:4]+'.pdb.gz')+' > '+os.path.join(tmpdir,'pdb'+code[:4]+'.pdb'))
        #    return [tmpdir,aln]           
        else:
            return [pdbdir,aln]

    def reshapesp(self): #should be useless after all the data are upgraded
        if runenv.hostn:
            return 0
        f=h5py.File(self.sppath+'.hdf5','r')
        ms=f['mdt'].shape
        f.close()
        if self.bl<ms[self.bp]:
            print "reshaping original mdt"
            env = environ()
            log.minimal()
            mlib=mdt.Library(env)
            featurelist=feature(self.features,mlib).get_featurelist()
            m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
            m.read_hdf5(self.sppath+'.hdf5')
            fl=len(featurelist)
            #m=m.reshape(featurelist,[0 for i in range(0,fl)], [-1 for i in range(0,fl)])
            m.write_hdf5(self.sppath+'.hdf5')

    def make_symetry(self,mdtb):
        ms=mdtb.shape
        if feature(self.features).issymetry():
            if (mdtb[1,0]==mdtb[0,1]).all():
                return 0
            for i in range(0,ms[0]):
                for j in range(i+1,ms[1]):
                    mdtb[i,j,...]=mdtb[i,j,...]+mdtb[j,i,...]
            for i in range(0,ms[0]):
                for j in range(0,i):
                    mdtb[i,j,...]=mdtb[j,i,...]          

    def get_stretch_dictionary(self,pm=''):
        sd=[]
        kwl=[]
        if re.search('([0-9\.]{2,10})',pm):
            rer=re.search('([0-9\.]{2,10})',pm)
            l=float(rer.group(1))
        else:
            l=0.59
        for i in range(len(self.fo.fclist)):
            if self.fo.fclist[i].startswith('d'):
                if pm.endswith('u'):
                    sd.append([{'type':'uniform','value':(self.fo.frlist[i]/self.fo.fdlist[i])/(l*1.414)}])                    
                else:
                    sd.append([{'type':'lowpass','turnpoint':4/(self.fo.frlist[i]/self.fo.fdlist[i]),'turnpower':2,'scale':10},
                    {'type':'lowcount','turnpoint':30,'turnpower':1,'xkernelwidth':[0.1/(self.fo.frlist[i]/self.fo.fdlist[i])],'lowvalue':0.01},
                    {'type':'uniform','value':(self.fo.frlist[i]/self.fo.fdlist[i])/(l*1.414)}]) #l is the actual kernel width
                kwl.append(0.2/(self.fo.frlist[i]/self.fo.fdlist[i]))    
            elif self.fo.fclist[i].startswith('h') or self.fo.fclist[i].startswith('g'):
                if self.fo.fclist[i].startswith('h'):
                    pp=True
                else:
                    pp=False
                sd.append([{'type':'uniform','value':(self.fo.frlist[i]/self.fo.fdlist[i])/(10*l*1.414),'periodic':pp}])
                kwl.append(10/(self.fo.frlist[i]/self.fo.fdlist[i]))
            elif 0:#self.fo.fclist[i].startswith('l'):
                sd.append([{'type':'lowpass','turnpoint':1/(self.fo.frlist[i]/self.fo.fdlist[i]),'turnpower':8,'scale':1000},
                {'type':'lowcount','turnpoint':30,'turnpower':1,'xkernelwidth':[0.1/(self.fo.frlist[i]/self.fo.fdlist[i])],'lowvalue':0.1},
                {'type':'uniform','value':(self.fo.frlist[i]/self.fo.fdlist[i])/(l*1.414)}]) #l is the actual kernel width
                kwl.append(10/(self.fo.frlist[i]/self.fo.fdlist[i]))                  
            elif self.fo.fclist[i][0] not in ['a','r','t']:
                sd.append([{'type':'uniform','value':(self.fo.frlist[i]/self.fo.fdlist[i])/(l*1.414)}])
                kwl.append(10/(self.fo.frlist[i]/self.fo.fdlist[i])) #this should not be used as we don't know much....
        return sd,kwl

    def ksmoothing(self,ma):
        gps=gpsmoothing()
        sd,kwl=self.get_stretch_dictionary()
        ma=np.transpose(ma,self.rs1)
        ms=ma.shape    
        dvdim=1
        for i in range(0,self.ifl):
            dvdim=dvdim*ms[i]
        ma=np.reshape(ma,(dvdim,-1))
        for i in range(0,dvdim):
            print i
            ma[i]=gps.ksmoothing_single_withscale(ma[i],sd)
        ma=np.reshape(ma,ms) 
        if self.permute:
            ma=np.transpose(ma,self.rs2)
        return ma
    
    def gpsmoothing(self,ma):        
        gps=gpsmoothing()
        #define prior
        #prior_mean=self.get_prior(ma)
        sd,kwl=self.get_stretch_dictionary()
        if self.permute:
            ma=np.transpose(ma,self.rs1)
        ms=ma.shape    
        dvdim=1
        for i in range(0,self.ifl):
            dvdim=dvdim*ms[i]
        ma=np.reshape(ma,(dvdim,-1))
        for i in range(0,dvdim):
            print i
            ma[i]=gps.gpsmoothing_single(ma[i],sd,kwl)
        ma=np.reshape(ma,ms) 
        if self.permute:
            ma=np.transpose(ma,self.rs2)
        return ma

    def get_typeindependent_distribution(self):
        ma=self.load2npy(self.sppath+'.hdf5')
        if self.permute:
            ma=np.transpose(ma,self.rs1)
        ms=ma.shape    
        self.make_symetry(ma)         
        dvdim=1
        for i in range(0,self.ifl):
            ma=ma.sum(axis=0)
        return ma
    
    def get_prior(self,mdta):
        pd=np.array(1,dtype=np.float32)
        ind1=[]
        ind2=[]
        for i in range(len(self.fo.fclist)):
            if i in self.rs1[0:self.ifl]:
                continue
            if self.fo.fclist[i].startswith('d'):
                increasefunction=True
            else:
                increasefunction=False
            #ma=self.get_typeindependent_distribution()
            ma=mdta
            il=range(len(self.fo.fclist))
            il.reverse()
            for j in il:
                if j!=i:
                    ma=ma.sum(axis=j)
            ma=ma/ma.sum()
            malog=np.log(ma)
            malogmin=malog[malog>-100000].min()
            malog[malog<-1000]=malogmin-10
            #maindex=malog>-1000
            #malog=malog[maindex]
            nmalog=np.zeros(len(ma))
            if len(ma)>12: #means distance...
                std=0.4/(self.fo.frlist[i]/self.fo.fdlist[i])
                nd=scipy.stats.norm.pdf(range(-len(ma)+1,0)+range(len(ma)),loc=0,scale=std)
                for k in range(len(ma)):
                    si=len(ma)-k-1
                    ei=2*len(ma)-k-1
                    nmalog[k]=(nd[si:ei]*malog).sum()/nd[si:ei].sum()
                spd=np.exp(nmalog)
                if 0:
                    def errfunc(scale):
                        scale=np.abs(scale)
                        us=scipy.interpolate.InterpolatedUnivariateSpline(range(0,int(self.fo.frlist[i])),np.abs(scale))
                        igd=us(self.bins[0])
                        igd=igd/sum(igd)
                        err=10*((ma-igd)**2).sum()+((malog-np.log(np.abs(igd))[maindex])**2).sum()
                        if increasefunction:
                            err=err+((scale[2:]-scale[1:-1])<0).sum()
                        return err
                    bw=self.fo.fdlist[0]/self.fo.frlist[0]
                    optres=scipy.optimize.fmin_powell(errfunc,ma[0::bw])
                    print optres
                    scale=np.abs(optres)
                    us=scipy.interpolate.InterpolatedUnivariateSpline(range(0,int(self.fo.frlist[0])),scale)
                    spd=us(self.bins[0])
            else:
                spd=ma
            pd=pd.reshape(ind1+[1])*(spd.reshape(ind2+[len(ma)]))
            ind1.append(len(ma))
            ind2.append(1)
        pd[...]=pd[...]/pd[...].sum()
        return pd

    def prepare_atompair_clustering(self):
        a1=self.load2npy(self.sppath+'.hdf5')
        a1[...]=a1+0.1
        ndim=self.fo.fdlist
        fnl=feature('a158').get_feature_names()
        data=np.zeros([ndim[1]*(ndim[1]+1)/2,ndim[0]])
        dnl=[]
        apa=np.zeros([ndim[1],ndim[1]],dtype=np.int32)-1
        k=0
        for i in range(ndim[1]):
            for j in range(i,ndim[1]):
                #pdb.set_trace()
                data[k,:]=a1[:,i,j]+a1[:,j,i]
                apa[i,j]=k
                apa[j,i]=k
                dnl.append(fnl[i]+'-'+fnl[j])
                k+=1
        for i in range(ndim[1]*(ndim[1]+1)/2):
            data[i,:]=data[i,:]/max(np.mean(data[i,:]),0.00000000000001)
        return [a1,data,dnl,apa,ndim]

    def generate_atompair_clusters(self,a1,data,dnl,apa,ndim,dist,linkage,r,c,minmember=4,libname='',figtype='.png',mincluster=5):
        apaoriginal=np.copy(apa)
        os.chdir(self.dir)
        if not os.path.isdir('apcluster'):
            os.mkdir('apcluster')
        os.chdir('apcluster')
        #clustering the distributions
        if not os.path.isfile(dist+'.npy'):
            da=scipy.spatial.distance.pdist(data,dist)
            np.save(dist,da)
        else:
            da=np.load(dist+'.npy')
        if not os.path.isfile(dist+'-'+linkage+'.npy'):
            la=scipy.cluster.hierarchy.linkage(da,linkage)
            np.save(dist+'-'+linkage+'.npy',la)
        else:
            la=np.load(dist+'-'+linkage+'.npy')
        noc=r*c
        type=dist+'-'+linkage+str(noc)+'-'+str(minmember)
        print type
        if not os.path.isfile(type+'.npy'):
            l1=scipy.cluster.hierarchy.fcluster(la,noc,'maxclust')
            nl=[]
            for i in range(1,noc+1):
                if (l1==i).sum()<minmember:
                    l1[l1==i]=9999
                else:
                    nl.append(i)   
            for i in range(len(nl)):
                l1[l1==nl[i]]=i
            l1[l1==9999]=len(nl)
            apa[...]=l1[apa]    
            hc=np.histogram(l1,noc)
            plt.plot(hc[0])
            plt.savefig('histogram.eps')
            plt.clf()
            if apa.max()!=l1.max() or apa.min()!=0:
                raise Exception('bugs in code')
            np.save(type+'.npy',apa)
        else:
            apa=np.load(type+'.npy')
        if not os.path.isfile(type+figtype) and apa.max()>5:
            cnl=[[] for i in range(noc)]
            plt.tick_params(labelsize=5)
            for i in range(noc):
                plt.subplot(r,c,i+1)
                print 'plotting '+type+str(i)
                for j in range(ndim[1]):
                    for k in range(i,ndim[1]):
                        if apa[j,k]==i:
                            plt.plot(np.arange(0,15,0.1),data[apaoriginal[j,k]])
                            #cnl[i].append(dnl[j])
                            plt.title(str(i),fontsize=7)
                            plt.xticks([])
                            plt.yticks([])
                            plt.title('max '+str(apa.max()))
            plt.tick_params(labelsize=5)
            print apa.max()
            print apa.min()
            plt.savefig(type+figtype)
            plt.clf()
        if len(libname)>0:
            print os.system('mv '+type+'.npy '+runenv.libdir+'/'+libname)
            print os.system('scp '+runenv.libdir+'/'+libname+' '+runenv.jobserver+':~/lib/')
        
class scaledsp(rawsp):
    def __init__(self, pm='',rspo=[],model=[],previouspm=''):
        if model:
            ns=copy.deepcopy(model)
            rspo=rawsp(model=ns)
            pm=model['pm']
        self.opm=pm
        self.pm=pm.replace('+','.')
        if 'i1' in self.pm:
            self.ifl=[0]
            #self.pm=self.pm[:-2]
        if rspo!=[]:
            self.rsp=rspo
            rawsp.__init__(self,rspo.pdbset,rspo.features,rspo.genmethod,tlib=rspo.tlib)
            if previouspm:
                self.sppath=self.sppath+'.'+previouspm
                self.spname=self.spname+'.'+previouspm
            self.ssppath=self.sppath+'.'+pm.replace('+','.')
            self.sspname=self.spname+'.'+pm.replace('+','.')
            self.array=[]
            self.justwait=False
        self.nors=400

    def get_task(self):
        pml=self.opm.split('+')
        taskchains=[]
        ssp=scaledsp(pm=pml[0],rspo=self.rsp)
        if not os.path.isfile(ssp.ssppath+'.hdf5'):
            taskchains.append(ssp.get_scale_task())
        if len(pml)>1:
            for i in range(1,len(pml)):
                nssp=scaledsp(pm=pml[i],rspo=self.rsp,previouspm='.'.join(pml[0:i]))
                if not os.path.isfile(nssp.ssppath+'.hdf5'):
                    taskchains.append(nssp.get_scale_task())
        if len(taskchains)>0:
            return taskchain(chains=taskchains)
        else:
            return 0
        
    def get_scale_task(self):
        if self.pm.startswith('gps') or self.pm.startswith('ks') or self.pm.startswith('sks'):
            self.task=task('','',afterprocessing=self,preparing=self,targetfile=self.dir+self.sspname+'.hdf5')
            self.rundir=self.sspname
            self.runpath='/netapp/sali/gqdong/'+self.rundir+'/'
            self.task.dir=self.dir+self.rundir+'/'
            self.task.rdirname=self.rundir
            return self.task
        else:
            return localtask(func=self.scalesp,inputlist=[])

    def prepare_task_input(self):
        #prepare the run directories
        os.chdir(self.dir)
        #os.system('rm -r '+self.rundir)      
        os.mkdir(self.sspname)
        os.chdir(self.dir+self.sspname)
        #save the input array
        ma=self.scalesp(pm='getma')
        ms=ma.shape
        trn=ms[0]*ms[1]*ms[1]/smoothingjobratio
        self.nors=max(min(int(trn),self.nors),1)
        nofrpn=np.round(float(ms[0])/self.nors)
        #self.nors=int(np.round(float(ms[0])/nofrpn+0.5))
        for i in range(int(self.nors)):
            if i==self.nors-1:
                ep=ms[0]
            else:
                ep=nofrpn*(i+1)
            sma=ma[nofrpn*i:ep]
            np.save(str(i+1),sma)
        gps=gpsmoothing()
        gps.nofrpn=nofrpn
        gps.lib=self.tlib
        gps.stretchdict,gps.priorkw=self.get_stretch_dictionary(self.pm) 
        gps.sm=self.pm
        fh=open('input.pickle','w')
        pickle.dump(gps,fh)
        fh.close()
        if self.pm.startswith('gps'):
            freememory=0.55+8*(ms[1]**2)/1000000000.0
            if freememory<1.9:
                freememory=1.9
        else:
            freememory=0.8
        nors=self.nors
        #write the input list for different runs, maybe not necessary
        inputlist=open('inputlist','w')
        inputlist.write(','.join(['runme.py' for i in range(1,nors+1)]))
        inputlist.close()
        os.system('cp '+runenv.scriptpath+' ./')
        makemdt=open('runme.py','w')
        makemdt.write('keeptry=5\nwhile keeptry:\n    try:\n        from SOAP.statsTable import *\n        keeptry=0\n    except:\n        keeptry-=1\n'
                      +'import sys \nso=gpsmoothing()\n'
                      +'so.runtask_cluster(sys.argv[1])')
        makemdt.flush()
        generate_job_submit_script(freememory,self.rundir,runtime,nors)
        return [] 
    
    def afterprocessing(self):
        os.chdir(self.dir+self.rundir)
        nal=[]
        for i in range(1,self.nors+1):
            nal.append(np.load(str(i)+'.smoothed.npy'))
        nma=np.vstack(nal)
        self.scalesp(pm='savema',sma=nma)
        print os.system('rm -r '+self.dir+self.sspname)
        return 1

    def scalesp(self,pm='',ssppath='',ifpos=None,nbpos=None, ref=None,sma=[]):
        #self.reshapesp()
        if pm=='':
            pm=self.pm
        if not pm:
            return 0
        if not ifpos:
            ifpos=self.ifl
        if not ssppath:
            ssppath=self.ssppath
        if os.path.isfile(ssppath+'.hdf5'):
            return 0
        if pm.startswith('pp'):
            print os.system('cp '+self.sppath+'.hdf5 '+ssppath+'.hdf5')
            f2=h5py.File(ssppath+'.hdf5','r+')
            mdtc=f2['mdt']
            mdtc[...]=mdtc[...]+float(pm[2:])
            f2.flush()
            f2.close()
            return 1
        elif pm.startswith('exp'):
            print os.system('cp '+self.sppath+'.hdf5 '+ssppath+'.hdf5')
            f2=h5py.File(ssppath+'.hdf5','r+')
            mdtc=f2['mdt']
            mdtc[...]=np.exp(mdtc[...])
            f2.flush()
            f2.close()
            return 1
        elif pm.startswith('svd'):
            s,v=self.get_svd()
            si=int(pm[3])
            pi=int(pm[4:])
            diagElem=[int((i>si and i<pi) and 1) for i in range(si.shape[1])]
            f2=h5py.File(ssppath+'.hdf5','r+')
            mdtc=f2['mdt']
            pa=self.get_SVD_subset(mdtc.shape,diagElem)
            mdtc[...]=pa[...]
            f2.flush()
            f2.close()
            return 1                        
        f=h5py.File(self.sppath+'.hdf5','r')
        mdtb=f['mdt']
        ms=mdtb.shape
        ma=np.zeros(ms)
        ma[...]=mdtb[...]
        f.close()
        if self.permute:
            ma=np.transpose(ma,self.rs1)
        ms=ma.shape   
        if pm!='savema':          
            self.make_symetry(ma)         
            dvdim=1
            for i in range(0,ifpos):
                dvdim=dvdim*ms[i]
            ma=np.reshape(ma,(dvdim,-1))
            #do the actual scaling or other...
            if pm=='getma': #return ma 
                return ma
            else:
                if pm.startswith('ks') or pm.startswith('gps'):
                    gps=gpsmoothing()
                    gps.sm=self.pm
                    gps.lib=self.tlib
                    gps.bincenters=self.bins
                    gps.stretchdict,gps.priorkw=self.get_stretch_dictionary(self.pm) 
                else:
                    nolog, sums,sume,ref=self.decode_pm(pm,ma)            
                for i in range(0,dvdim):
                    if pm.startswith('ks'):
                        ma[i]=gps.ksmoothing_single_withscale(ma[i])
                    elif pm.startswith('gps'):
                        gps.ind=i
                        ma[i]=gps.gpsmoothing_single(ma[i])
                    else:
                        self.applysf_single(nolog, ma[i],ref,sums, sume)
        else:
            ma=sma
        ma=np.reshape(ma,ms) 
        if self.permute:
            ma=np.transpose(ma,self.rs2)
        print os.system('cp '+self.sppath+'.hdf5 '+ssppath+'.hdf5')
        f2=h5py.File(ssppath+'.hdf5','r+')
        mdtc=f2['mdt']
        mdtc[...]=ma[...]
        f2.flush()
        f2.close()
        return 1
    
    def get_sf(self,sfm,ma):
        ref=sf(sftype=sfm, data=ma).get_sf()
        ref=np.exp(ref)
        return np.squeeze(ref)

    def get_SVD_subset(self,ms,ratio):
        s,v=self.get_svd()
        pa=np.dot(np.dot(s,np.eye(ratio)),v)
        return pa[...].reshape(ms)        
   
    def decode_pm(self,pm,ma):
        if pm.endswith('i1'):
            pm=pm[:-2]
        if pm.endswith('nolog'):
            nolog=True
            pm=pm[:-5]
        else:
            nolog=False
        if re.search('-',pm):
            pms=pm.split('-')
            ref=self.get_sf(pms[0],ma)
            nmm=pms[1]
            sums,sume=self.decode_nmm(nmm,ma)
            self.normalize_sinlge_dist(ref,sums,sume) # the reference is always normalized using the same mehtod            
        else:
            ref=[]
            nmm=pm
            sums,sume=self.decode_nmm(nmm,ma)
        return nolog, sums,sume,ref
   
    def decode_nmm(self,nmm,ma):
        mas=ma.shape
        if re.search('nb([0-9]*)to([0-9]*)',nmm):
            rer=re.search('nb([0-9]*)to([0-9]*)',nmm) 
            sume=int(rer.group(2))+1
            sums=int(rer.group(1))
        elif nmm=='nbsum':
            sume=ma.shape[1]
            sums=0
        elif nmm=='npend':
            sume=ma.shape[1]
            sums=sume-1
        elif nmm.startswith('np'):
            sums=int(nmm[2:])
            sume=sums+1
        return sums,sume
        
    def applysf_single(self,nolog, mai,ref,sums, sume):
        self.normalize_sinlge_dist(mai,sums,sume,nolog)
        if len(ref)>0:
            mai[...]=mai/ref
        if not nolog:
            self.log_nozero(mai)
            
    def log_nozero(self,sd):
        sd[...]=sd+(sd==0)*np.exp(-20.123456789)
        sd[...]=-np.log(sd)
    
    def normalize_sinlge_dist(self,sd,sums,sume,nolog):
        sdsum=sd[sums:sume].mean()
        if sdsum==0:
            print "Can not normalize array by sum because of zero:"+self.pm
            return 0
        sd[...]=sd/(sdsum)
        
    def load(self):
        if os.path.isfile(self.ssppath+'.npy'):
            self.loadnpz()
            return 0
        if not os.path.isfile(self.ssppath+'.hdf5'):
            self.scalesp()
        f=h5py.File(self.ssppath+'.hdf5','r')
        mdtb=f['mdt']
        ms=mdtb.shape
        print ms
        mdta=np.zeros(np.array(list(ms)))
        mdta[...]=mdtb[...]
        f.close()
        self.array=mdta
        return 0
        
    def savenpz(self):
        if len(self.array)==0:
            self.load()
        np.save( self.ssppath,self.array)
        
    def loadnpz(self):
        self.array=np.load(self.ssppath+'.npy')
        
    def add_ref(self,ref,refname='',savetype='',apca=[],ratio=[],svdilist=[],permute=False):
        #normalize raw distribution by making the value at pos to be 1. This version is only suitable for distance distribution.
        if len(self.array)==0:
            self.load()
        #pdb.set_trace()
        mdtb=np.copy(self.array)
        ms=mdtb.shape
        if self.permute:
            mdtb=np.transpose(mdtb,self.rs1)
        if len(svdilist)>0:
            ratiorlist=[0 for i in range(ms[0])]
            for sr,si in zip(ratio,svdilist):
                ratiolist[si]=sr
            mdtb=self.get_SVD_subset(ms,ratiolist)
        elif len(apca)==0:
            rl=ref.shape
            refshape=[1 for i in range(0,len(ms)-len(rl))]+list(rl)
            ref.resize(refshape)            
            mdtb=(mdtb+ref)*ratio[0]
        else:
            rl=ref[0].shape
            refshape=[1 for i in range(0,len(ms)-len(rl))]+list(rl)
            for sref in ref:
                try:
                    sref.resize(refshape)
                except:
                    pdb.set_trace()
            for i in range(apca.shape[0]):
                for j in range(apca.shape[1]):
                    mdtb[i,j,...]=(mdtb[i,j]+ref[apca[i,j]])*ratio[apca[i,j]]
        #pdb.set_trace()
        if len(svdilist)==0 and self.permute:
            mdtb=np.transpose(mdtb,self.rs2) 
        if savetype=='hdf5':
            self.write_hdf5(self.ssppath+'.'+refname+'.hdf5',mdtb,permute=permute)
        elif savetype=='lib':
            self.write_lib(self.ssppath+'.'+refname,mdtb)
        elif savetype=='npy':
            np.save(self.ssppath+'.'+refname+'.npy',mdtb)
        return mdtb
    
    def write_hdf5(self,path, mdtarray, permute=False):
        env = environ()
        env.libs.topology.read('${LIB}/top_heav.lib')
        env.libs.parameters.read('${LIB}/par.lib')        
        mlib=mdt.Library(env)
        mlib.atom_classes.read('$(LIB)/atmcls-mf.lib')
        mlib.bond_classes.read('${LIB}/bndgrp.lib')
        featurelist=feature(self.features,mlib).get_featurelist()
        m=mdt.Table(mlib,features=featurelist,bin_type=mdt.Float,shape=[-1]*len(featurelist))
        mdl = model(env)
        mdl.build_sequence('CCG/CCG')
        aln = alignment(env)
        aln.append_model(mdl, align_codes='tmp')
        csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(self.genmethod)        
        m.add_alignment(aln,chain_span_range=(-csmax,-csmin,csmin,csmax),residue_span_range=(-kcut1,-kcut,kcut,kcut1),bond_span_range=(bs,bsm))#,startpos=0,endpos=0,scale=0,refyes=False,refmdt='',io=io)
        ma=m.get_array_view()
        ma[...]=mdtarray
        if permute:
            d,a2,a1=featurelist
            m = m.reshape(features=(a1,a2,d), offset=[0,0,0], shape=m.shape[::-1])
        m.write_hdf5(path,gzip=True)
        
    def write_lib(self,path=[],mdtarray=[]):
        env = environ()
        log.minimal()
        mlib=mdt.Library(env)
        featurelist=feature(self.features,mlib).get_featurelist()
        m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
        dim=m.shape
        if self.bl==dim[self.bp]:
            dim=np.array(dim)+1
        entry_header=[0.0,0.75,14.75,0.5,0.0,0.0]
        sp=float(self.fo.fslist[0])
        ep=float(self.fo.frlist[0])
        nofb=int(self.fo.fdlist[0])
        binsize=(ep-sp)/nofb
        entry_header[3]=binsize
        entry_header[1]=binsize*1.5
        entry_header[2]=ep-binsize*0.5
        header=[10,29,1,31,2,35,0]
        header[1]=nofb-1
        header[-2]=header[1]+6
        headerstrlist=["%5d" % item for item in header]
        header='R'+''.join(headerstrlist)
        if path:
            potfile=open(path+'.lib','w')
            print path+'.lib'
        else:
            potfile=open(self.ssppath+'.lib','w')
        if mdtarray!=[]:
            mdtt=mdtarray
        else:
            self.load()
            mdtt=self.array
        if self.permute:
            mdtt=np.transpose(mdtt,self.rs1)
        potfile.write('MOD5')
        potfile.write('\n')
        symbollist=[]
        fp=self.rs1[0]
        for typ1 in range(dim[fp]-1):
            symbollist.append(m.features[1].bins[typ1].symbol.ljust(4))
        for typ1 in range(dim[fp]-1):
            for typ2 in range(typ1,dim[fp]-1):
                # entry header
                entryheader = header+'   '+symbollist[typ1]+'    '+symbollist[typ2]
                potfile.write(entryheader)
                for number in entry_header:
                    potfile.write("%10.4f" % number)
                for dist in range(1,dim[self.rs1[-1]]-1):
                    potfile.write("%10.4f" % mdtt[typ1,typ2,dist] )
                potfile.write('\n')

    def get(self):
        if not os.path.isfile(self.ssppath+'.hdf5'):
            self.scalesp()
        else:
            return 0
    
    def poscutoff(self,cutoff):
        refname='pc'+str(cutoff)
        print os.system('cp '+self.ssppath+'.hdf5 '+self.ssppath+'.'+refname+'.hdf5')
        f=h5py.File(self.ssppath+'.'+refname+'.hdf5','r+')
        mdtb=f['mdt']
        mdtb[:,:,cutoff:]=0
        f.flush()
        f.close()

    def get_svd(self):
        if os.path.isfile(self.ssppath+'svdv.npy'):
            return [np.load(self.ssppath+'svdu.npy'),np.load(self.ssppath+'svdv.npy')]
        mdta=self.load2npy(self.ssppath+'.hdf5')
        k=0
        #for fc in self.fo.fclist:
        #    if fc.startswith('d'):
        #        break
        #    k+=1
        if k!=0:
            raise Exception('can only run svd for distance with index 0 at the moment')
        nd=1
        for i in range(1, len(mdta.shape)):
            nd*=mdta.shape[i]
        u,s,v=numpy.linalg.svd(mdta.reshape([mdta.shape[0],nd]),full_matrices=0)
        us=np.dot(u,np.diag(s))
        np.save(self.ssppath+'svdu.npy',us)
        np.save(self.ssppath+'svdv.npy',v)
        return [us,v]

class gpsmoothing(object):
    def __init__(self):
        self.data_variance_power=-1
        self.prior_covariance_turnpoint=[160.0]
        self.prior_covariance_turnpower=[4.0]
        self.prior_covariance_stretchratio=[10]
        self.prior_covariance_cutoff=30
        self.prior_baseline_variance=1
        self.prior_covariance_correlationlegth_inv2=[0.001]
        self.data_baseline_variance=1
        #data
        self.raw_data=[]
        self.data_mean=[]
        self.data_covariance=[]
        self.data_covariance_inv=[]
        self.prior_mean=[]
        self.prior_covariance=[]
        self.prior_covariance_inv=[]
        self.posterior_covariance=[]
        self.posterior_mean=[]
        self.stretchdict=[]
        self.priorkw=[]
        self.sm='gps'
        self.nofrpn=0
        self.lib=''
        self.ind=None
        self.apo=None
        self.bincenters=[]

    def prepare_mean(self,data):
        ndata=np.array(data)
        ndata=data.astype(np.float32)/data.sum()
        ndata.resize([ndata.size,1])
        #ndata[np.isnan(ndata)]=0
        #ndata[ndata>1000]=0
        #ndata[ndata<-1000]=0
        return ndata
    
    def calc_prior_mean(self):
        self.prior_mean=self.prepare_mean(self.prior_mean)
    
    def get_kernel_mask(self,ksize,width=[1]):
        ka=np.array(1)
        ind1=[]
        ind2=[]
        sumarray=np.zeros(ksize)
        for i in range(len(ksize)):
            ds=ksize[i]
            dw=width[i]
            nd=scipy.stats.norm.pdf(range(-ds+1,0)+range(ds),loc=0,scale=dw)
            ka=ka.reshape(ind1+[1])*(nd.reshape(ind2+[ds*2-1]))
            ind1.append(ds*2-1)
            ind2.append(1)
        if len(ksize)==1:
            for i in range(ds):
                si=ds-i-1
                ei=2*ds-i-1
                sumarray[i]=ka[si:ei].sum() #the sum array for each data  point
        elif len(ksize)==4:
            for i in range(ksize[0]):
                s0=ksize[0]-i-1
                e0=2*ksize[0]-i-1
                for j in range(ksize[1]):
                    s1=ksize[1]-j-1
                    e1=2*ksize[1]-j-1
                    for k in range(ksize[2]):
                        s2=ksize[2]-k-1
                        e2=2*ksize[2]-k-1
                        for l in range(ksize[3]):
                            s3=ksize[3]-l-1
                            e3=2*ksize[3]-l-1
                            sumarray[i,j,k,l]=ka[s0:e0,s1:e1,s2:e2,s3:e3].sum()
        return [ka,sumarray]
    
    def ksmoothing_single_usingmask(self,data,ka,sumarray):
        nd=np.zeros(data.shape)
        if len(data.shape)==1:
            for i in range(len(data)):
                si=len(data)-i-1
                ei=2*len(data)-i-1
                nd[i]=(ka[si:ei]*data).sum()/sumarray[i]
        elif len(data.shape)==4:
            ksize=data.shape
            for i in range(ksize[0]):
                s0=ksize[0]-i-1
                e0=2*ksize[0]-i-1
                for j in range(ksize[1]):
                    s1=ksize[1]-j-1
                    e1=2*ksize[1]-j-1
                    for k in range(ksize[2]):
                        s2=ksize[2]-k-1
                        e2=2*ksize[2]-k-1
                        for l in range(ksize[3]):
                            s3=ksize[3]-l-1
                            e3=2*ksize[3]-l-1
                            nd[i,j,k,l]=(ka[s0:e0,s1:e1,s2:e2,s3:e3]*data).sum()/sumarray[i,j,k,l]
        return nd
    
    def ksmoothing_single_uniform(self,data,width):
        ksize=data.shape
        axisscale=[]
        for i in range(len(ksize)):
            ds=ksize[i]
            dw=width[i]
            axisscale.append(np.arange(1.0,ds+1.0)*(1.0/dw))
        na=self.ksmoothing_single_usingscale(data, axisscale)
        return na

    def get_stretched_axis(self,data,stretchdict):
        #calc the strecthed axis to taking into account low count, and other things
        sal=[]
        self.periodic=np.zeros(len(data.shape))
        for i in range(len(data.shape)):
            sd=stretchdict[i]
            sddata=data
            for j in range(len(data.shape)-1,-1,-1):
                if j!=i:
                    sddata=sddata.sum(axis=j)
            sddata=np.squeeze(sddata)
            sas=np.ones(len(sddata))
            for d in sd:
                if 'periodic' in d and d['periodic']==True:
                    self.periodic[i]=1
                if d['type']=='lowpass':
                    sas=sas*(d['scale']*(d['turnpoint']**d['turnpower'])/(d['turnpoint']**d['turnpower']+np.arange(len(data),dtype=np.float32)**d['turnpower'])+1)
                #self.ksmoothing_single_uniform(sddata,)
                elif d['type']=='lowcount':
                    sdd=self.ksmoothing_single_uniform(sddata,d['xkernelwidth'])
                    sddr=np.sqrt(sdd)
                    sas=sas*(sdd**d['turnpower']/(sdd**d['turnpower']+d['turnpoint']**d['turnpower'])+d['lowvalue'])
                elif d['type']=='uniform':
                    sas=sas*d['value']
            sal.append(sas.cumsum())
        return sal
    
    def get_kernel_mask_withscale_singlepos(self,ksize,scale=[],pos=[]):
        #get the weight array around a single point, so we can do smoothing by weighted sum over the values
        ka=np.array(1)
        ind1=[]
        ind2=[]
        for i in range(len(ksize)):
            ds=ksize[i]
            if self.periodic[i]:
                points=self.get_periodic_points(scale[i],pos[i])
            else:
                points=scale[i]-scale[i][pos[i]]
            nd=scipy.stats.norm.pdf(points,loc=0.0,scale=1.0)
            ka=ka.reshape(ind1+[1])*(nd.reshape(ind2+[len(nd)]))
            ind1.append(len(nd))
            ind2.append(1)
        return ka
        
    def ksmoothing_single_usingscale(self,data,axisscale):
        na=np.zeros(data.shape)
        for i in range(data.size):
            ci=np.unravel_index(i,data.shape)
            ka=self.get_kernel_mask_withscale_singlepos(data.shape,axisscale,ci)
            na[ci]=(data*ka).sum()/ka.sum()
        return na
    
    def ksmoothing_single_withscale(self,data):
        axisscale=self.get_stretched_axis(data,self.stretchdict)        
        na=self.ksmoothing_single_usingscale(data, axisscale)
        return na

    def get_periodic_points(self,scale,pos):
        distm=np.abs(scale-scale[pos])
        #gaplen=(scale[i][-1]-scale[i][-2]+scale[i][1]-scale[i][0])/2
        distl=np.abs(scale[pos]+scale[-1]-scale)
        distr=np.abs(scale[-1]+scale-scale[pos])
        adist=np.vstack([distl,distm,distr])
        return adist.min(axis=0)

    def get_skernel_mask_withscale_singlepos(self,ksize,scale=[],pos=[]):
        ka=np.array(1.0)
        ind1=[]
        ind2=[]
        for i in range(len(ksize)):
            ds=ksize[i]
            if pos[i]==0:
                width=scale[i][0]
            else:
                width=scale[i][pos[i]]-scale[i][pos[i]-1]
            if self.periodic[i]:
                points=self.get_periodic_points(np.arange(1,ds+1),pos[i])
            else:
                points=np.arange(ds)-pos[i]
            nd=scipy.stats.norm.pdf(points,loc=0,scale=1.0/width)
            ka=ka.reshape(ind1+[1])*(nd.reshape(ind2+[len(nd)]))
            ind1.append(len(nd))
            ind2.append(1)
        return ka

      
    def gpsmoothing_single(self,data):
        #pdb.set_trace()
        self.raw_data=data
        self.raw_data[data==0]=2
        rawdatamean=self.raw_data.mean()
        self.data_mean=np.squeeze(self.raw_data)*200/rawdatamean
        saxis=self.get_stretched_axis(data,self.stretchdict)
        if 'mm' in self.sm:
            if self.apo==None:
                self.apo=atmsprop(self.lib)
            self.prior_mean=self.apo.get_repulsion(self.ind,self.bincenters[0])#self.ksmoothing_single_uniform(self.data_mean,self.priorkw)            
        else:
            self.prior_mean=np.zeros(self.data_mean.shape)#self.ksmoothing_single_uniform(self.data_mean,self.priorkw)
        #pdb.set_trace()
        self.calc_prior_covariance(saxis)
        #self.prior_covariance=self.prior_covariance
        self.data_mean.resize(self.data_mean.size)        
        self.prior_mean.resize(self.prior_mean.size)
        smootheddatamean=self.ksmoothing_single_uniform(self.data_mean,self.priorkw)
        #smootheddatamean=np.ones(len(smootheddatamean))*1
        #barrier=np.mean(smootheddatamean)+np.std(smootheddatamean)
        #smootheddatamean[smootheddatamean>barrier]=smootheddatamean[smootheddatamean>barrier]*10
        smootheddatamean=smootheddatamean/10
        if 'a' in self.sm:
            if re.search('([0-9\.]{1,10})', self.sm.split('a')[-1]):
                covariance=float(re.search('([0-9\.]{1,10})', self.sm.split('a')[-1]).group(1))
            else:
                covariance=1./smootheddatamean.mean()
            self.data_covariance=np.eye(len(smootheddatamean))*(covariance)
            self.data_covariance_inv=np.eye(len(smootheddatamean))*(1/covariance)
        else:
            self.data_covariance=np.diag(1./smootheddatamean)
            self.data_covariance_inv=np.diag(smootheddatamean)
        self.calc_posterior_mean()
        ra=np.squeeze(self.posterior_mean)*rawdatamean/200
        ra[ra>100000000000000]=1
        return ra

    def get_periodic_array(self, axc):
        ua=np.abs(axc-axc.T)
        al=np.abs(axc+axc[-1]-axc.T)
        ar=np.abs(axc-axc[-1]-axc.T)
        uas=ua.shape
        adist=np.vstack([ua.reshape([1,ua.size]),ar.reshape([1,ar.size]),al.reshape([1,al.size])])
        return adist.min(axis=0).reshape(uas)
        
    def calc_prior_covariance(self,saxis):
        pp=self.prior_mean.shape
        tn=self.prior_mean.size
        covarray=np.array(1)
        ind1=[]
        ind2=[]
        for i in range(len(pp)):
            axc=saxis[i]
            axc.resize([axc.size,1])
            if self.periodic[i]:
                points=self.get_periodic_array(axc)
            else:
                points=axc-axc.T
            sp=np.exp(-(points)**2)
            covarray=covarray.reshape(ind1+[1]+ind1+[1])*(sp.reshape(ind2+[pp[i]]+ind2+[pp[i]]))
            ind1.append(pp[i])
            ind2.append(1)
        #covarray=np.copy(np.transpose(covarray,range(0,len(pp)*2,2)+range(1,len(pp)*2,2)))
        covarray.resize([tn,tn])
        self.prior_covariance=covarray*self.prior_baseline_variance
        #self.prior_covariance_inv=numpy.linalg.inv(self.prior_covariance)
            
    def calc_data_mean(self):
        self.data_mean=self.prepare_mean(self.raw_data)
        
    def calc_data_covariance(self):
        dcv=self.data_baseline_variance*self.raw_data.astype(np.float32)**self.data_variance_power
        dcv.resize([dcv.size,1])
        #dcv[dcv>1]=1000
        self.data_covariance=np.diag(np.squeeze(dcv))
        self.data_covariance_inv=np.diag(np.squeeze(dcv**-1))
        
    def calc_posterior_mean(self):
#self.posterior_mean=np.dot(self.posterior_covariance,np.dot(self.prior_covariance_inv,self.prior_mean)+np.dot(self.data_covariance_inv,self.data_mean))
        mi=np.eye(len(self.data_mean))
        m2=np.dot(self.prior_covariance,self.data_covariance_inv)
        t1=mi+m2        
        t2=np.dot(m2,self.data_mean)
        self.posterior_mean=numpy.linalg.solve(t1,t2)
        #pdb.set_trace()
#res=np.dot(mi+np.dot(self.prior_covariance,self.data_covariance_inv),self.posterior_mean)-(self.prior_mean+np.dot(self.prior_covariance,np.dot(self.data_covariance_inv,self.data_mean)))
        self.posterior_mean.resize(self.raw_data.shape)
        #pdb.set_trace()

    def calc_posterior_covariance(self):
        self.posterior_covariance=numpy.linalg.inv(self.prior_covariance_inv+self.data_covariance_inv)

    def calc_smoothness(self):
        pass
    
    def runtask_cluster(self,runnum):
        print runnum
        self.runpath=get_path_to_thisfile()
        fh=open('input.pickle')
        nso=pickle.load(fh)
        fh.close()
        nso.smoothing_2d_npa(runnum)
        report_job_runstatus(self.runpath, True, runnum, '.smoothed',inputname='runme.py')
    
    def smoothing_2d_npa(self,npaname):
        orn=npaname
        if not npaname.endswith('npy'):
            npaname=npaname+'.npy'
        ma=np.load(npaname)
        print 'start smoothing'
        self.smoothing_2d(ma,int(orn))
        print 'finish smoothing'
        #save the resulting array
        np.save(orn+'.smoothed',ma)
        
    def smoothing_2d(self,ma,orn):
        ds=ma.shape
        for i in range(ds[0]):
            self.ind=(orn-1)*self.nofrpn+i
            if self.sm.startswith('ks'):
                ma[i]=self.ksmoothing_single_withscale(ma[i])
            elif self.sm.startswith('gps'):
                ma[i]=self.gpsmoothing_single(ma[i])

class atmsprop(object):
    def __init__(self,libfile):
        self.vdrd={'C':1.8,'N':1.64,'O':1.46,'S':1.77,'P':1.88,'F':1.56,'Cl':1.74,'Br':1.98,'I':2.09}
        self.libfile=libfile
        self.atomNameList=self.read_libfile()
        self.atomTypeList=[self.get_atomtype(an) for an in self.atomNameList]
        self.atomrl=[self.vdrd[at] for at in self.atomTypeList]
        
    def get_atomtype(self,atomname):
        if 'Cl' in atomname:
            return 'Cl'
        for key in self.vdrd:
            if key in atomname:
                return key
        raise exception('atom type not known')
    
    def read_libfile(self):
        fh=open(runenv.basedir+'lib/'+self.libfile)
        fhc=fh.read()
        fh.close()
        if 'DBLGRP' in fhc:
            agl=[l[8:-1] for l in fhc.split('\n') if l.startswith('DBLGRP')]
        elif 'ATMGRP' in fhc:
            agl=[l[8:-1] for l in fhc.split('\n') if l.startswith('ATMGRP')]
        return [a for a in agl if len(a)>0]
    
    def get_repulsion(self,ind,bins):
        atm1ind=ind/len(self.atomTypeList)
        atm2ind=ind%len(self.atomTypeList)
        radius=self.atomrl[atm1ind]+self.atomrl[atm2ind]
        pv=np.zeros(len(bins))
        for i in range(len(bins)):
            if bins[i]>=radius:
                break
            pv[i]=np.exp(-0.5*((bins[i]-radius)/0.05)**2)/(0.05*2.51)
        return pv


