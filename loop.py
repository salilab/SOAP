"""
   SOAP loop model.

"""
from env import *

class MyLoop(loopmodel):
    # This routine picks the residues to be refined by loop modeling
    def __init__(self, env, sequence, alnfile=None, knowns=[], inimodel=None,
                 deviation=None, library_schedule=None, csrfile=None,
                 inifile=None, assess_methods=None, loop_assess_methods=None,
                 refinepot=['$(LIB)/atmcls-mf.lib','$(LIB)/dist-mf.lib'],
                 loops=[],calcrmsds='111'):
        loopmodel.__init__(self, env, sequence, alnfile, knowns, inimodel,
                 deviation, library_schedule, csrfile,
                 inifile, assess_methods, loop_assess_methods)
        self.loops=loops
        self.refinepotential=refinepot
        self.load_native_model()
        self.calc_rmsds=calcrmsds
        self.loop.env.schedule_scale = physical.values(default=1.0,
                                                       nonbond_spline=0.02)#0.6
        edat = self.loop.env.edat
        edat.contact_shell=12.00
        edat.dynamic_sphere=True#True
        edat.dynamic_lennard=False#False
        edat.dynamic_coulomb=False#False
        edat.relative_dielectric=1.0
        edat.dynamic_modeller=True#True
        #self.loop.library_schedule


    def read_potential(self):
        return group_restraints(self.env, classes='$(LIB)/atmcls-mf.lib',
                                parameters=runenv.serverUserPath+'sps/soap_loop/dist.lib')


    def loop_restraints(self, atmsel, aln):
        dih_lib_only = True
        mnch_lib = 1
        res_at = 1
        self.restraints.clear()
        for typ in ('bond', 'angle', 'improper', 'dihedral'):
            self.restraints.make(atmsel, aln=aln, restraint_type=typ,
                                 spline_on_site=self.spline_on_site,
                                 dih_lib_only=dih_lib_only,
                                 mnch_lib=mnch_lib, restraint_sel_atoms=res_at)
        for typ in ('omega', 'chi1', 'chi2', 'chi3', 'chi4'):
            self.restraints.make(atmsel, aln=aln,
                                 restraint_type=typ+'_dihedral',
                                 spline_on_site=self.spline_on_site,
                                 dih_lib_only=dih_lib_only, mnch_lib=mnch_lib,
                                 restraint_sel_atoms=res_at, spline_range=4.0,
                                 spline_dx=0.3, spline_min_points=5)
        # Generate restraints on all non-standard residues:
        self.nonstd_restraints(aln)
        self.special_restraints(aln)        
        
    def select_loop_atoms(self):
        s=selection()
        if not isinstance(self.loops[0],list):
            self.add_loop2selection(self.loops,s)
        else:
            for loop in self.loops:
                self.add_loop2selection(loop,s)
        return s
   
    def add_loop2selection(self,loop,s):
        try:
            s.add(selection(self.residue_range(str(loop[2])+':'+loop[0],str(loop[3])+':'+loop[0])))
        except:
            lind=self.residues[str(loop[2])+':'+loop[0]].index
            s.add(selection(self.residues[lind:(lind+loop[1])]))
    
    def loop_model_analysis(self, atmsel, ini_model, filename, out, id1, num):
        """Energy evaluation and assessment, and write out the loop model"""
        self.user_after_single_loop_model(out)
        if self.loop.write_selection_only:
            self.select_loop_atoms().write(file=filename)
        else:
            self.write(file=filename)
        if self.accelrys:
            # Accelrys wants their analysis *after* the model is written, so
            # that the written-out model keeps the original template-derived
            # Biso, rather than getting the energy profile Biso:
            self.loop_model_analysis_accelrys(ini_model, id1, num)

        # Do model assessment if requested
        self.assess(atmsel, self.loop.assess_methods, out)
        
    def user_after_single_loop_model(self,out):
        if not self.calc_rmsds:
            return 0
        ormdl=Mymodel(self.env)
        ormdl.read(file=self.inimodel)
        self.ormdl=ormdl
        self.ors=ormdl.select_loop_atoms(self.loops)
        #self.load_native_model()
        s2=self.select_loop_atoms()
        aln=alignment(self.env)
        aln.append_model(self.ormdl, align_codes='c1', atom_files=self.inimodel)
        aln.append_model(self, align_codes='c2')
        # We can now use the calculated RMS, DRMS, etc. from the returned 'r' object:
        rt=str(self.calc_rmsds)
        if rt[-1]=='1':
            r = self.ors.superpose(self, aln,fit=False,refine_local=False)
            out['rmsd']=min(r.rms,out['rmsd'])
        if len(rt)>=2 and rt[-2]=='1':
            r=self.ors.only_mainchain().superpose(self, aln,fit=False,refine_local=False)
            out['mcrmsd']=min(r.rms,out['rmsd'])
        if len(rt)>=3 and rt[-3]=='1':
            r=self.ors.only_sidechain().superpose(self, aln,fit=False,refine_local=False)
            out['scrmsd']=r.rms
        print r.rms
        return r.rms
    
    def load_native_model(self):
        ormdl=Mymodel(self.env)
        ormdl.read(file=self.inimodel)
        ormdl.write(self.sequence+'.native.pdb')
        self.ormdl=ormdl
        self.ors=ormdl.select_loop_atoms(self.loops)
       
    #def read_potential(self):    
    #    return group_restraints(self.env, classes=self.refinepotential[0],parameters=self.refinepotential[1])   

    def single_loop_model(self, atmsel, ini_model, num, id1, sched,
                          parallel=False):
        """Build a single loop model"""
        self.tracefile = self.new_loop_trace_file(id1, num)

        if parallel:
            self.read_top_par()
            self.read_ini_loop_model(ini_model)
            aln = self.read_loop_alignment(ini_model)
            self.create_loop_topology(aln, ini_model)
        else:
            self.read_ini_loop_model(ini_model)

        atmsel.randomize_xyz(deviation=50)

        if parallel:
            self.group_restraints = self.read_potential()
            self.rd_restraints()

        filename = self.get_loop_model_filename(self.sequence, id1, num,
                                                self.pdb_ext)
        out = {'name':filename, 'loopnum':id1, 'num':num, 'failure':None,'mcrmsd':999.0,'rmsd':999.0}

        self.loop.outputs.append(out)

        actions = self.get_loop_actions()
        try:
            # Refine without the rest of the protein:
            self.env.edat.nonbonded_sel_atoms = 2
            self.optimize_loop(atmsel, sched, actions)
            # Refine in the context of the rest of the protein:
            self.env.edat.nonbonded_sel_atoms = 1
            self.optimize_loop(atmsel, sched, actions)

            (out['molpdf'], out['pdfterms']) = atmsel.energy()

            self.to_iupac()
        except (ModellerError, OverflowError):
            detail = sys.exc_info()[1]
            if len(str(detail)) > 0:
                out['failure'] = detail
            else:
                out['failure'] = 'Optimization failed'
            del detail
        else:
            self.loop_model_analysis(atmsel, ini_model, filename, out, id1, num)
        del self.tracefile
        return out
    
    def get_loop_actions(self):
        """Get actions to carry out during loop optimization"""
        act = []
        act.append(checkRMSD(20,self))
        return act

    def multiple_loop_models(self, atmsel, ini_model, num, sched):
        """Build all loop models for a given initial model"""
        for id1 in range(self.loop.starting_model, self.loop.ending_model + 1):
            r = self.single_loop_model(atmsel, ini_model, num, id1, sched)
            
               
class sprefine(object):
    def __init__(self,dslist='101',bm='loop.100.8',criteria='bestrmsd',spmodel=None):
        #use statistical potential to refine loops/strutures, and analyze the refine results to judge the statistical potential.
        if spmodel: #needfix???
            if spmodel['scoretype']!='sprefine':
                raise Exception('You can not call sprefinescore with this type of model, wrong data type')
            bm=spmodel['bm']
            dslist=model['dslist']
            criteria=model['criteria']
        self.bm=bm
        self.dslist=dslist
        self.criteria=criteria
        #self.codelist=decoysets(self.dslist).get_nativelist()  not applicalbel for loops    
        self.refpot=['$(LIB)/atmcls-mf.lib','$(LIB)/dist-mf.lib']
        bml=bm.split('.')
        self.refinetype=bml[0]
        self.refinemethod=bml[1:-1][0]# [restratins,optimizationmethod,initialmodel,randomizem,numofmodelstobuild]
        self.nofloops=int(bml[-1])
        self.startNum=0
        self.endNum=self.nofloops
        self.slavenumber=int(bm[-1])
        self.codedict={}
        self.pe='local'
        self.calc_rmsds=True
        self.envinitialized=False
        self.j=[]
        self.nors=3000
        self.cv=None
        
    def get_refine_dict(self):
        if self.refinetype=='loop':
            if isinstance(self.dslist,str):
                if not self.dslist.endswith('pickle'):
                    self.dslist=self.dslist+'.pickle'
                self.initialize_dslist()
            elif isinstance(self.dslist,list):
                self.codedict=decoysets(dslist=self.dslist).get_loopdict()
        elif self.refinetype=='native':
            codelist=copy.deepcopy(decoysets(dslist=self.dslist).codelist)
            for code in codelist:
                self.codedict[code]={'toberefined':[code]}
        elif self.refinetype.startswith('decoys'):
            rtr=self.refinetype[6:]
            rtrl=rtr.split('-')
            rmsd=float(rtrl[0])
            nofdecoy=int(rtrl[1])
            ds=decoysets(dslist=self.dslist)
            codelist=copy.deepcopy(ds.codelist)
            for code,i in zip(codelist,range(codelist)):
                self.codedict[code]={}
                ssa=ds.sa[ds.pos[i]:ds.pos[i+1]]
                ri=np.searchsorted(ssa['rmsd'],rmsd)
                if len(ssa)<=nofdecoy:
                    tsp=0
                    tep=len(ssa)
                else:
                    tsp=int(ri-nofdecoy/2)
                    tep=int(ri+nofdecoy/2)
                    if tsp<0:
                        tsp=0
                        tep=nofdecoy
                    elif tep>len(ssa):
                        tep=len(ssa)
                        tsp=tep-nofdecoy
                self.codedict[code]['sa']=np.copy(ssa[tsp:tep])
                self.codedict[code]['toberefined']=copy.deepcopy(ds.dnlist[tsp:tep])
        
    def initialize_dslist(self):
        if runenv.hostn:
            fh=open(self.dslist+'.pickle')
        else:
            fh=open(runenv.basedir+'loop/'+self.dslist+'.pickle')
        self.codedict=pickle.load(fh)
        if 'runenv' in self.codedict:
            self.startNum=self.codedict['runenv'][0]
            self.endNum=self.codedict['runenv'][1]
            del self.codedict['runenv']
        fh.close()
    
    def get_ds_part(self,sample):
        newcodedict={}
        for key in self.codedict:
            if key in sample:
                newcodedict[key]=self.codedict[key]
        nlr=sprefine(dslist=self.dslist,bm=self.bm,criteria=self.criteria)
        nlr.codedict=newcodedict
        return nlr

    def initialize_runenv(self):
        if self.pe=='local':
            self.envinitialized=True
            self.j = job()
            for i in range(0,self.slavenumber):
                self.j.append(local_slave())
        elif self.pe=='cluster':
            self.envinitialized=True
            self.j  = sge_qsub_job(maxslave=self.slavenumber, options="-l arch=lx24-amd64")
            #self.j.append(local_slave())
    
    def refine(self,codes):
        if self.refinetype=='loop':
            self.model_loops_list(self.codedict)
        else:
            pass
    
    def refine_single_code(self):
        pass
      
    def model_loops(self, input):
        log.none()
        if input[-7:]!='.pickle':
            input=input+'.pickle'
        self.input=input[:-7]
        print input
        print input[:-7]
        fh=open(input,'rb')
        loopdict=cPickle.load(fh)
        fh.close()
        self.model_loops_list(loopdict)
        fl=os.listdir('./')
        fl=[item for item in fl if item[-3:]=='pdb']
        for item in fl:
            print os.system('mv '+item+' loop'+str(self.input)+'.'+item)
                        
    def model_loops_list_old(self,loopdict):
        resultdict={}
        for key in loopdict:
            print key
            loops=loopdict[key]
            for loop in loops:
                print loop
                loopname=key+str(loop[1])+loop[0]+str(loop[2])
                print os.system('rm -r '+'loop.'+loopname)
                #print os.mkdir('loop.'+loopname)
                #print os.chdir('loop.'+loopname)                  
                try:                  
                    res=self.loop_model_single_protein(modelfile='pdb'+key+'.ent',loops=[loop], modelname=loopname)
                    resultdict[loopname]=res
                    #print os.system('mv '+loopname+'.output.pickle loop'+str(self.input)+'.'+loopname+'.output.pickle')
                except Exception,e:
                    print >> sys.stderr, e
                print os.system('rm '+loopname+'.lrsr')
                print os.system('rm '+loopname+'.DL*')                    
                #print os.chdir('..')
        return resultdict
            
    def loop_model_single_protein(self,modelfile,loops,modelname):
        numofmodels=self.nofloops
        print "number of loops to be calculated: "+str(numofmodels)
        try:
            rs=-2-int(self.dslist)
        except:
            rs=-2
        env = environ(rand_seed=rs)
        # directories for input atom files
        env.io.atom_files_directory = [os.path.join(runenv.basedir,'loop/',self.dslist),scratchdir,runenv.opdbdir,'.']
        # Create a new class based on 'loopmodel' so that we can redefine
        # select_loop_atoms (necessary)
        m = MyLoop(env,
                   inimodel=modelfile, # initial model of the target
                   sequence=modelname,loop_assess_methods=(assess.DOPE),
                   refinepot=self.refpot,loops=loops)
        #m.load_native_model()
        m.loop.starting_model= self.startNum           # index of the first loop model
        m.loop.ending_model  = self.endNum
        print "refine"
        print self.refinemethod
        if self.refinemethod=='fast':
            m.loop.md_level = refine.fast
        elif self.refinemethod=='none':
            m.loop.md_level= None
        elif self.refinemethod=='super':
            m.loop.md_level= refine.very_fast
        elif self.refinemethod=='slow':
            m.loop.md_level= refine.slow
        elif self.refinemethod=='sgmd':
            print "sgmd"
            m.loop.md_lvel=sgmd
        m.loop.write_selection_only=True
        if len(self.j)>0:
            m.use_parallel_job(self.j)   
        m.make()
        resultlist=[]
        nativescore=m.ors.assess_dope()
        loopoutput=m.loop.outputs
        #pdb.set_trace()
        for i in range(0,len(loopoutput)):
            try:
                resultlist.append((loopoutput[i]['rmsd'],loopoutput[i]['mcrmsd'],loopoutput[i]['DOPE score']))
            except Exception,e:
                print >> sys.stderr, e 
                continue
        #rmsddict[modelname+'.native.pdb']=0.0
        #dopedict[modelname+'.native.pdb']=nativescore
        #fh=open('rmsd.pickle','wb')
        #pickle.dump(rmsddict,fh)
        #fh.close()
        #fh=open('dope.pickle','wb')
        #pickle.dump(dopedict,fh)
        #fh.close()
        #fh=open(modelname+'.output.pickle','wb')
        #pickle.dump(m.loop.outputs,fh)
        #fh.close()
        return resultlist

    def prepare_task_input(self):
        self.initialize_dslist()
        self.dir=runenv.runs
        self.task.get_tasksn()
        self.rundir=str(self.task.rsn)
        self.runpath=runenv.serverUserPath+self.rundir+'/'
        #self.task=task('','',afterprocessing=self,preparing=self)        
        self.task.dir=self.dir+self.rundir+'/'
        self.task.rdirname=self.rundir
        #return self.task
        os.chdir(self.dir)
        os.system('rm -r '+self.rundir)
        os.mkdir(self.rundir)
        os.chdir(self.dir+self.rundir)
        freememory=1
        self.freememory=freememory
        loopdict=self.codedict
        #separate the dslist into individual small sets
        nors=self.nors
        dictlist=[]
        trl=0
        looprtdict={}
        for key in loopdict:
            clrl=0
            for loop in loopdict[key]:
                clrl=clrl+loop[1]**2
            looprtdict[key]=clrl
            trl=trl+clrl
        art=trl/float(nors)
        k=0
        cl=0
        ndict={}
        for key in loopdict:
            ndict[key]=loopdict[key]
            cl=cl+looprtdict[key]
            #if cl>art:
            dictlist.append(ndict)
            ndict={}
            cl=0
            k=k+1
        if len(ndict)>0:
            dictlist.append(ndict)    
        print k
        tnt=k*self.nofloops
        numofruns=min(3000,tnt/20)
        nrpl=max(1,numofruns/k)
        nrpn=self.nofloops/nrpl
        for i in range(0,k*nrpl):
            fh=open(str(i+1)+'.pickle','wb')
            cd=dictlist[i/nrpl]
            if nrpl>1:
                cd['runenv']=[(i%k)*nrpn,(i%k+1)*nrpn]
            pickle.dump(cd,fh)
            fh.close()
        self.nrpl=nrpl
        self.nors=k*nrpl
        nors=k*nrpl
        #write the input list for different runs
        inputlist=open('inputlist','w')
        inputlist.write(','.join([str(i)+'.pickle' for i in range(1,nors+1)]))
        inputlist.close()
        makemdt=open('runme.py','w')
        makemdt.write('from SOAP.loop import *\nimport sys \n \nspopt=sprefine(sys.argv[1],\''+self.bm+'\',\''+self.criteria+'\')\nspopt.runpath=\''+self.runpath+'\'\n\nspopt.initialize_dslist()'+'\nspopt.initialize_runenv()'+'\nspopt.assess_cluster_node()\n')
        makemdt.flush()
        dp=self.slavenumber
        generate_job_submit_script(freememory,self.rundir,runtime,nors,parallel=self.slavenumber)        
        return 0

    def runtask_cluster(self,nors):
        self.nors=nors
        self.task=task('','',afterprocessing=self,preparing=self)
        self.prepare_task_input()
        return self.task.run_task_cluster()    

    def get_task(self,cv=None):
        if cv!=None:
            self.cv=cv
            self.refpot[1]=os.path.join(runenv.serverUserPath+'lib',cv.rundir+'.lib')
        self.task=task('','',afterprocessing=self,preparing=self)
        return self.task
    
    def model_loops_list(self,loopdict):
        resultdict={}
        for key in loopdict:
            print key
            loops=loopdict[key]
            for i,loop in enumerate(loops):
                print loop
                loopname=key+str(i)
                print os.system('rm -r '+'loop.'+loopname)                 
                #try:
                if 1:
                    res=self.loop_model_single_protein(modelfile='pdb'+key+'.ent',loops=loop, modelname=loopname)
                    resultdict[loopname]=res
                #except Exception,e:
                #    print >> sys.stderr, e
                    #pdb.set_trace()
                print os.system('rm '+loopname+'.lrsr')
                print os.system('rm '+loopname+'.DL*')                    
        return resultdict          

    def assess_cluster_node(self):
        os.chdir(scratchdir)
        result=self.model_loops_list(self.codedict)
        #if self.criteria=='bestrmsd':
        #    result=self.assess_bestrmsd(result)
        pickle.dump(result,open(self.dslist+'.pickle','w'))
        report_job_runstatus(self.runpath, True, self.dslist, '.pickle',inputname='runme.py',temppath=scratchdir)
        
    def afterprocessing(self):
        os.chdir(self.dir+self.rundir)
        rd={}
        print self.nrpl
        print self.nors
        for i in range(self.nors):
            if os.path.isfile(str(i+1)+'.pickle'):
                fh=open(str(i+1)+'.pickle','rb')
                res=pickle.load(fh)
                fh.close()
                if res.keys()[0] in rd:
                    rd[res.keys()[0]].extend(res.values()[0])
                else:
                    rd[res.keys()[0]]=res.values()[0]
        #print os.system('rm -r '+self.dir+self.rundir)
        #del self.clusters
        self.result=rd
        #pdb.set_trace()
        
        del self.task
        gc.collect()
        return self.analyze_loop_modeling()
    
    def assess(self,refpot=[]):
        if refpot:
            self.refpot=refpot
        if runenv.hostn==0:
            self.runtask_cluster(500)
            result=self.result
        else:
            result=self.model_loops_list(self.codedict)
        if self.criteria=='bestrmsd':
            return self.analyze_loop_modeling_results(result)
            #return self.assess_bestrmsd(result)
            
    def assess_local(self,refpot=[]):
        if refpot:
            self.refpot=refpot
        result=self.model_loops_list(self.codedict)
        if self.criteria=='bestrmsd':
            return self.analyze_loop_modeling_results(result)        

    def analyze_loop_modeling(self):
        mn=9999999999
        result=self.result
        for key in result:
            if len(result[key])<mn:
                mn=len(result[key])
        print "min n "+str(mn)
        #mn=1000
        ra=np.zeros([len(result.keys()),mn,3])
        k=0
        for key in result:
            ra[k,:,:]=np.array(result[key])[:mn,:]
            k+=1
        self.averageRMSD=ra[:,:,0].min(axis=1).mean()+ra[:,:,1].min(axis=1).mean()
        if self.cv!=None:
            self.cv.resultsummary=self.averageRMSD
        return self.averageRMSD

    def analyze_loop_modeling_results(self,result):
        mn=9999999999
        for key in result:
            if len(result[key])<mn:
                mn=len(result[key])
        print "min n "+str(mn)
        #mn=1000
        ra=np.zeros([len(result.keys()),mn,3])
        k=0
        for key in result:
            ra[k,:,:]=np.array(result[key])[:mn,:]
            k+=1
        interval=max(1,mn/1000)
        minrmsd=[]
        minmcrmsd=[]
        minscored=[]
        for i in range(interval,mn,interval):
            print i
            minrmsd.append(ra[:,:i,0].min(axis=1).mean())
            minmcrmsd.append(ra[:,:i,1].min(axis=1).mean())
            minind=ra[:,:i,2].argmin(axis=1)
            minscored.append(np.median(ra[np.arange(ra.shape[0]),minind,np.zeros(ra.shape[0],dtype=np.int)]))
            print minrmsd[-1]
            print minmcrmsd[-1]
        x=range(interval,mn,interval)
        phs=plt.plot(x,minrmsd)
        phs2=plt.plot(x,minmcrmsd)
        #phs3=plt.plot(x,minscored)
        plt.xlabel('Number of models')
        plt.ylabel('Average RMSD')
        plt.legend((phs[0],phs2[0]),('All atom','Main chain'),loc=0)        
        plt.savefig(self.rundir+self.refinemethod+'.eps')
        #return minrmsd[-1]
        
    def assess_bestrmsd(self,result):
        brmsdlist=[]
        for code in result:
            rmsddict=result[code][0]
            rmsdlist=[]
            for key in rmsddict:
                if rmsddict[key]>0:
                    rmsdlist.append(rmsddict[key])
            brmsdlist.append(np.array(rmsdlist).min())
        return brmsdlist
