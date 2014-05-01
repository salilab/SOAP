"""
   SOAP RefineScore module 

"""

from env import *

    
class sprefinescore(object):
    def __init__(self,dslist='',bm='fast',criteria='bestrmsd',pn=8, model=[]):
        if spmodel: #needfix???
            if spmodel['scoretype']!='sprefine':
                raise Exception('You can not call sprefinescore with this type of model, wrong data type')
            bm=spmodel['bm']
            dslist=model['dslist']
            criteria=model['criteria']
        self.bm=bm
        self.dslist=dslist
        self.criteria=criteria
        self.codelist=decoysets(self.dslist).get_nativelist()
        
    def get_refprotein_list(self):
        pass
    
    def get_ds_part(self):
        pass
    
    def prepare_task_input(self):
        os.chdir(self.dir)
        os.system('rm -r '+self.scorename)
        os.mkdir(self.scorename)
        os.chdir(self.dir+self.scorename)
        nors=len(self.codelist)
        #write the input list for different runs
        inputlist=open('inputlist','w')
        inputlist.write(','.join(['runme.py' for i in range(1,nors+1)]))
        inputlist.close()
        #copy scripts to this rundir
        os.system('cp '+runenv.scriptpath+' ./')
        #generate the lib file for modeller refinement
        self.ssp.write_lib_dist150()
        os.system('cp '+self.ssppath+'.lib ./')
        #write the run python script file
        makemdt=open('runme.py','w')
        makemdt.write('from SOAP import *\nimport sys \n \nrsp=rawsp(e,\''+self.pdbset+'\',\''
                      +self.features+'\',\''+self.genmethod+'\')\n'+'ssp=scaledsp(e,\''+self.pm+'\',rsp)\n'
                      +'so=sprefinescore(e,\''+self.refinem+'\',ssp,dslist='+str(self.dslist)+')\n'
                      +'so.refine(sys.argv[1])')
        makemdt.flush()
        #write the job submission script
        runmdt=open('runme.sh','w')
        submdt=file(runenv.scriptdir+'md_refine.sh','r')
        submdtf=submdt.read()
        submdtf=submdtf.replace('undefinedrundir',self.scorename)
        submdtf=submdtf.replace('tundefined',str(1)+'-'+str(nors))
        submdtf=submdtf.replace('undefinedlist','('+' '.join(self.codelist)+')')
        runmdt.write(submdtf)
        runmdt.flush()
    
    def process_task_output(self):
        os.chdir(self.dir+self.scorename)
        self.score=np.zeros(len(self.codelist),dtype=[('codeposn','i2'),('nativetotalscore','f4'),('nativescscore','f4')
                                              ,('nativespscore','f4'),('decoytotalscore','f4'),('decoyscscore','f4')
                                              ,('decoyspscore','f4'),('rmsd','f4')])
        print self.score
        for i in range(0,len(self.codelist)):
            code=self.codelist[i]
            if os.path.isfile(self.scorepath+'/'+code+'.refine.score'):
                fh=open(self.scorepath+'/'+code+'.refine.score')
                fc=fh.read()
                self.score['codeposn'][i]=i
                rer=re.search(code+' score: ([0-9.-]*),([0-9.-]*),([0-9.-]*)',fc)
                self.score['nativetotalscore'][i]=float(rer.group(1))
                self.score['nativescscore'][i]=float(rer.group(2))
                self.score['nativespscore'][i]=float(rer.group(3))
                rer=re.search(code+'.refine.pdb score: ([0-9.-]*),([0-9.-]*),([0-9.-]*)',fc)
                self.score['decoytotalscore'][i]=float(rer.group(1))
                self.score['decoyscscore'][i]=float(rer.group(2))
                self.score['decoyspscore'][i]=float(rer.group(3))
                rer=re.search('armsd ([0-9.]*)',fc)
                self.score['rmsd'][i]=float(rer.group(1))
        self.scorestats=np.zeros(1,dtype=[('codeposn','i2'),('nativetotalscore','f4'),('nativescscore','f4')
                                              ,('nativespscore','f4'),('decoytotalscore','f4'),('decoyscscore','f4')
                                              ,('decoyspscore','f4'),('rmsd','f4')])
        for item in ['nativetotalscore','nativescscore','nativespscore','decoytotalscore','decoyscscore','decoyspscore','rmsd']:
            self.scorestats[item]=self.score[item].mean()
        self.scorestats['codeposn']=len(self.codelist)

    def start_get_score(self):
        self.prepare_task_input()
        sj=task(self.scorepath+'/',self.scorename)
        sj.run_task_cluster()
     
    def calc_score_cluster(self):
        self.prepare_task_input()
        sj=task(self.scorepath+'/',self.scorename)
        sj.run_task_cluster()
        self.process_task_output()
        self.save()
        
    def calc_score_local(self):
        self.prepare_task_input()
        sj=task(self.scorepath+'/',self.scorename)
        sj.run_task_local(self.refine,self.codelist)

    def refine(self,code):
        rer=re.search('md([0-9]+)',self.refinem)
        genrn=int(rer.group(1))
        if genrn==0:
            self.mdrefine(code)
        if genrn==1:
            temps=[300.0]
            samplits=20000
            self.mdrefine(code,templ=temps,samplits=samplits)
        if genrn==2:
            print 0

    def mdrefine(self,code,templ=[],samplits=0, step=0):
	log.verbose()	
	env = environ()
	env.io.atom_files_directory = [runenv.ddbdir,'./']
	env.libs.topology.read(file='$(LIB)/top_heav.lib')
	env.libs.parameters.read(file='$(LIB)/par.lib')
	env.edat.contact_shell = 15.0 
	env.edat.excl_local=(True, True, True, True)
	env.edat.dynamic_coulomb =  False 
	env.edat.dynamic_lennard = False 
	env.edat.dynamic_sphere = False   
	env.edat.dynamic_modeller = True	
	aln=alignment(env)
        pir( code).make_pir_single_forautomodel()
	aln.append(file= code+'.pir')
	mdl=automodel.automodel(env, alnfile=code+'.pir', knowns=code, sequence=code+'.refine.ini')
        mdl.clear_topology()
	mdl.generate_topology(aln[code+'decoys'],patch_default=True)
	print "generating topology finished"
	mdl.patch_ss_templates(aln)
        mdl.transfer_xyz(aln, cluster_cut=-1.0)
	mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
        mdl.write(code+'.refine.native.pdb')
	atmsel = selection(mdl)
	rsr = mdl.restraints
	rsr.clear()
	#charm restraint
	rsr.make(atmsel, restraint_type='stereo',spline_on_site=True)
	gprsr = group_restraints(env, classes='$(LIB)/atmcls-mf.lib', parameters=self.sspname+'.lib')
	mdl.group_restraints = gprsr
	# Generate restraints on non standard residues:
	mdl.nonstd_restraints(aln)	
	# save restraints
	mdl.restraints.condense()
	scoref=open(code+'.refine.score','w')
	(totalscore,terms) = atmsel.energy()
        stereoscore=terms[physical.bond]+terms[physical.angle]+terms[physical.dihedral]+terms[physical.improper]
	scoref.write(code+' score: '+str(totalscore)+','+str(stereoscore)+','+str(terms[physical.nonbond_spline])+'\n')
        timestep=10.0           
        equil_its=200
        equil_equil=20
        equil_temps=(150.0, 250.0, 400.0, 700.0, 1000.0,1500.0)
        sampl_its=2000
        sampl_equil=200
        sampl_temps=(1500.0,1000.0, 800.0, 600.0, 500.0, 430.0, 370.0, 320.0, 300.0)
        if timestep:
            teimstep=step
        if templ:
            equil_temps=templ.reverse()
            sampl_temps=templ
        if samplits:
            sampl_its=samplits
        md = optimizers.molecular_dynamics(cap_atom_shift=0.39, md_time_step=2.0,md_return='FINAL')
        init_vel = True
        # First run for equilibration, the second for sampling:
        print 'starting refinement'
        k=0
        for (its, equil, temps) in ((equil_its, equil_equil, equil_temps),(sampl_its, sampl_equil, sampl_temps)):
            for temp in temps:
                print temp
                md.optimize(atmsel, max_iterations=its, equilibrate=equil,temperature=temp, init_velocities=init_vel,
                            actions=[optimizers.actions.write_structure(100, 'md'+str(k)+'.%05d.pdb',last=True)])
                init_vel=False
                k=k+1
        filelist=os.listdir('./')
        mdfilelist=[]
        for item in filelist:
            if item[:2]=='md' and item[-3:]=='pdb':
                mdfilelist.append(item)
        self.cluster_models(mdfilelist,mdl,env,cluster_cut=1.5)
        mdl.write(code+'.refine.pdb')
        atmstl=selection(mdl)
        (totalscore,terms) = atmsel.energy()
        stereoscore=terms[physical.bond]+terms[physical.angle]+terms[physical.dihedral]+terms[physical.improper]
	scoref.write(code+'.refine.pdb score: '+str(totalscore)+','+str(stereoscore)+','+str(terms[physical.nonbond_spline])+'\n')
        mdl2=model(env)
        mdl2.read(code)
        r = atmsel.superpose(mdl2, aln)
        scoref.write('armsd '+str(r.rms))
        scoref.close()

    def cluster_models(self,mdfilelist,mdl, env,cluster_cut=1.5):
        aln=alignment(env)
        m=model(env)
        for mf in mdfilelist:
            m.read(mf)
            aln.append_model(mdl=m,align_codes=mf[:-4],atom_files=mf[:-4])  
        aln.malign3d(gap_penalties_3d=(0, 3), fit=False)
        aln.append_model(mdl, align_codes='cluster',
                         atom_files='cluster.opt')
        mdl.transfer_xyz(aln,cluster_cut=cluster_cut)
        mdl.write(file='cluster.ini')
        atmsel=selection(mdl)
        cg = optimizers.conjugate_gradients()            
        cg.optimize(atmsel, max_iterations=5000, output=mdl.optimize_output, min_atom_shift=0.01)

class refinemodel(model):
    def decode_refinem(self,refinem):
        rsrm=refinem[0]
        optm=refinem[1]
        self.initialm=refinem[2]
        self.randomm=refinem[3]
        self.starting_model=1
        self.ending_model=int(refinem[4])

    def get_cgschedule(self):
        pass
    
    def get_mdschedule(self):
        pass
    
    def select_atoms(self):
        """Select atoms to be optimized in the model building procedure. By
           default, this selects all atoms, but you can redefine this routine
           to select a subset instead."""
        return selection(self)
    
    def homcsr(self, exit_stage):
        """Construct the initial model and restraints"""
        # Check the alignment
        aln = self.read_alignment()

        self.check_alignment(aln)

        # make topology and build/read the atom coordinates:
        self.make_initial_model(aln)

        # exit early?
        if exit_stage == 2:
            return

        # make and write the stereochemical, homology, and special restraints?
        if self.create_restraints:
            self.mkhomcsr(selection(self), aln)
            self.restraints.condense()
            self.restraints.write(self.csrfile)

    def mkhomcsr(self, atmsel, aln):
        """Construct typical comparative modeling restraints"""
        rsr = self.restraints
        rsr.clear()
        if self.library_restraints is not None:
            self.build_library_restraints(atmsel, rsr, self.library_restraints)
        else:
            self.build_charmm_restraints(atmsel, rsr, aln)

        # Generate homology-derived distance restraints:
        self.distance_restraints(atmsel, aln)

        # Generate restraints on non standard residues:
        self.nonstd_restraints(aln)

        # Special restraints have to be called last so that possible cis-proline
        # changes are reflected in the current restraints:
        self.special_restraints(aln)
        
    def make_initial_model(self, aln):
        """Make initial model topology and build/read the atom coordinates"""
        self.generate_method(self, aln)
        self.write(file=self.inifile)
        self._check_model_hetatm_water()    

    def randomize_initial_structure(self, atmsel):
        """Get and randomize the initial structure"""
        self.read_initial_model()
        if self.rand_method:
            self.rand_method(atmsel)
    
    def single_model(self, atmsel, num, parallel=False):
        """Build a single optimized model from the initial model"""
        self.tracefile = self.new_trace_file(num)

        if parallel:
            self.read_top_par()
            self.create_topology(self.read_alignment())

        self.randomize_initial_structure(atmsel)

        if parallel:
            self.rd_restraints()

        if not hasattr(self.library_schedule, "make_for_model"):
            raise TypeError("""
library_schedule should now be a schedule object, not an integer as in
older versions of Modeller""")
        sched = self.library_schedule.make_for_model(self)
        sched = sched * runenv.schedule_scale
        fh = open(self.schfile, "w")
        sched.write(fh)
        fh.close()
                
    def single_model_pass(self, atmsel, num, sched):
        """Perform a single pass of model optimization"""
        actions = self.get_optimize_actions()
        for (numstep, step) in enumerate(sched): #slow, normal, fast, very_fast, fastest
            molpdf = step.optimize(atmsel, output=self.optimize_output,
                                   max_iterations=self.max_var_iterations,
                                   actions=actions)
            self.write_int(numstep + 1, num)
            # Also check for molpdf being NaN (depends on Python version; on 2.3
            # x86 it evaluates as equal to everything; with 2.4 x86 it is
            # not greater or smaller than anything)
            if molpdf > self.max_molpdf \
               or (molpdf == 0. and molpdf == 1.) \
               or (not molpdf >= 0. and not molpdf < 0):
                log.error('single_model',
                          "Obj. func. (%.3f) exceeded max_molpdf (%.3f) " \
                                      % (molpdf, self.max_molpdf))
        actions = self.get_refine_actions()
        self.refine(atmsel, actions)

    def refine(self, atmsel, actions): #very_fast, fast,slow, very_slow, slow_large
        """Refine the optimized model with MD and CG"""
        # Save the current model:
        if self.fit_in_refine != 'NO_FIT':
            self.write(file='TO_BE_REFINED.TMP')

        # Possibly skip selecting hot atoms only and optimize all atoms:
        if self.refine_hot_only:
            self.initial_refine_hot(atmsel)

        # Do simulated annealing MD:
        if self.md_level:
            self.md_level(atmsel, actions)

        # Possibly skip 'HOT CG' after MD:
        if self.refine_hot_only:
            self.final_refine_hot(atmsel)

        # Get a final conjugate gradients refined structure:
        cg = conjugate_gradients()
        cg.optimize(atmsel, max_iterations=200, output=self.optimize_output,
                    actions=actions)

        # Evaluate gross changes between the initial and final refined model:
        if 'NO_FIT' not in self.fit_in_refine:
            aln = alignment(self.env)
            mdl2 = read_model(file='TO_BE_REFINED.TMP')
            casel = selection(self).only_atom_types('CA')
            casel.superpose(mdl2, aln)
            casel = selection(self)
            casel.superpose(mdl2, aln)
            modfile.delete('TO_BE_REFINED.TMP')
                
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
        self.runpath='/netapp/sali/gqdong/'+self.rundir+'/'
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
        os.system('cp '+runenv.scriptpath+' ./')
        makemdt=open('runme.py','w')
        makemdt.write('from SOAP import *\nimport sys \n \nspopt=sprefine(sys.argv[1],\''+self.bm+'\',\''+self.criteria+'\')\nspopt.runpath=\''+self.runpath+'\'\n\nspopt.initialize_dslist()'+'\nspopt.initialize_runenv()'+'\nspopt.assess_cluster_node()\n')
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
            self.refpot[1]=os.path.join('/netapp/sali/gqdong/lib',cv.rundir+'.lib')
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
        if hcn=='bell.ucsf.edu':
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
