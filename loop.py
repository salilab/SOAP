"""
   SOAP loop model.

"""
from env import *
from modeller import physical
from modeller.schedule import schedule, step
from modeller.optimizers import conjugate_gradients as CG
from modeller.automodel.autosched import *
log.none()

class MyLoop(loopmodel):
    """
    SOAP loop modeling class
    """
    # This routine picks the residues to be refined by loop modeling
    def __init__(self, env, sequence, alnfile=None, knowns=[], inimodel=None,
                 deviation=None, library_schedule=None, csrfile=None,
                 inifile=None, assess_methods=None, loop_assess_methods=None,
                 refinepot=['$(LIB)/atmcls-mf.lib','$(LIB)/dist-mf.lib'],
                 loops=[],calcrmsds='111',nonbond_spine=0.1,contact_shell=12.0,
                 deviations=50,energytrace=False,assess_trace=True):
        loopmodel.__init__(self, env, sequence, alnfile, knowns, inimodel,
                 deviation, library_schedule, csrfile,
                 inifile, assess_methods, loop_assess_methods)
        self.loops=loops
        self.refinepotential=refinepot
        #self.load_native_model()
        self.calc_rmsds=calcrmsds
        self.deviations=deviations
        self.energytrace=energytrace
        self.assess_trace=assess_trace
        self.loop.env.schedule_scale = physical.values(default=1.0,
                                                       nonbond_spline=nonbond_spine)#0.6
        edat = self.loop.env.edat
        edat.contact_shell=contact_shell
        edat.dynamic_sphere=True#True
        edat.dynamic_lennard=False#False
        edat.dynamic_coulomb=False#False
        edat.relative_dielectric=1.0
        edat.dynamic_modeller=True#True
        #self.loop.library_schedule
        self.loop.library_schedule=loopschedule()
        self.rmsd_calc_initialized=False

    def initilialize_rmsd_calcualtion(self):
        ormdl=Mymodel(self.env)
        ormdl.read(file=self.inimodel)
        self.ormdl=ormdl
        self.ors=ormdl.select_loop_atoms(self.loops)
        #self.load_native_model()
        s2=self.select_loop_atoms()
        self.s2=s2
        aln=alignment(self.env)
        aln.append_model(self.ormdl, align_codes='c1', atom_files=self.inimodel)
        aln.append_model(self, align_codes='c2')
        self.aln=aln
        self.rmsd_calc_initialized=True

    def read_potential(self):
        return group_restraints(self.env, classes=self.refinepotential[0],
                                parameters=self.refinepotential[1])


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
        if loop[0]=='':
            loop[0]=self.chains[0].name
        try:
            s.add(selection(self.residue_range(str(loop[2])+':'+loop[0],str(loop[3])+':'+loop[0])))
        except:
            lind=self.residues[str(loop[2])+':'+loop[0]].index
            s.add(selection(self.residues[lind:(lind+loop[1])]))
    
    def loop_model_analysis(self, atmsel, ini_model, filename, out, id1, num):
        """Energy evaluation and assessment, and write out the loop model"""
        self.user_after_single_loop_model(out,last=True)
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
        
    def user_after_single_loop_model(self,out,last=False):
        if not self.calc_rmsds:
            return 0
        if not self.rmsd_calc_initialized:
            self.initilialize_rmsd_calcualtion()        
        # We can now use the calculated RMS, DRMS, etc. from the returned 'r' object:
        rt=str(self.calc_rmsds)
        if rt[-1]=='1':
            r = self.ors.superpose(self, self.aln,fit=False,refine_local=False)
            out['rmsd']=min(r.rms,out['rmsd']) if 'rmsd' in out  else r.rms
        if len(rt)>=2 and rt[-2]=='1':
            r=self.ors.only_mainchain().superpose(self, self.aln,fit=False,refine_local=False)
            mcrmsd=r.rms
            out['mcrmsd']=min(r.rms,out['mcrmsd']) if 'mcrmsd' in out  else r.rms
        if len(rt)>=3 and rt[-3]=='1':
            r=self.ors.only_sidechain().superpose(self, self.aln,fit=False,refine_local=False)
            out['scrmsd']=r.rms
        if self.energytrace and last:
            if not 'trace' in out:
                out['trace']=[]
            dopescore=self.loop.assess_methods[0](self.s2)[1] if last else 9999999
            soapscore=self.loop.assess_methods[1](self.s2)[1] if last else 9999999
            dopescore=999999999 if np.isnan(dopescore) else dopescore
            soapscore=999999999 if np.isnan(soapscore) else soapscore
            out['trace'].append((mcrmsd,dopescore,soapscore))                
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

        atmsel.randomize_xyz(deviation=self.deviations)

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
        if self.assess_trace:
            act.append(checkRMSD(20,self))
        return act

    def multiple_loop_models(self, atmsel, ini_model, num, sched):
        """Build all loop models for a given initial model"""
        for id1 in range(self.loop.starting_model, self.loop.ending_model + 1):
            r = self.single_loop_model(atmsel, ini_model, num, id1, sched)
            
class ORLoop(loopmodel):
    """
    Original loop modeling class
    """
    # This routine picks the residues to be refined by loop modeling
    # This routine picks the residues to be refined by loop modeling
    def __init__(self, env, sequence, alnfile=None, knowns=[], inimodel=None,
                 deviation=None, library_schedule=None, csrfile=None,
                 inifile=None, assess_methods=None, loop_assess_methods=None,
                 refinepot=['$(LIB)/atmcls-mf.lib','$(LIB)/dist-mf.lib'],assess_trace=False,
                 loops=[],calcrmsds='111',nonbond_spine=0.1,contact_shell=12.0,deviations=50,energytrace=False):
        loopmodel.__init__(self, env, sequence, alnfile, knowns, inimodel,
                 deviation, library_schedule, csrfile,
                 inifile, assess_methods, loop_assess_methods)
        self.loops=loops
        #self.load_native_model()
        self.calc_rmsds=calcrmsds
        self.energytrace=energytrace
        self.assess_trace=assess_trace
        self.rmsd_calc_initialized=False

    def initilialize_rmsd_calcualtion(self):
        ormdl=Mymodel(self.env)
        ormdl.read(file=self.inimodel)
        self.ormdl=ormdl
        self.ors=ormdl.select_loop_atoms(self.loops)
        #self.load_native_model()
        s2=self.select_loop_atoms()
        self.s2=s2
        aln=alignment(self.env)
        aln.append_model(self.ormdl, align_codes='c1', atom_files=self.inimodel)
        aln.append_model(self, align_codes='c2')
        self.aln=aln
        self.rmsd_calc_initialized=True
       
    def select_loop_atoms(self):
        s=selection()
        if not isinstance(self.loops[0],list):
            self.add_loop2selection(self.loops,s)
        else:
            for loop in self.loops:
                self.add_loop2selection(loop,s)
        return s
   
    def add_loop2selection(self,loop,s):
        if loop[0]=='':
            loop[0]=self.chains[0].name
        try:
            s.add(selection(self.residue_range(str(loop[2])+':'+loop[0],str(loop[3])+':'+loop[0])))
        except:
            lind=self.residues[str(loop[2])+':'+loop[0]].index
            s.add(selection(self.residues[lind:(lind+loop[1])]))

    def loop_model_analysis(self, atmsel, ini_model, filename, out, id1, num):
        """Energy evaluation and assessment, and write out the loop model"""
        self.user_after_single_loop_model(out,last=True)
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
        
    def user_after_single_loop_model(self,out,last=False):
        if not self.calc_rmsds:
            return 0
        if not self.rmsd_calc_initialized:
            self.initilialize_rmsd_calcualtion()        
        # We can now use the calculated RMS, DRMS, etc. from the returned 'r' object:
        rt=str(self.calc_rmsds)
        if rt[-1]=='1':
            r = self.ors.superpose(self, self.aln,fit=False,refine_local=False)
            out['rmsd']=min(r.rms,out['rmsd']) if 'rmsd' in out  else r.rms
        if len(rt)>=2 and rt[-2]=='1':
            r=self.ors.only_mainchain().superpose(self, self.aln,fit=False,refine_local=False)
            mcrmsd=r.rms
            out['mcrmsd']=min(r.rms,out['mcrmsd']) if 'mcrmsd' in out  else r.rms
        if len(rt)>=3 and rt[-3]=='1':
            r=self.ors.only_sidechain().superpose(self, self.aln,fit=False,refine_local=False)
            out['scrmsd']=r.rms
        if self.energytrace and last:
            if not 'trace' in out:
                out['trace']=[]
            dopescore=self.loop.assess_methods[0](self.s2)[1] if last else 9999999
            soapscore=self.loop.assess_methods[1](self.s2)[1] if last else 9999999
            dopescore=999999999 if np.isnan(dopescore) else dopescore
            soapscore=999999999 if np.isnan(soapscore) else soapscore
            out['trace'].append((mcrmsd,dopescore,soapscore))                
        print r.rms
        return r.rms

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

        atmsel.randomize_xyz(deviation=4)

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
        if self.assess_trace:
            act.append(checkRMSD(20,self))
        return act
               
class sprefine(object):
    """
    SOAP loop refinement benchmark class
    """
    def __init__(self,dslist='101',bm='loop.100.8',criteria='bestrmsd',nonbond_spline=0.1,contact_shell=12.0,deviations=50,assess='',mcrmsdonly=False,spmodel=None,refineProtocal=None,trace=True,saveStructure=False):
        self.mcrmsdonly=mcrmsdonly
        if refineProtocal!=None:
            self.refineProtocal=refineProtocal
            dslist=refineProtocal['dslist']
            bm=refineProtocal['bm']
            nonbond_spline=refineProtocal['nonbond_spline']
            contact_shell=refineProtocal['contact_shell']
            deviations=refineProtocal['deviations']
            if 'assess' in refineProtocal:
                assess=refineProtocal['assess']
            if 'report' in refineProtocal and refineProtocal['report']=='mcrmsd':
                self.mcrmsdonly=True
            if 'trace' in refineProtocal:
                trace=refineProtocal['trace']
            if 'save' in refineProtocal:
                saveStructure=refineProtocal['save']
        #use statistical potential to refine loops/strutures, and analyze the refine results to judge the statistical potential.
        if spmodel: #needfix???
            if spmodel['scoretype']!='sprefine':
                raise Exception('You can not call sprefinescore with this type of model, wrong data type')
            bm=spmodel['bm']
            dslist=model['dslist']
            criteria=model['criteria']
        self.bm=bm
        self.dslist=dslist
        self.trace=trace
        self.assess_method=assess
        self.criteria=criteria
        self.saveStructure=saveStructure
        #self.codelist=decoysets(self.dslist).get_nativelist()  not applicalbel for loops    
        self.refpot=['$(LIB)/atmcls-mf.lib','$(LIB)/dist-mf.lib']
        bml=bm.split('.')
        self.refinetype=bml[0]
        self.refinemethod=bml[1:-1][0]# [restratins,optimizationmethod,initialmodel,randomizem,numofmodelstobuild]
        self.nofloops=int(bml[-1])
        self.startNum=0
        self.endNum=self.nofloops
        self.slavenumber=0
        self.codedict={}
        self.pe='local'
        self.calc_rmsds=True
        self.envinitialized=False
        self.j=[]
        self.nors=5000
        self.cv=None
        self.nonbond_spline=nonbond_spline
        self.contact_shell=contact_shell
        self.deviations=deviations
        
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
            self.codedict=pickle.load(fh)
        else:
            self.codedict={}
            if not isinstance(self.dslist,list):
                self.dslist=[self.dslist]
            for ds in self.dslist:
                fh=open(runenv.basedir+'loop/'+ds+'.pickle')
                self.codedict.update(pickle.load(fh))
        
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
                res=self.loop_model_single_protein(modelfile='pdb'+key+'.ent',loops=[loop], modelname=loopname)
                resultdict[loopname]=res
                    #print os.system('mv '+loopname+'.output.pickle loop'+str(self.input)+'.'+loopname+'.output.pickle')
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
            rs=-2        # directories for input atom files

        env = environ(rand_seed=rs)
        env.io.atom_files_directory = [os.path.join(runenv.loopStructurePath,self.dslist),scratchdir,runenv.opdbdir,'.']
        # Create a new class based on 'loopmodel' so that we can redefine
        # select_loop_atoms (necessary)
        if self.assess_method=='SOAP':
            from modeller import soap_loop
            loop_assess_methods=(assess.DOPE, soap_loop.Scorer())            
            energytrace=True
        else:
            loop_assess_methods=(assess.DOPE)
            energytrace=False
        if self.refinetype=='original':
            from modeller import soap_loop
            loop_assess_methods=(assess.DOPE, soap_loop.Scorer())
            m = ORLoop(env,
                   inimodel=modelfile, # initial model of the target
                   sequence=modelname,loop_assess_methods=loop_assess_methods,assess_trace=self.trace,
                   refinepot=self.refpot,loops=loops, nonbond_spine=self.nonbond_spline,contact_shell=self.contact_shell,deviations=self.deviations,energytrace=energytrace)
        else:
            m = MyLoop(env,
                   inimodel=modelfile, # initial model of the target
                   sequence=modelname,loop_assess_methods=loop_assess_methods,assess_trace=self.trace,
                   refinepot=self.refpot,loops=loops, nonbond_spine=self.nonbond_spline,contact_shell=self.contact_shell,deviations=self.deviations,energytrace=energytrace)
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
        #nativescore=m.ors.assess_dope()
        loopoutput=m.loop.outputs
        #pdb.set_trace()
        for i in range(0,len(loopoutput)):
            try:
                lastone='trace' if self.assess_method=='SOAP' else 'DOPE score'
                resultlist.append((loopoutput[i]['rmsd'],loopoutput[i]['mcrmsd'],loopoutput[i][lastone]))
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
        if self.cv!=None:
            self.cv.generate_potential()
            self.refpot[1]=os.path.join(runenv.serverUserPath,'lib',self.cv.rundir+'.lib')
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
        if self.assess_method=='SOAP':
            freememory=2.6
        else:
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
        numofruns=min(self.nors,tnt/20)
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
        makemdt=open('runme.py','w')#nonbond_spline=0.1,contact_shell=12.0,deviations=50
        makemdt.write('from SOAP.loop import *\nimport sys \n \nspopt=sprefine(sys.argv[1],"'+self.bm+'","'+self.criteria+'",'+str(self.nonbond_spline)+','+str(self.contact_shell)+','+str(self.deviations)+')\nspopt.runpath=\''+self.runpath+'\'\nspopt.refpot[1]="'+self.refpot[1]+'"\n\nspopt.assess_method="'+self.assess_method+'"\n\nspopt.saveStructure='+str(self.saveStructure)+'\n\nspopt.trace='+str(self.trace)+'\n\nspopt.initialize_dslist()'+'\n'+'\nspopt.assess_cluster_node()\n')
        makemdt.flush()
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
            om=copy.deepcopy(cv.model)
            om['refineProtocal']=self.refineProtocal
            self.model=om
            if self.cv.modelinlog(self.model):
                return 0
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
                if self.saveStructure:
                    fl=[f for f in os.listdir('./') if f.startswith(loopname+'.BL')]
                    print os.system('mkdir '+self.dslist)
                    for f in fl:
                        print os.system('mv '+f+' '+self.dslist+'/'+f)
                    print os.system('cd '+self.dslist+';gzip *;cd ..')
        return resultdict          

    def assess_cluster_node(self):
        os.chdir(scratchdir)
        result=self.model_loops_list(self.codedict)
        #if self.criteria=='bestrmsd':
        #    result=self.assess_bestrmsd(result)
        pickle.dump(result,open(self.dslist+'.pickle','w'))
        report_job_runstatus(self.runpath, True, self.dslist, '',inputname='runme.py',temppath=scratchdir)
        
    def afterprocessing(self,returnDetails=False):
        os.chdir(self.dir+self.rundir)
        rd={}
        print self.nrpl
        print self.nors
        for i in range(self.nors):
            if os.path.isfile(str(i+1)+'.pickle'):
                fh=open(str(i+1)+'.pickle','rb')
                res=pickle.load(fh)
                fh.close()
                if not res.keys()[0] in rd:
                    rd[res.keys()[0]]=[]
                #if self.assess_method=='SOAP':
                #    for item in res.values()[0]:
                #        rd[res.keys()[0]].extend(item[-1])
                rd[res.keys()[0]].extend(res.values()[0])
        #print os.system('rm -r '+self.dir+self.rundir)
        #del self.clusters
        self.result=rd
        #pdb.set_trace()
        try:
            del self.task
            gc.collect()
        except:
            traceback.print_exc()
        return self.analyze_loop_modeling(returnDetails)
    
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

    def analyze_loop_modeling(self,returnDetails=False):
        if self.assess_method=='SOAP':            
            result=self.result
            mcrmsds=[]
            mcrmsds2=[]
            dopemins=[]
            soapmins=[]
            ral=[]
            from itertools import chain
            if returnDetails:
                rrl=[]
                rl=[]
                dl=[]
                sl=[]
            rastart=-1
            for key in result:
                mcrmsds.append(min([item[1] for item in result[key]]))
                if returnDetails:
                    srrl=[result[key][0][1]]
                    sa0=np.array(result[key][0][-1][rastart:])
                    srl=[sa0[:,0].min()]
                    md=sa0[:,1].min()
                    ms=sa0[:,2].min()
                    sdl=[sa0[np.argmin(sa0[:,1]),0]]
                    ssl=[sa0[np.argmin(sa0[:,2]),0]]
                    if result[key][0]==result[key][1]:
                        ress=result[key][2::2]
                    else:
                        ress=result[key][1:]
                    for item in ress:
                        sa=np.array(item[-1][rastart:])
                        srrl.append(min(item[1],srrl[-1]))
                        srl.append(min(sa[:,0].min(),srl[-1]))                        
                        dmi=np.argmin(sa[:,1])
                        smi=np.argmin(sa[:,2])
                        if sa[dmi,1]<md:
                            md=sa[dmi,1]
                            sdl.append(sa[dmi,0])
                        else:
                            sdl.append(sdl[-1])    
                        if sa[smi,2]<ms:
                            ms=sa[smi,2]
                            ssl.append(sa[smi,0])
                        else:
                            ssl.append(ssl[-1])                            
                    rl.append(srl)
                    rrl.append(srrl)
                    dl.append(sdl)
                    sl.append(ssl)
                flist=list(chain.from_iterable([item[-1][rastart:] for item in result[key]]))
                ra=np.array(flist)
                ral.append(ra)
                mcrmsds2.append(ra[:,0].min())
                dopeMinInd=np.argmin(ra[:,1])
                soapMinInd=np.argmin(ra[:,2])
                dopemins.append(ra[dopeMinInd,0])
                soapmins.append(ra[soapMinInd,0])
                print ra.shape
            print mcrmsds
            print mcrmsds2
            print 'mcrmsd',np.mean(mcrmsds)
            print 'mcrmsd2',np.mean(mcrmsds2)
            print 'dope',np.mean(dopemins)
            print 'soap',np.mean(soapmins)
            self.averageRMSD=np.mean(mcrmsds)
            if returnDetails:
                rlen=min([len(r) for r in rrl])
                return np.array([r[:rlen] for r in rrl]),np.array([r[:rlen] for r in rl]),np.array([r[:rlen] for r in dl]),np.array([r[:rlen] for r in sl])
        else:
            mn=9999999999
            result=self.result
            mcrmsds=[]
            for key in result:
                mcrmsds.append(min([item[1] for item in result[key]]))
                if len(result[key])<mn:
                    mn=len(result[key])
            print "min n "+str(mn)
            #mn=1000
            ra=np.zeros([len(result.keys()),mn,3])
            k=0
            for key in result:
                ra[k,:,:]=np.array(result[key])[:mn,:]
                k+=1            
            self.averageRMSD=ra[:,:,1].min(axis=1).mean()
            if not self.mcrmsdonly:
                self.averageRMSD+=ra[:,:,0].min(axis=1).mean()
            print 'mcrmsd',np.mean(mcrmsds) 
        print self.averageRMSD
        if self.cv!=None:
            self.cv.resultsummary=-self.averageRMSD
            self.cv.resultarray[-1]=-self.averageRMSD
            self.cv.write2logshelve(self.model)
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

    def generate_refine_decoys(self,dsname=''):
        #dpath=os.path.join(runenv.decoysbasedir,'loop',dsname)
        dspath=os.path.join(runenv.localDecoysPath,'loop',dsname)
        os.makedirs(dspath)
        os.chdir(self.dir+self.rundir)
        rd={}
        print self.nrpl
        print self.nors
        for i in range(self.nors):
            if os.path.isfile(str(i+1)+'.pickle'):
                inputdict=pickle.load(gzip.open(str(i+1)+'.pickle.gz'))
                del inputdict['runenv']
                fh=open(str(i+1)+'.pickle','rb')
                res=pickle.load(fh)
                fh.close()
                l=inputdict.values()[0][0]
                code=inputdict.keys()[0]
                dsname=code+str(l[1])+('A' if l[0]=='' else l[0])+str(l[2])

                sdspath=os.path.join(dspath,dsname)
                if not os.path.isdir(sdspath):
                    os.makedirs(sdspath)
                    rmsddict={}
                    scoredict={}
                else:
                    rmsddict=pickle.load(open(os.path.join(sdspath,'rmsd.pickle')))
                    scoredict=pickle.load(open(os.path.join(sdspath,'score.pickle')))
                env=environ()
                env.io.atom_files_directory = [os.path.join(runenv.loopStructurePath,self.dslist),scratchdir,runenv.opdbdir,'.']
                basem=Mymodel(env)
                basem.load_model_file(os.path.join(runenv.opdbdir,code[1:3],'pdb'+code+'.ent.gz'))
                #basem.gen_atomlist()
                basem.gen_chainlist()
                cc=l[0] if len(l[0])>0 else 'A'
                basem.gen_base_file(cc,l[2:],os.path.join(sdspath,dsname+'.base'),'backup')
                fl=os.listdir(str(i+1))
                pdb.set_trace()
                for f,r in zip(fl,res.values()[0]):  
                    dpn='.'.join(f.split('.')[1:])
                    dn=dpn.split('.')[0]
                    m2=Mymodel(env)
                    m2.load_model_file(str(i+1)+'/'+f)
                    nl=basem.transfer_model_atoms(m2)
                    gzip.open(os.path.join(sdspath,dpn),'w').write('\n'.join(basem.get_part(cc,l[2:])))
                    rmsddict[dn]=[r[1],r[0]]
                    scoredict[dn]=r[2]
                    pdb.set_trace()
                                
                print os.system('touch '+os.path.join(sdspath,'needcombinewithbase'))
                pickle.dump(rmsddict,open(os.path.join(sdspath,'rmsd.pickle'),'w'))
                pickle.dump(scoredict,open(os.path.join(sdspath,'score.pickle'),'w'))
                pickle.dump([('mcrmsd','f4')],open(os.path.join(sdspath,'extrapar.pickle'),'w'))
        

def sgmd(atmsel, actions):
    """Very slow MD annealing"""
    sgmdrefine(atmsel, actions, cap=1.39, timestep=4.0,
           equil_its=300, equil_equil=20,
           equil_temps=(150.0, 250.0, 400.0, 700.0, 1000.0, 1300.0),
           sampl_its=1000, sampl_equil=200,
           sampl_temps=(1300.0, 1000.0, 800.0, 600.0, 500.0, 430.0, 370.0,
                        320.0, 300.0))

def sgmdrefine(atmsel, actions, cap, timestep, equil_its, equil_equil,
           equil_temps, sampl_its, sampl_equil, sampl_temps):
    from modeller.optimizers import molecular_dynamics

    mdl = atmsel.get_model()
    md = molecular_dynamics(cap_atom_shift=cap, md_time_step=timestep,
                            md_return='FINAL', output=mdl.optimize_output,
                            actions=actions,friction=0,guide_factor=0,guide_time=0)
    init_vel = True
    # First run for equilibration, the second for sampling:
    for (its, equil, temps) in ((equil_its, equil_equil, equil_temps),(sampl_its, sampl_equil, sampl_temps)):
        for temp in temps:
            md.optimize(atmsel, max_iterations=its, equilibrate=equil,
                        temperature=temp, init_velocities=init_vel)
            init_vel=False

def loopschedule():
    return schedule(4,
           [ step(CG, None, mk_scale(default=1.00, nonbond=0.0, spline=1.00)),
             step(CG, None, mk_scale(default=2.00, nonbond=0.01, spline=0.01)),
             step(CG, None, mk_scale(default=1.00, nonbond=0.10, spline=0.10)),
             step(CG, None, mk_scale(default=1.00, nonbond=0.50, spline=0.50)),
             step(CG, None, physical.values(default=4.00)) ])

class checkRMSD(actions.action):
    """calcualte the RMSDs of intermediate loops"""
    def __init__(self, skip, loopObj,first=False,
                 last=False, start=0):
        actions.action.__init__(self, skip, first, last)
        self.num = start
        self.loopObj=loopObj

    def __call__(self, opt):
        self.loopObj.user_after_single_loop_model(self.loopObj.loop.outputs[-1])
        self.num = self.num + 1
