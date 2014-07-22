"""
   SOAP RefineScore module

"""

from env import *


class sprefinescore(object):
    """
    Benchmark SOAP based on refinement of structures.
    """
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
        self.codelist=DecoySets(self.dslist).get_nativelist()

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
