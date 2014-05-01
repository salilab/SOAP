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
                                parameters='/netapp/sali/gqdong/sps/soap_loop/dist.lib')


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
            