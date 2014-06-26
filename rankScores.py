"""
   SOAP RankScore module, calculating scores on decoys.
   
   In SOAP, recovery function scores are treated as a separate scoring term, which enables a substaintial~1000000x increase in evaluation speed.

"""

from env import *
from recoveryFunction import *
from statsTable import *
from decoys import decoysets

                
class dsscore(object):
    """
    Base class for decoy sets scores
    """
    def __init__(self,scorename='',dir='',dslist=[]):
        self.scorename=scorename
        self.dir=dir
        self.scorepath=self.dir+self.scorename
        self.dslist=[]
        self.score=[]
        self.method=[]
        
    def __add__(self,other):
        ndss=self
        ndss.score=self.score+other.score
        ndss.scorename=self.scorename+other.scorename.replace('_'.join(self.dslist),'')
        ndss.scorepath=self.dir+self.scorename
        ndss.method=self.method+other.method
        return ndss
            
    def save(self):
        np.save(self.scorepath+'.numpy',self.score)
        #fh=open(self.scorepath+'.pickle','wb')
        #pickle.dump(self.dnlist,fh)
        #fh.close()
        #sio.savemat(self.scorepath+'.mat',{'score':self.score})
        
    def load(self):
        fh=open(self.scorepath+'.pickle','rb')
        no=cPickle.load(fh)
        fh.close()
        return no
    
    def get(self,nors=2000):
        if os.path.isfile(self.scorepath+'.pickle'):
            return self.load()
        else:
            self.calc_score_cluster(nors)
            return self
        
    def get_score(self,nors=2000):
        if len(self.score)>0:
            return self.score
        elif os.path.isfile(self.scorepath+'.numpy.npy'):
            self.score=np.load(self.scorepath+'.numpy.npy').astype(np.float32)
        else:
            self.calc_score_cluster(nors)
        return self.score

class sfdsscore(dsscore,sf): #the score comes from scaling function
    """
    Decoy sets scores calculated using recovery function
    """
    def __init__(self, dslist=[],bm='ss1', sfm='', features='', sftype='',par='',parvalue=[],type='',**argvs):
        dsscore.__init__(self,scorename='')
        sf.__init__(self,sftype=sftype,par=par,parvalue=parvalue,features=features,type=type, **argvs)
        self.dslist=dslist
        self.bm=bm
        self.scorename=sfm+'.'+bm+'.'+features+'.'+'_'.join(dslist)
        self.dir=runenv.basedir+'sfscore/'
        self.scorepath=self.dir+self.scorename
        self.idist=[]
        #self.method=[{'scoreclass':'sf','bm':bm, 'features':features, 'sfm':sfm,'dslist':self.dslist}]
        if sfm and sfm!='optimized':
            self.get_sfv()
        if len(self.dslist)>0:
            self.get_dsidist()
        
    def get_dsidist(self,dslist=[],bm='', features=''):
        if dslist:
            self.dslist=dslist
        if bm:
            self.bm=bm
        if features:
            self.features=features
        idist=rawsp('_'.join(self.dslist),self.features,self.bm,decoy=True,routine='calc_sp_i').get_i()#.astype(np.float32)     
        sdl=1
        ms=idist.shape
        for msi in ms[1:]:
            sdl=sdl*msi
        try:
            if idist.shape!=[ms[0],sdl]:
                self.idist=np.array(idist.reshape([ms[0],sdl]))
        except:
            pdb.set_trace()
        self.score=np.zeros(ms[0],dtype=np.float32)
    
    def filter(self,fa,newscore):
        self.idist=self.idist[fa,:]
        self.score=newscore
    
    def get_dsidistpath(self):
        distfeatures=re.search('(d[0-9]+)',self.features).group(1)
        return rawsp('_'.join(self.dslist),distfeatures,self.bm,decoy=True).sppath+'.npy'

    def get_dsidistname(self):
        distfeatures=re.search('(d[0-9]+)',self.features).group(1)
        return rawsp('_'.join(self.dslist),distfeatures,self.bm,decoy=True).spname+'.npy'
        
    def get_score(self,indexlist=[],sfv=[]):
        if len(sfv)==0:
            self.get_sf()
        else:
            self.sfv=sfv
        sfv=self.sfv.ravel()
        sfv=sfv.astype(np.float32)
        if len(indexlist)==0:
            np.dot(self.idist,sfv,self.score)
        else:
            for index in indexlist:
                try:
                    np.dot(self.idist[index[0]:index[1],:],sfv,self.score[index[0]:index[1]])
                except:
                    pdb.set_trace()
        if np.isnan(self.score).sum()>0:
            print sfv
            print "nan in score"
            raise NanInScore()
        return self.score
        
class spdsscore(dsscore,scaledsp): #the score from distributions
    """
    Decoy sets scores calculated using stats tables
    """    
    def __init__(self,dslist=[],bm='',ssp=scaledsp(),model=[]):
        dsscore.__init__(self)
        if model: #needfix???
            #if model['type']!='dsscore':
            #    raise Exception('You can not call spdsscore with this type of model, wrong data type')
            bm=model['bm']
            dslist=model['dslist']
            #del model['scoretype']
            del model['bm']
            del model['dslist']
            ssp=scaledsp(model=model)
        self.bm=bm
        self.dslist=dslist
        scaledsp.__init__(self,ssp.pm,ssp.rsp)
        #pdb.set_trace()
        self.scaledsp=ssp
        self.scorename=ssp.sspname+'.'+'_'.join(dslist)+'.'+bm
        self.scorepath=self.dir+self.scorename
        self.method=[{'scoreclass':'sp','pdbset':self.pdbset,  'features':self.features, 'genmethod':self.genmethod, 'pm':self.pm, 'bm':self.bm, 'dslist':self.dslist}]
        self.nors=2000
        self.test=False
        self.runpath='./'
        self.justwait=False
        self.isatompairclustering()
        self.isdecomposition()
        
    def isatompairclustering(self):
        genmethod=self.bm
        self.atompairclustering=False;
        self.apca=0
        #pdb.set_trace()
        if re.search('(ap[0-9]{1,2})',genmethod):
            rer=re.search('(ap[0-9]{1,2})',genmethod)
            apname=rer.group(1)
            self.atompairclustering=True
            self.apca=np.load(runenv.libdir+apname+'.npy')

    def isdecomposition(self):
        genmethod=self.bm
        self.decomposition=False;
        self.svdu=None
        self.svdv=None
        if genmethod.endswith('svd'):
            self.decomposition=True

    def prepare_task_input(self):
        os.chdir(self.dir)
        os.system('rm -r '+self.scorename)
        os.mkdir(self.scorename)
        os.chdir(self.dir+self.scorename)
        if self.decomposition:
            self.svdu,self.svdv=self.get_svd()
            print os.system('cp '+self.ssppath+'svdu.npy ./')
            print os.system('cp '+self.ssppath+'svdv.npy ./')
        if self.tlib:
            os.system('cp '+runenv.libdir+self.tlib+' ./')        
        freememory=feature(self.features).get_runmem()
        if self.features=='epad':
            freememory=3.0
        elif self.features=='plop':
            freememory=2.0
        elif self.features in ['soap_loop','soap_pep','soap_protein_od']:
            freememory=14
        #if freememory>8:
        #    self.nors=300
        if 'l' in self.fo.fclist:
            residuepower=1
        else:
            residuepower=2
        self.nors=2000
        pir(decoysets( self.dslist).pirpath).sep_pir('./',self.nors,permin=scoreminnumber,residuepower=residuepower)
        nors=int(os.popen('ls | grep -c pdbs').read())
        self.nors=nors
        self.runpath=runenv.serverUserPath+self.scorename+'/'
        #write the input list for different runs
        inputlist=open('inputlist','w')
        inputlist.write(','.join(['pdbs'+str(i) for i in range(1,nors+1)]))
        inputlist.close()
        if os.path.isfile(self.ssppath+'.hdf5'):
            os.system('ln -s '+self.ssppath+'.hdf5 ./'+self.sspname+'.hdf5')
        makemdt=open('runme.py','w')
        makemdt.write('keeptry=5\nwhile keeptry:\n    try:\n        from SOAP.rankScores import *\n        keeptry=0\n    except:\n        keeptry-=1\n'
                      +'import sys \n \nrsp=rawsp(\''+self.pdbset+'\',\''
                      +self.features+'\',\''+self.genmethod+'\')\n'+'ssp=scaledsp(\''+self.pm+'\',rsp)\n'
                      +'so=spdsscore('+str(self.dslist)+',\''+self.bm+'\',ssp)\n'
                      +'so.ssppath=\''+runenv.serverUserPath+self.scorename+'/'+self.sspname+'\'\n'
                      +'so.test='+str(self.test)+'\n'
                      +'so.calc_score(sys.argv[1])')
        makemdt.flush()
        if self.features in ['firedock','zrank','rdock']:
            freememory=2.0
        generate_job_submit_script(freememory,self.scorename,runtime,nors)
        return []

    def calc_score_cluster(self):
        sj=self.get_task()
        tasklist(tasklist=[sj]).monitor2end()
    
    def calc_score_local(self,nors):
        self.prepare_task_input(nors)
        sj=task(self.scorepath+'/',self.scorename)
        sj.run_task_local(self.calc_score)
        os.system('cat *score > allbenchmark')
        self.sort_score()
        self.save()
        
    def calc_score_nmrfile(self,fn):
        fh=open(fn)
        fc=fh.read()
        fh.close()
        os.mkdir(fn[:-4])
        os.chdir(fn[:-4])
        fl=fc.split('ENDMDL')
        for f,i in zip(fl[:-1],range(len(fl)-1)):
            fh=open(str(i)+'.pdb','w')
            fh.write(f)
            fh.close()
        self.calc_score_dir('./')

    def calc_score_dir(self,path):
        os.chdir(path)
        fl=os.listdir(path)
        self.calc_score_fromlist(fl)

    def calc_score_fromlist(self, dnlist):
        """Methods to benchmark the statistical potential, run on current directory"""
        bm=self.bm
        features=self.features
        mdtfile=self.ssppath+'.hdf5'    
        env = environ()
        log.minimal()
        io=env.io
        mlib=mdt.Library(env)
        mlib.bond_classes.read('${LIB}/bndgrp.lib')
        featurelist=feature(self.features,mlib).get_featurelist()
        m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
        print mdtfile
        m.read_hdf5(mdtfile)
        ms=m.shape
        a=alignment(env)
        ffile=open('score','w')
        csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(bm)  
        fp=0;
        tn=0;
        xp=0
        for k in range(0,len(dnlist)):
            mdl=model(env)
            mdl.read(file=dnlist[k])
            aln=alignment(env)
            aln.append_model(mdl,align_codes=dnlist[k],atom_files=dnlist[k])
            tn=tn+1
            try:
                sr=0;
                alncode=aln[0].code
                io.atom_files_directory=[scratchdir,'./']
                sr=m.open_alignment(aln,io=io)
                #pdb.set_trace()
                score=sr.sum(chain_span_range=(-csmax,-csmin,csmin,csmax),residue_span_range=(-kcut1,-kcut,kcut,kcut1),bond_span_range=(bs,bsm),disulfide=ssp)
                print score
                ffile.write(aln[0].code+', '+str(score)+'\n')
                fp+=1
                print alncode
            except Exception,e:
                traceback.print_exc()
                if re.search('fork failure',str(e)):
                    ffile.close()
                    print os.system('rm '+sfile)
                    raise Exception('fork failure, quiting')
                xp=1
            aln=0
            aln=alignment(env)
            mdl=[]
        print tn
        print fp
        print xp
        ffile.flush()

    def calc_special_score(self,env,fd,naln):
        pdbfile=fd+naln[0].atom_file
        if self.features=='dope':
            if not 'currentmodel' in self.__dict__.keys():
                mdl=model(env)                
                mdl.read(pdbfile+'.pdb')
            else:
                mdl=self.currentmodel
            return selection(mdl).assess_dope()
        elif self.features.startswith('soap'):
            if not 'currentmodel' in self.__dict__.keys():
                mdl=model(env)                
                mdl.read(pdbfile+'.pdb')
            else:
                mdl=self.currentmodel
            return selection(mdl).energy()[0]        
        else:
            raise Exception('Do not know how to calculate scores using potential '+self.features)

    def calc_score(self,runnumber,reportstatus=True):
        """Methods to benchmark the statistical potential, run on current directory"""
        ddbfile='pdbs'+str(runnumber)
        tmpdir=scratchdir
        bm=self.bm
        env = environ()
        if self.atompairclustering and self.decomposition:
            dnlist=pir(ddbfile).get_codelist()
            idist=self.calc_individual_sp(ddbfile,dnlist,'apdecompscore')
            np.save(str(runnumber)+'.score.npy',idist)
            if reportstatus:
                report_job_runstatus(self.runpath, True, runnumber, '.score.npy',inputname='runme.py')
            return idist
        elif self.atompairclustering:
            dnlist=pir(ddbfile).get_codelist()
            idist=self.calc_individual_sp(ddbfile,dnlist,'apscore')
            np.save(str(runnumber)+'.score.npy',idist)
            if reportstatus:
                report_job_runstatus(self.runpath, True, runnumber, '.score.npy',inputname='runme.py')
            return idist
        elif self.decomposition:
            dnlist=pir(ddbfile).get_codelist()
            idist=self.calc_individual_sp(ddbfile,dnlist,'decompscore')
            np.save(str(runnumber)+'.score.npy',idist)
            if reportstatus:
                report_job_runstatus(self.runpath, True, runnumber, '.score.npy',inputname='runme.py')
            return idist            
        if self.features=='soap_loop':
            from modeller import soap_loop
            pair = soap_loop.Scorer()
            env.edat.sphere_stdv=-1
            env.edat.contact_shell=6
            env.edat.dynamic_sphere=False
            env.edat.energy_terms.append(pair)
        elif self.features=='soap_pep':
            from modeller import soap_peptide
            pair = soap_peptide.Scorer('soap_pep.hdf5')
            env.edat.sphere_stdv=-1
            env.edat.contact_shell=6
            env.edat.dynamic_sphere=False
            env.edat.energy_terms.append(pair)
        elif self.features=='soap_pp':
            from modeller import soap_pp
            atom = soap_pp.AtomScorer()
            pair = soap_pp.PairScorer()
            env.edat.sphere_stdv=-1
            env.edat.contact_shell=15
            env.edat.dynamic_sphere=False
            env.edat.energy_terms.append(atom)            
            env.edat.energy_terms.append(pair)
        elif self.features=='soap_protein_od':
            from modeller import soap_protein_od
            pair = soap_protein_od.Scorer()
            env.edat.sphere_stdv=-1
            env.edat.contact_shell=7
            env.edat.dynamic_sphere=False
            env.edat.energy_terms.append(pair)            
        features=self.features
        mdtfile=self.ssppath+'.hdf5'    
        print 'Benchmarking: '+ ddbfile
        log.minimal()
        io=env.io
        mlib=mdt.Library(env)
        mlib.bond_classes.read('${LIB}/bndgrp.lib')
        featurelist=feature(self.features,mlib).get_featurelist()
        if self.pdbset!='Otherpot':
            m=mdt.Table(mlib,features=featurelist,shape=[-1]*len(featurelist))
            print mdtfile
            m.read_hdf5(mdtfile)
        a=alignment(env)
        f=modfile.File(self.runpath+ddbfile,'r')
        csmin, csmax,kcut, kcut1, bs,bsm, errsv,ssp=decode_genmethod(bm) 
        fp=0;
        tn=0;
        xp=0
        scorestr=''
        aln=alignment(env)
        tmml=[]
        fdl=[]
        while aln.read_one(f):
            tn=tn+1
            try:
                sr=0;
                alncode=aln[0].code
                codelist=alncode.split('.')
                spdbdir=runenv.ddbdir+codelist[0]+'/'+codelist[1]+'/'
                if self.features in ['firedock','zrank','rdock','dfire','epad','rosetta','plop','posescore','rankscore']:
                    naln=aln
                    fd=spdbdir
                    if len(fdl)==0 or fd!=fdl[-1]:
                        fdl.append(fd)
                        tmml.append([])
                    if self.features in ['dfire','epad']:
                        if not os.path.isfile('.'.join(alncode.split('.')[:2])+'.seq'):
                            aln.write(file='.'.join(alncode.split('.')[:2])+'.seq', alignment_format='FASTA')
                        tmml[-1].append(alncode)
                    else:
                        tmml[-1].append(alncode)
                else:
                    errcount=1
                    while errcount<10:
                        try:
                            fd,naln=self.get_pdbfile(spdbdir,codelist,alncode,env,aln,tmpdir)                    
                            io.atom_files_directory=[fd]
                            if self.features=='dope' or self.features.startswith('soap'):
                                score=self.calc_special_score(env,fd,naln)
                            else:
                                sr=m.open_alignment(naln,io=io)
                                score=sr.sum(chain_span_range=(-csmax,-csmin,csmin,csmax),residue_span_range=(-kcut1,-kcut,kcut,kcut1),bond_span_range=(bs,bsm),disulfide=ssp)
                            print score
                            scorestr=scorestr+'>'+aln[0].code+' score: '+str(score)+'\n'
                            fp+=1
                            print alncode
                            break
                        except Exception,e:
                            time.sleep(0.1)
                            print open(os.path.join(tmpdir,alncode+'.pdb')).read()
                            print os.system('rm -f '+os.path.join(tmpdir,'*'))
                            errcount+=1
                            if errcount>=10:
                               raise e 
            except Exception,e:
                traceback.print_exc()
                xp=1
                if re.search('fork failure',str(e)):
                    raise Exception('fork failure, quiting')
                else:
                    if not self.test:
                        raise e                 
            aln=0
            aln=alignment(env)
        print tn
        print fp
        print xp
        if self.features in ['firedock','zrank','rdock']:
            scorestr=self.calc_dock_score(runnumber,tmml,fdl)
            tn=fp
        elif self.features in ['posescore','rankscore']:
            scorestr=self.calc_poserank_score(runnumber,tmml,fdl)
            tn=fp
        elif self.features=='dfire':
            scorestr=self.calc_dfire_score(runnumber,tmml,fdl)
            tn=fp
        elif self.features=='rosetta':
            scorestr=self.calc_rosetta_score(runnumber,tmml,fdl)
            tn=fp
        elif self.features=='plop':
            scorestr=self.calc_plop_score(runnumber,tmml,fdl)
            tn=fp              
        elif self.features=='epad':
            scorestr=self.calc_epad_score(runnumber,tmml,fdl)
            tn=fp
        if self.test or (tn==fp and len(scorestr)>0):
            sfile=str(runnumber)+'.score'
            try:
                ffile=open(sfile,'w')
                ffile.write(os.path.join(self.runpath,scorestr))
                ffile.flush()
                runsuccess=True
            except Exception, e:
                print traceback.print_exc()
                runsuccess=False
                print 'FatalError: can not save run output file'
            print 'Benchmarking finished: '+ ddbfile
        else:
            runsuccess=False
            print "FatalError: encountered when benchmarking using mdt "+ddbfile
        if reportstatus:
            report_job_runstatus(self.runpath, runsuccess, runnumber, '.score',inputname='runme.py')

    def calc_poserank_score(self,runnumber,tmml,fdl):
        wd=os.getcwd()
        scorestr=''
        os.chdir(scratchdir)
        for i in range(len(fdl)):
            fd=fdl[i]
            for sml in tmml[i]:
                smll=sml.split('.')
                print os.system('cp '+os.path.join(runenv.ddbdir,smll[0],smll[1],sml+'.pdb.gz')+' '+scratchdir)
                print os.system('gunzip '+os.path.join(scratchdir,sml+'.pdb.gz'))
                print os.system('/netapp/sali/openeye/wrappers/v2012.Oct.1/python/examples/oechem/convert.py '+sml+'.pdb '+sml+'.mol2')
                if self.features =='posescore':
                    print os.system('/netapp/sali/haofan/imp/bin/imppy.sh python '+runenv.serverUserPath+'/gst/test_PLscore.py '+os.path.join(runenv.ddbdir,smll[0],smll[1],smll[0]+'.'+smll[1]+'.base.pdb')+' '+sml)
                elif self.features =='rankscore':
                    print os.system('/netapp/sali/haofan/imp/bin/imppy.sh python '+runenv.serverUserPath+'/gst/rankscore.py '+os.path.join(runenv.ddbdir,smll[0],smll[1],'base.pdb')+' '+sml)
                #print os.system('rm '+os.path.join(scratchdir,sml+'.pdb*'))
                fh3=open(os.path.join(scratchdir,sml+'.score'))
                fc=fh3.read()
                fh3.close()
                #rer=re.findall(' DFIRE-score :\s*([0-9\.\-\s]+)',fc)
                #for item in rer:
                try:
                    scorestr=scorestr+'>'+sml+' score: '+fc.split()[-1]+'\n'
                except:
                    scorestr=scorestr+'>'+sml+' score: '+'999.0'+'\n'
        print os.system('rm -r '+scratchdir)
        os.chdir(wd)
        return scorestr

    def calc_dfire_score(self,runnumber,tmml,fdl):
        wd=os.getcwd()
        scorestr=''
        for i in range(len(fdl)):
            fd=fdl[i]
            for sml in tmml[i]:
                if not os.path.isfile(os.path.join(scratchdir,sml+'.pdb')):
                    print os.system('cp '+os.path.join(fd,sml+'.pdb.gz')+' '+scratchdir)
                    print os.system('gunzip '+os.path.join(scratchdir,sml+'.pdb.gz'))
                print os.system('/netapp/sali/haofan/score/rapdf/know -p /netapp/sali/haofan/score/rapdf/pmf1_dfire151_bin.dat -f '+os.path.join(scratchdir,sml+'.pdb')+' > '+os.path.join(scratchdir,sml+'.score'))
                #print os.system('rm '+os.path.join(scratchdir,sml+'.pdb*'))
                fh3=open(os.path.join(scratchdir,sml+'.score'))
                fc=fh3.read()
                fh3.close()
                #rer=re.findall(' DFIRE-score :\s*([0-9\.\-\s]+)',fc)
                #for item in rer:
                scorestr=scorestr+'>'+sml+' score: '+fc.split()[-1]+'\n'
        print os.system('rm -r '+scratchdir)    
        return scorestr

    def calc_rosetta_score(self,runnumber,tmml,fdl):     
        wd=os.getcwd()
        scorestr=''
        for i in range(len(fdl)):
            fd=fdl[i]
            for sml in tmml[i]:
                print os.chdir(scratchdir)
                if not os.path.isfile(os.path.join(scratchdir,sml+'.pdb')):
                    print os.system('cp '+os.path.join(fd,sml+'.pdb.gz')+' '+scratchdir)
                    print os.system('gunzip '+os.path.join(scratchdir,sml+'.pdb.gz'))
                if self.genmethod=='fpd':
                    print os.system(runenv.serverUserPath+'rosetta_source/bin/FlexPepDocking.linuxgccrelease -database '+runenv.serverUserPath+'/rosetta_source/rosetta_database/ -flexpep_score_only -s '+os.path.join(scratchdir,sml+'.pdb')+' -out:file:scorefile '+os.path.join(scratchdir,sml+'.score'))                    
                else:
                    print os.system(runenv.serverUserPath+'/rosetta_source/bin/score.linuxgccrelease -database '+runenv.serverUserPath+'/rosetta_source/rosetta_database/ -s '+os.path.join(scratchdir,sml+'.pdb')+' -out:file:scorefile '+os.path.join(scratchdir,sml+'.score'))
                #print os.system('rm '+os.path.join(scratchdir,sml+'.pdb*'))
                fh3=open(os.path.join(scratchdir,sml+'.score'))
                fc=fh3.read()
                fh3.close()
                #rer=re.findall(' DFIRE-score :\s*([0-9\.\-\s]+)',fc)
                #for item in rer:
                if self.genmethod=='fpd':
                    scorestr=scorestr+'>'+sml+' score: '+fc.split('\n')[2].split()[1]+'\n'
                else:
                    scorestr=scorestr+'>'+sml+' score: '+fc.split('\n')[1].split()[1]+'\n'
        print os.chdir(wd)
        print os.system('rm -r '+scratchdir)    
        return scorestr

    def calc_plop_score(self,runnumber,tmml,fdl):
        wd=os.getcwd()
        scorestr=''
        for i in range(len(fdl)):
            fd=fdl[i]
            for sml in tmml[i]:
                if not os.path.isfile(os.path.join(scratchdir,sml+'.pdb')):
                    print os.system('cp '+os.path.join(fd,sml+'.pdb.gz')+' '+scratchdir)
                    print os.system('gunzip '+os.path.join(scratchdir,sml+'.pdb.gz'))
                ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent sgbnp penalty yes rbuffer 2.0\nminim res all res\n\nenergy calc')
                open(os.path.join(scratchdir,'ploop'),'w').write(ploop)
                print os.system('/netapp/home/ck/plop/plop '+os.path.join(scratchdir,'ploop')+' > '+os.path.join(scratchdir,sml+'.score'))
                #print os.system('rm '+os.path.join(scratchdir,sml+'.pdb*'))
                fh3=open(os.path.join(scratchdir,sml+'.score'))
                fc=fh3.read()
                fh3.close()
                #rer=re.findall(' DFIRE-score :\s*([0-9\.\-\s]+)',fc)
                #for item in rer:
                if re.search('TOTALE\s*([0-9.-]*)',fc):
                    scorestr=scorestr+'>'+sml+' score: '+re.findall('TOTALE\s*([0-9.-]*)',fc)[-1]+'\n'
                else:
                    scorestr=scorestr+'>'+sml+' score: 9999999.999998\n'
                print scorestr
        print os.system('rm -r '+scratchdir)    
        return scorestr

    def calc_plop_score_old(self,runnumber,tmml,fdl):
        wd=os.getcwd()
        scorestr=''
        for i in range(len(fdl)):
            fdname=tmml[i][0].split('.')[1]
            #fdname=fd.split('/')[-2]
            rer=re.search('([0-9]{1,2})([A-Z\-]{1})([0-9]{1,4})',fdname[4:]).groups()
            if fdname=='1a2y7A26':
                decoyfile=open(runenv.serverUserPath+'/loop/'+fdname).read()
            else:
                decoyfile=open(runenv.serverUserPath+'/loop/'+rer[0]+'res_decoys/'+fdname[:4]+rer[0]+'-'+rer[2]).read()
            dfl=decoyfile.split('ENDMDL')
            dflines=[]
            dfdict={}
            k=0
            for line in dfl[0].split('\n'):
                if line.startswith('ATOM'):
                    dflines.append(line)
                    dfdict[line[:26]]=k
                    k+=1
            cc=rer[1]
            if cc=='-':
                cc='_' 
            loopdef=cc+':'+rer[2]+' '+cc+':'+str(int(rer[0])+int(rer[2])-1)
            for sml in tmml[i]:
                #if not os.path.isfile(os.path.join(scratchdir,sml+'.pdb')):
                    #    print os.system('cp '+os.path.join(fd,sml+'.pdb.gz')+' '+scratchdir)
                #    print os.system('gunzip '+os.path.join(scratchdir,sml+'.pdb.gz'))
                modelnum=sml.split('.')[2]
                if modelnum=='native':
                    open(os.path.join(scratchdir,sml+'.pdb'),'w').write(''.join(dfl[0]))
                else:
                    cm=dfl[int(modelnum)].split('\n')
                    #if not modelnum in ''.join(cm[:5]):
                    #    raise Exception('errmor in model number'+sml)
                    for line in cm[2:]:
                        if not (line.startswith('ATOM') or len(line)>50):
                            continue
                        if line[:26] in dfdict:
                            dflines[dfdict[line[:26]]]=line
                        else:
                            continue
                    open(os.path.join(scratchdir,sml+'.pdb'),'w').write('\n'.join(dflines))
                if self.genmethod=='dvd2':
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent vdgbnp penalty no rbuffer 2.0\nenergy lipocalc group heavy '+loopdef+'\n\n\nenergy calc')
                elif self.genmethod=='ds2':
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent sgbnp penalty no rbuffer 2.0\nenergy lipocalc group heavy '+loopdef+'\n\n\nenergy calc')
                elif self.genmethod=='dsp2':
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent sgbnp penalty yes rbuffer 2.0\nenergy lipocalc group heavy '+loopdef+'\n\n\nenergy calc')
                elif self.genmethod=='dvdm':
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent vdgbnp penalty no rbuffer 2.0\nminim res all res\n\nenergy calc')
                elif self.genmethod=='dsm':
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent sgbnp penalty no rbuffer 2.0\nminim res all res\n\nenergy calc')
                elif self.genmethod.startswith('dspm'):
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent sgbnp penalty yes rbuffer 2.0\nminim res all res\n\nenergy calc')
                elif self.genmethod.startswith('dvdmc'):
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent vdgbnp penalty no rbuffer 2.0\nenergy lipocalc group heavy '+loopdef+'\nminim res all res\n\nenergy calc')
                elif self.genmethod.startswith('dsmc'):
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent sgbnp penalty no rbuffer 2.0\nenergy lipocalc group heavy '+loopdef+'\nminim res all res\n\nenergy calc')
                elif self.genmethod.startswith('dspmc'):
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent sgbnp penalty yes rbuffer 2.0\nenergy lipocalc group heavy '+loopdef+'\nminim res all res\n\nenergy calc')
                elif self.genmethod=='n2':
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent vdgbnp penalty no rbuffer 2.0\nenergy lipocalc group heavy '+loopdef+'\nenergy calc')
                elif self.genmethod=='nndefaultmin':
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent vdgbnp penalty no rbuffer 2.0\nenergy lipocalc group heavy '+loopdef+'\n\nminim res all res\nenergy calc')
                else:
                    if 'sgb' in self.genmethod:
                        gb='sgb'
                    else:
                        gb='vdgb'
                    if 'penalty' in self.genmethod:
                        pn='yes'
                    else:
                        pn='no'                        
                    ploop='file datadir /netapp/home/ck/plop/data\n\nload pdb '+os.path.join(scratchdir,sml+'.pdb \nenergy para solvent '+gb+'np penalty '+pn+' rbuffer 2.0\nenergy lipocalc group heavy '+loopdef+'\n\nenergy calc')
                open(os.path.join(scratchdir,'ploop'),'w').write(ploop)
                print os.system('/netapp/home/ck/plop/plop '+os.path.join(scratchdir,'ploop')+' > '+os.path.join(scratchdir,sml+'.score'))
                #print os.system('rm '+os.path.join(scratchdir,sml+'.pdb*'))
                fh3=open(os.path.join(scratchdir,sml+'.score'))
                fc=fh3.read()
                fh3.close()
                #rer=re.findall(' DFIRE-score :\s*([0-9\.\-\s]+)',fc)
                #for item in rer:
                if re.search('TOTALE\s*([0-9.-]*)',fc):
                    scorestr=scorestr+'>'+sml+' score: '+re.findall('TOTALE\s*([0-9.-]*)',fc)[-1]+'\n'
                else:
                    scorestr=scorestr+'>'+sml+' score: 9999999.999998\n'
                print scorestr
        print os.system('rm -r '+scratchdir)    
        return scorestr

    def calc_epad_score(self,runnumber,tmml,fdl):
        wd=os.getcwd()
        scorestr=''
        os.chdir(scratchdir)
        print os.system('cp '+runenv.serverUserPath+'sps/epad/EPAD/*sh ./')
        print os.system('cp -r '+runenv.serverUserPath+'sps/epad/EPAD/bin ./')
        print os.system('cp -r '+runenv.serverUserPath+'sps/epad/EPAD/config ./')
        print os.system('cp -r '+runenv.serverUserPath+'sps/epad/EPAD/SEQ ./')
        #os.mkdir(os.path.join(scratchdir,'SEQ'))
        os.mkdir(os.path.join(scratchdir,'DECOYS'))        
        #print os.system('mv *.seq '+os.path.join(scratchdir,'SEQ'))
        #fl=os.listdir(os.path.join(scratchdir,'SEQ'))
        #fl=[f[:-4] for f in fl if f.endswith('.seq')]
        for i in range(len(fdl)):
            os.chdir(scratchdir)
            fd=fdl[i]
            cname='.'.join(tmml[i][0].split('.')[:2])
            print os.system('cp '+os.path.join(wd,cname+'.seq')+' ./SEQ/')
            cdir=os.path.join(scratchdir,'DECOYS',cname)
            os.makedirs(cdir)
            pdblist=[]
            for sml in tmml[i]:
                pdblist.append(sml+'.pdb')
                if os.path.isfile(os.path.join(scratchdir,sml+'.pdb')):
                    print os.system('cp '+os.path.join(fd,sml+'.pdb')+' '+cdir)
                else:
                    print os.system('cp '+os.path.join(fd,sml+'.pdb.gz')+' '+cdir)
            open(os.path.join(scratchdir,'DECOYS',cname+'.lst'),'w').write('\n'.join(pdblist))
            print os.system('gunzip '+os.path.join(cdir,'*.pdb.gz'))
            print os.system('./generateEPADfeature.sh '+cname)
            print os.system('./runEPAD_example.sh '+cname)
            scorelines=open(cname+'.EPAD').read().split('\n')
                #rer=re.findall(' DFIRE-score :\s*([0-9\.\-\s]+)',fc)
            for item in scorelines:
                if len(item)<6:
                    continue
                scorestr=scorestr+'>'+item.split()[0]+' score: '+item.split()[2]+'\n'
        print os.system('rm -r '+scratchdir)
        os.chdir(wd)
        return scorestr

    def calc_dock_score(self,runnumber,tmml,fdl):
        wd=os.getcwd()
        rdir='/scratch/'+runenv.userName+'/'+str(time.time())
        print os.system('mkdir -p '+rdir)
        os.chdir(rdir)
        scorestr=''
        for fd,sml in zip(fdl,tmml):
            print os.system('rm *')
            fh2=open(fd+'firedockinput.pickle')
            fid=pickle.load(fh2)
            fh2.close()
            print os.system('cp '+os.path.join(fd,fid['rpdb'])+' ./')
            print os.system('cp '+os.path.join(fd,fid['lpdb'])+' ./')            
            if self.features=='firedock':
                sml=[str(i)+' '+sml[i] for i in range(len(sml))]
                inputfile='\n'.join(sml)
                fh=open('transfile','w')
                fh.write(inputfile)
                fh.close()                        
                print os.system('/netapp/sali/dina/FiberDock1.1/buildFiberDockParams.pl '+fid['rpdb']+' '+fid['lpdb']+' U U '+fid['complextype']+' transfile out 0 50 0.8 0 glpk outp > /dev/null')
                print os.system('/netapp/sali/dina/FiberDock1.1/FiberDock outp')
                fh3=open('out.ref')
                fc=fh3.read()
                fh3.close()
                rer=re.findall('\n(.*)\t\|([0-9\.\-\s]+)\|.*',fc)
                for item in rer:
                    scorestr=scorestr+'>'+item[0].strip()+' score: '+item[1].strip()+'\n'
            elif self.features=='zrank':
                print os.system('cp '+os.path.join(fd,fid['orpdbAB'])+' ./')
                print os.system('cp '+os.path.join(fd,fid['olpdbAB'])+' ./')               
                fl=''
                for tm, i in zip(sml, range(len(sml))):
                    print os.system(runenv.serverUserPath+'bin/pdb_trans64 '+tm+' <'+fid['lpdb']+' >ltemp.pdb')
                    print os.system('cat '+fid['rpdb']+' ltemp.pdb > '+str(i)+'.pdb')
                    fl=fl+str(i)+'.pdb\n'
                fh=open('pdblist','w')
                fh.write(fl)
                fh.close()
                print os.system(runenv.serverUserPath+'bin/zrank_linux_64bit/zrank pdblist >/dev/null')
                fh=open('pdblist.zr.out')
                for line in fh:
                    if len(line)>4:
                        item=line.split()
                        scorestr=scorestr+'>'+item[0].strip()+' score: '+item[1].strip()+'\n'
                fh.close()
            elif self.features=='rdock':
                fl=''           
                for tm, i in zip(sml, range(len(sml))):
                    print os.system(runenv.serverUserPath+'bin/pdb_trans64 '+tm+' <'+fid['lpdb']+' >ltemp.pdb')
                    if fid['lpdb'][0:4] in ['1ACB','2J7P','1HE8','1YVB','1I4D','1ZM4']:
                        print os.system('cat  ltemp.pdb '+fid['rpdb']+' > '+str(i)+'.pdb')
                    else:
                        print os.system('cat '+fid['rpdb']+' ltemp.pdb > '+str(i)+'.pdb')
                    fl=fl+str(i)+'.pdb\n'
                fh=open('pdblist','w')
                fh.write(fl)
                fh.close()
                print os.system(runenv.serverUserPath+'rosetta_source/bin/score.linuxgccrelease -database '+runenv.serverUserPath+'rosetta_source/rosetta_database -s '+fid['lpdb']+' -score:weights docking') 
                print os.system(runenv.serverUserPath+'rosetta_source/bin/score.linuxgccrelease -database '+runenv.serverUserPath+'rosetta_source/rosetta_database -s '+fid['rpdb']+' -score:weights docking')
                fh=open('default.sc')
                fh.readline()
                sepscore=0
                for line in fh:
                    if len(line)>20:
                        item=line.split()
                        sepscore+=float(item[1].strip())
                fh.close()
                print os.system('rm default.sc')
                print os.system(runenv.serverUserPath+'rosetta_source/bin/score.linuxgccrelease -database '+runenv.serverUserPath+'rosetta_source/rosetta_database -in:file:l pdblist -score:weights docking >/dev/null')
                fh=open('default.sc')
                fh.readline()
                for line in fh:
                    if len(line)>20:
                        item=line.split()
                        scorestr=scorestr+'>'+item[1].strip()+' score: '+str(float(item[1].strip())-sepscore)+'\n'
                fh.close()
        print os.system('rm -r '+rdir)
        os.chdir(wd)
        return scorestr

    def complete_pdb(self,pdbfile):
        env=environ()
        env.io.hydrogen=True
        env.libs.topology.read('${LIB}/top_heav.lib')
        env.libs.parameters.read('${LIB}/par.lib')
        m = complete_pdb(env, pdbfile)
        m.write(file=pdbfile,no_ter=True)
        
    def afterprocessing(self, benchmarkfile=''):
    #Sort benchmark scores according to the name. The order will be a standard for combining scores and other analysis...
        time.sleep(1)
        if self.justwait:
            return 1
        elif self.atompairclustering and self.decomposition:
            allscore=self.cat_spi(name='score')
            for j in range(allscore.shape[2]):
                for i in range(allscore.shape[1]):
                    sdir=runenv.basedir+self.pdbset+'/'+self.features+'-'+self.genmethod+'/'
                    scorepath=sdir+self.pdbset+'.'+self.features+'.'+self.genmethod+'.'+self.pm.replace('+','.')+'.'+'_'.join(self.dslist)+'.'+self.bm[:-3]+'n'+str(i)+'svd'+str(j)
                    if not os.path.isdir(sdir):
                        os.mkdir(sdir)
                    np.save(scorepath+'.numpy',allscore[:,i].squeeze())
            self.score=allscore
            self.save()
            print os.system('rm -r '+self.dir+self.scorename)
            return 1            
        elif self.atompairclustering:
            #process atom pair features
            allscore=self.cat_spi(name='score')
            for i in range(allscore.shape[1]):
                sdir=runenv.basedir+self.pdbset+'/'+self.features+'-'+self.genmethod+'/'
                scorepath=sdir+self.pdbset+'.'+self.features+'.'+self.genmethod+'.'+self.pm.replace('+','.')+'.'+'_'.join(self.dslist)+'.'+self.bm+'n'+str(i)
                if not os.path.isdir(sdir):
                    os.mkdir(sdir)
                np.save(scorepath+'.numpy',allscore[:,i].squeeze())
            self.score=allscore
            self.save()
            print os.system('rm -r '+self.dir+self.scorename)
            return 1
        elif self.decomposition:
            #process atom pair features
            allscore=self.cat_spi(name='score')
            for i in range(allscore.shape[1]):
                sdir=runenv.basedir+self.pdbset+'/'+self.features+'-'+self.genmethod+'/'
                scorepath=sdir+self.pdbset+'.'+self.features+'.'+self.genmethod+'.'+self.pm.replace('+','.')+'.'+'_'.join(self.dslist)+'.'+self.bm+str(i)
                if not os.path.isdir(sdir):
                    os.mkdir(sdir)
                np.save(scorepath+'.numpy',allscore[:,i].squeeze())
            self.score=allscore
            self.save()
            print os.system('rm -r '+self.dir+self.scorename)
            return 1
        else:
            os.chdir(self.scorepath+'/')
            scorelist=' '.join([str(i)+'.score' for i in range(1,self.nors+1)])
            os.system('cat '+scorelist+' > allbenchmark')
            if not benchmarkfile:
                benchmarkfile=self.scorepath+'/allbenchmark'
            sfile=open(benchmarkfile,'r')
            #build a dictionary uisng the benchmark score
            scorelist=[]
            dnlist=[]
            i=0;
            scorere=re.compile('>(\S*)\sscore:\s([-.0-9a-z]*)')
            for scoreline in sfile:
                scoreg=scorere.search(scoreline)
                if scoreg:
                    scorelist.append([scoreg.group(1),scoreg.group(2)])
                    i+=1
            if len(scorelist)==0:
                raise Exception('The score files is empty something is wrong with the code...')
            #sortedlist=sorted(scorelist, key=itemgetter(0))
            scorea=np.zeros([len(scorelist)])
            for i in range(0,len(scorelist)):
                if 'nan' in scorelist[i][1]:
                    scorea[i]=999999999999.0
                else:
                    try:
                        scorea[i]=float(scorelist[i][1])
                    except:
                        scorea[i]=999999999999.0
                        continue
                        #pdb.set_trace()
                dnlist.append(scorelist[i][0])
            self.score=scorea
            self.scorelist=scorelist
            self.dnlist=dnlist
            self.save()
            print os.system('rm -r '+self.dir+self.scorename)
            return 1
    
    def read_singlescorefile(self,scorefile):
        sfile=open(scorefile,'r')
        #build a dictionary uisng the benchmark score
        scorelist=[]
        #dnlist=[]
        i=0;
        scorere=re.compile('>(\S*)\sscore:\s([-.0-9]*)')
        for scoreline in sfile:
            scoreg=scorere.search(scoreline)
            if scoreg:
                scorelist.append(float(scoreg.group(2)))
                #dnlist.append(scoreg.group(1))
                i+=1
        return np.array(scorelist)
    
    def get_task(self):
        if os.path.isfile(self.scorepath+'.numpy.npy'):
            return 0
        else:
            if re.search('(ap[0-9]+n[0-9]+)(svd)([0-9]+)$',self.bm):
                rer=re.search('(ap[0-9]+)(n[0-9]+)(svd)([0-9]+)$',self.bm)
                nbm=self.bm[:-(len(rer.group(2))+len(rer.group(3))+len(rer.group(4)))]+'svd'
                if any2str([nbm,self.dslist,self.scaledsp.ssppath]) in runenv.spdsscoreTasks:
                    return runenv.spdsscoreTasks[any2str([nbm,self.dslist,self.scaledsp.ssppath])]                
                nsp=spdsscore(bm=nbm,dslist=self.dslist,ssp=self.scaledsp)
                newtask=nsp.get_task()
                runenv.spdsscoreTasks[any2str([nbm,self.dslist,self.scaledsp.ssppath])]=newtask
                return newtask                
            if re.search('(ap[0-9]+n[0-9]+)$',self.bm):
                rer=re.search('(ap[0-9]+)(n[0-9]+)$',self.bm)
                nbm=self.bm[:-len(rer.group(2))]
                if any2str([nbm,self.dslist,self.scaledsp.ssppath]) in runenv.spdsscoreTasks:
                    return runenv.spdsscoreTasks[any2str([nbm,self.dslist,self.scaledsp.ssppath])]                
                nsp=spdsscore(bm=nbm,dslist=self.dslist,ssp=self.scaledsp)
                newtask=nsp.get_task()
                runenv.spdsscoreTasks[any2str([nbm,self.dslist,self.scaledsp.ssppath])]=newtask
                return newtask
            if re.search('(svd[0-9]+)$',self.bm):
                rer=re.search('(svd)([0-9]+)$',self.bm)
                nbm=self.bm[:-len(rer.group(2))]
                if any2str([nbm,self.dslist,self.scaledsp.ssppath]) in runenv.spdsscoreTasks:
                    return runenv.spdsscoreTasks[any2str([nbm,self.dslist,self.scaledsp.ssppath])]
                nsp=spdsscore(bm=nbm,dslist=self.dslist,ssp=self.scaledsp)
                newtask=nsp.get_task()
                runenv.spdsscoreTasks[any2str([nbm,self.dslist,self.scaledsp.ssppath])]=newtask
                return newtask   
            if len(self.dslist)==1:  
                sj=task(self.scorepath+'/',self.scorename,afterprocessing=self,preparing=self,targetfile=self.scorepath+'.numpy.npy')
                return sj
            else:
                tl=[]
                for pset in self.dslist:
                    nspo=spdsscore([pset],self.bm,self.scaledsp)
                    tl.append(nspo.get_task())
                return tasklist(tasklist=tl,afterprocessing=dummy)
                
    def get_score(self,nors=2000):
        if len(self.score)>0:
            return self.score
        else:
            if len(self.dslist)==1:
                if os.path.isfile(self.scorepath+'.numpy.npy'):
                    self.score=np.load(self.scorepath+'.numpy.npy').astype(np.float32)
                else:
                    self.calc_score_cluster()
                    self.get_score()
            else:
                scorelist=[]
                self.calc_score_cluster()
                for pset in self.dslist:
                    nspo=spdsscore([pset],self.bm,self.scaledsp)
                    nspo.get_score()
                    scorelist.append(nspo.score)
                self.score=np.hstack(scorelist)
        return self.score
    
