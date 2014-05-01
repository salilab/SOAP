"""
   SOAP Model Selection module

"""
from env import *
from crossValidate import *

#import objgraph
#import multiprocessing.util as util
#util.log_to_stderr(util.SUBDEBUG)
debug=True
#import errno


def cvplot(cl):
    print 'plotting started'
    os.chdir(cl)
    if not os.path.isfile('figinput.pickle'):
        return 0
    cvcluster=mypickle().load('figinput')
    cvcluster.plot()
    print os.system('rm '+cl+'figinput*')
    print 'ploting finished'
    del cvcluster
    gc.collect()
    return 0

#from multiprocessing import Pool
#p = Pool(1,maxtasksperchild=1000)#,maxtasksperchild=1


class sps(object):
    def __init__(self,env=env(),modellist=[],evalPotFunc=None):
        self.env=env
        self.models=modellist
        self.modelso=copy.deepcopy(self.models)
        self.cvlist=[]
        self.tasktempdirs=[]
        self.logpathlist=[]
        self.evalPotFunc=evalPotFunc
        
    def append(self,m):
        self.models.append(copy.deepcopy(m))
        self.modelso.append(copy.deepcopy(m))
     
    def initialize(self):
        tl=[]
        for model in self.models:
            tcl=[]
            tcl.append(localtask(func=self.initialize_model_pir,inputlist=[model]))
            self.initialize_model_distributions(model,tcl)
            self.initialize_model_scores(model,tcl)
            tco=taskchain(chains=tcl)
            tl.append(tco)
        to=tasklist(tasklist=tl)
        
        to.monitor2end()
        
    def get_alltasks(self):
        tl=[]
        for model in self.models:
            tl.append(self.get_singlemodel_taskchain(model))
        return tasklist(tasklist=tl)
    
    def get_singlemodel_taskchain(self,model):
        tcl=[]
        tcl.append(localtask(func=self.initialize_model_pir,inputlist=[model]))
        self.initialize_model_distributions(model,tcl)
        self.initialize_model_scores(model,tcl)
        cv=k2cvcluster(model=model, initialize=False)
        self.cvlist.append(cv)
        tcl.append(cv.get_task())
        return taskchain(chains=tcl)


    def initialize_model_pir(self,model):
        pirlist=[]
        for scorer in model['scorers']:
            if scorer['type']=='scaledsp':
                pirlist.append(scorer['pdbset'])
        pirlist=squeeze_list(pirlist)
        for pird in pirlist:
            pir().get_pir(pird)
        return 1
        
    def initialize_model_distributions(self,model,taskchain):
        rawsplist=[]
        scaledsplist=[]
        for scorer in model['scorers']:
            if scorer['type']=='scaledsp':
                ns=copy.deepcopy(scorer)
                scaledsplist.append(scaledsp(model=ns))
                rawsplist.append(rawsp(model=ns))
            elif scorer['type']=='sf' and model['bmtype']['type']=='dsscore':
                if 'bm' in scorer:
                    bmt=scorer['bm']
                else:
                    bmt=model['bmtype']['bm']
                rawsplist.append(rawsp(pdbset='_'.join(model['bmtype']['dslist']),features=scorer['features'],genmethod=bmt,decoy=True,routine='calc_sp_i'))
        scaledsplist=squeeze_list(scaledsplist)
        rawsplist=squeeze_list(rawsplist)
        listoftasks=[]
        for rsp in rawsplist:
            r=rsp.get_task()
            if r!=0:
                listoftasks.append(r)
        if len(listoftasks)>0:
            taskchain.append(tasklist(tasklist=listoftasks))
        listoftasks=[]
        for ssp in scaledsplist:
            r=ssp.get_task()
            if r!=0:
                listoftasks.append(r)
        if len(listoftasks)>0:
            taskchain.append(tasklist(tasklist=listoftasks))
            
    def initialize_model_scores(self,model,taskchain):
        listoftasks=[]
        scorelist=[]
        searches=[]
        for item in model['searches']:
            if item['key']!='ratio':
                searches.append(item['object'])
        for scorer in model['scorers']:
            if scorer['type']=='scaledsp':
                if scorer in searches:
                    continue
                if model['bmtype']['type']=='dsscore':    
                    ns=copy.deepcopy(scorer)
                    ns['scoretype']='dsscore'
                    ns['dslist']=model['bmtype']['dslist']
                    if 'bm' in scorer:
                        bmt=scorer['bm']
                    else:
                        bmt=model['bmtype']['bm']
                    ns['bm']=bmt
                    scorelist.append(spdsscore(model=ns))
        scorelist=squeeze_list(scorelist)
        for score in scorelist:
            sj=score.get_task()
            if sj!=0:
                listoftasks.append(sj)
        if len(listoftasks)>0:
            taskchain.append(tasklist(tasklist=listoftasks))
            
    def cv_start(self,testperc,testsample):
        #self.loadlog()
        self.cvlist=[]
        listoftasks=[]
        for model in self.models:
            so=scorer(model=model)
            opt=optimizer(scorer=so)
            cv=k2cvcluster(spsfo=opt)
            self.cvlist.append(cv)
            listoftasks.append(cv.cross_validation_with_test(testperc,testsample))
        self.cvtask=tasklist(tasklist=listoftasks)
        
    def cv(self):
        ml=[]
        for mdl in self.models:
            nmr=modelinlog(mdl)
            if len(nmr)>0:
                print nmr[1]
                self.logpathlist.append(nmr)
            else:
                ml.append(mdl)
        self.models=ml
        alltasks=self.get_alltasks()
        #finished getting tasks, now start running
        #pdb.set_trace()
        runenv.spdsscoreTasks={}
        alltasks.monitor2end()
        self.show_summary()
        
    
    def cvlocal(self,testperc=0,testsample=[]):
        self.cvlist=[]
        for opt in self.optlist:
            self.cvlist.append(k2cvlocal(spsfo=opt))
        self.cvresult=[]
        for cv in self.cvlist:
            self.cvresult.append(cv.cross_validation_with_test(testperc,testsample))
        return self.cvresult

    def assess_basescore(self):
        for scorer in self.scorerlist:
            print scorer.assess_basescore()

    def write_optimal_potential(self,type='lib'):
        for i in range(0,len(self.scorerlist)):
            scorer=self.scorerlist[i]
            scorer.build_model(self.cvresult[i][-1][0][0][0])
            scorer.write_potential(type)
            
    def show_summary(self):
        for cv in self.cvlist:
            print cv.statsstr
        for cv in self.cvlist:
            print cv.bestmodelresult
        for cv in self.cvlist:
            print cv.resultstr
        
    def plot(self,cl=[]):
        clusterslist=[]
        for i in range(len(self.cvlist)):
            if (len(cl)>0 and (i in cl)) or len(cl)==0:
                cv=self.cvlist[i]
                #os.chdir(cv.clusters.figpath)
                clusterslist.append(cv.logpath)
        #p.map_async(cvplot, clusterslist)

    def show(self):
        plt.show()

    def loadlog(self):
        fh=open(runenv.basedir+'log.pickle','r')
        self.resultdict=pickle.load(fh)
        fh.close()

    def get_cvrundir(self):
        rdl=[]
        for cv in self.cvlist:
            rdl.append(cv.rundir)
        return rdl
  
class spss(object):
    def __init__(self, model=[],logfoldername='',evalPotFunc=None):
        self.dsearches=model['dsearches']
        self.env=env()
        self.model=model
        self.modelstr=any2str(model)
        self.complete_model()
        self.complete_rrf(model['searches'])
        self.dsearchpos=0
        self.originaldsearches=copy.deepcopy(self.dsearches)
        path=self.setup_logdir()
        self.runsnlist=[]
        self.resultdict={}
        self.evalPotFunc=evalPotFunc
        if len(path)>0:
            print path+'/spss.pickle'
            print 'skipping'
            fh=open(path+'/spss.pickle')
            a=pickle.load(fh)
            fh.close()
            for key in self.__dict__.keys():
                self.__dict__[key]=a.__dict__[key]
            self.logdir=a.logdir#'/bell2/gqdong/statpot/results/DecoysRUS_casp58_baker1/NativeSelection+0.1xcc/bndsep6--d30#15a158as158.d30#15-0.86009,-1308.715_0.915+0.004_0.861+0.005_0.906/'
            self.logdir2=a.logdir2
            return None
        self.write2log(self.get_valuesets())
        self.currentcutoffpercinitial=0.09
        self.currentcutoffperc=0.09
        self.currentcutoffpercratio=0.1
        self.maximumround=10
        self.numofevals=0
        self.num2evalall=10
        if 'refineProtocal' in model:
            self.refineProtocal=model['refineProtocal']
            del model['refineProtocal']
        else:
            self.refineProtocal={}
    
    def modelinlog(self,nm):
        os.chdir(self.baselogdir)
        with FileLock("log.shelve", timeout=100, delay=2) as lock:
            print("Lock acquired.")
            if not os.path.isfile(self.baselogdir+'log.shelve'):
                nmr=[]
            else:
                resultdictlog=shelve.open(self.baselogdir+'log.shelve')
                if nm['str'] in resultdictlog:
                    nmr=resultdictlog[nm['str']]
                    ko=k2cvcluster()
                    rundirname=[item for item in os.listdir(os.path.join(self.baselogdir,'runs')) if item.startswith(nmr[2])][0]
                    ko.load_fromlogdir(os.path.join(self.baselogdir,'runs',rundirname))
                    nmr[0]=ko.resultarray
                else:
                    nmr=[]
                resultdictlog.close()
            print("lock released")
        return nmr
    
    def modelexist(self):
        os.chdir(self.baselogdir)
        path=''
        with FileLock("spsslog.shelve", timeout=100, delay=2) as lock:
            print("Lock acquired.")
            if os.path.isfile(self.baselogdir+'spsslog.shelve'):
                resultdictlog=shelve.open(self.baselogdir+'spsslog.shelve')
                if self.modelstr in resultdictlog:
                    path=resultdictlog[self.modelstr]
                resultdictlog.close()
            print("lock released")
        return path

    def setup_logdir(self):
        bdir=runenv.basedir+'results/'
        tdir=bdir+'_'.join(self.model['bmtype']['dslist'])
        if not os.path.isdir(tdir):
            os.mkdir(tdir)
        tdir=tdir+'/'+self.model['bmtype']['criteria']
        if not os.path.isdir(tdir):
            os.mkdir(tdir)
        self.baselogdir=tdir+'/'
        path=''#self.modelexist()
        if len(path)>0:
            return path
        bdir2='/bell1/home/gqdong/Dropbox/spresults/'
        tdir=bdir2+'_'.join(self.model['bmtype']['dslist'])
        if not os.path.isdir(tdir):
            os.mkdir(tdir)
        tdir=tdir+'/'+self.model['bmtype']['criteria']
        if not os.path.isdir(tdir):
            os.mkdir(tdir)
        logfoldername=sys.argv[0][:-3]
        os.chdir(self.baselogdir)
        fl=os.listdir('./')
        existlogfoldernum=len([item for item in fl if item.startswith(logfoldername)])
        logfoldername=logfoldername+str(existlogfoldernum+1)
        self.logdir=self.baselogdir+logfoldername+'/'
        self.logdir2=tdir+'/'+logfoldername+'/'
        os.mkdir(self.logdir2)
        os.mkdir(self.logdir)
        os.chdir(self.logdir)
        return ''
        
    def complete_model(self):
        defaultmodel={'sm':['fmins','pbp.powell','powell'],'cvk':20,'cvm':'parallel','clustermethod':{'clusteringperc':0.999,'clustercutoff':2.5,'scoreratio':[1,1,0,3,1],'pickwhichone':'median'}}
        for key in defaultmodel:
            if not key in self.model:
                self.model[key]=defaultmodel[key]
            
    def complete_rrf(self,searches):
        return 0
        for so in searches:
            if so['key']=='parvalue':
                par=so['object']['par']
                so['object']['parvalue']=list(np.ones(len(par)))
                so['pos']=range(len(par))
            
    def get_model_list(self,dsearch):
        nml=[]
        mkl=[]
        mil=[]
        for i in range(len(dsearch['valueset'])):
            self.set_model(dsearch,i)
            modelkey=self.get_current_keyvalues()
            mkl.append(any2str(modelkey))
            if any2str(modelkey) in self.resultdict:
                continue
            nm=copy.deepcopy(self.model)
            convert2old(nm)
            nm['str']=any2str(nm)
            nmresult=self.modelinlog(nm)
            if len(nmresult)>0:
                self.resultdict[any2str(modelkey)]=nmresult
                continue
            nml.append(nm)
            mil.append(i)
        return [nml,mkl,mil]        
    
    def set_model(self,dsearch,i):
        sobj=dsearch['object']
        skey=dsearch['key']
        sval=dsearch['valueset']
        if isinstance(skey,list):
            for j in range(len(skey)):
                sobj[j][skey[j]]=sval[i][j]
        else:
            sobj[skey]=sval[i]
                
    def eval_single_parlist(self,dsearch):
        #get the list of models
        [nml,mkl,mil]=self.get_model_list(dsearch)
        #evaluate the list of models
        if len(nml)==0:
            print "no model to evaluate skip this round"
            return 0
        pickedindex=self.eval_modellist(nml,mkl,mil,'SinglePar\n')
        self.update_dsearches(pickedindex,dsearch)

    def eval_modellist(self,nml,mkl,mil,logprefix):
        spso=sps(modellist=nml,evalPotFunc=self.evalPotFunc)
        spso.cv()
        #spso.log()
        self.numofevals+=len(nml)
        for cv,i in zip(spso.cvlist,mil):
            self.resultdict[str(mkl[i])]=[cv.resultarray,cv.resultstr,cv.rundir,cv.bestmodel,cv.bestmodelresult] #the format is the same as the dictionary
        #hijack the result here to
        if self.evalPotFunc!=None:#use the eval function to further evaluate the devired potential
            rl=[]
            tcl=[]
            for key,value in self.resultdict.items():
                cvo=k2cvcluster(logpath=value[5])                
                cvo.generate_potential()
                rl.append([cvo,value[0]])
                tcl.append(self.evalPotFunc(cvo))#get refine task
            tasklist(tcl).monitor2end()
            for cvo, ra in rl:
                ra[-1]=cvo.resultsummary
        self.runsnlist=self.runsnlist+spso.get_cvrundir()
        pickedindex=self.analyze_results(mkl,logprefix)
        clusterslist=[]
        #for cv,i in zip(range(len(spso.cvlist)),mil):
        #    if i in list(pickedindex):
        #        clusterslist.append(cv)
        #spso.plot(clusterslist)
        del spso
        gc.collect()
        #if debug and objgraph.count('scorer')>0:
            #print 'some scorer object is still in memory'
            #pdb.set_trace()
        return pickedindex

    def show_current_result(self,spso,mkl):
        for cv,mk in zip(spso.cvlist,mkl):
            try:
                print cv.resultstr+' '+str(mk)
            except:
                pass
    
    def get_current_keyvalues(self):
        vl=[]
        for dsearch in self.dsearches:
            if not isinstance(dsearch['key'],list):
                vl.append(dsearch['object'][dsearch['key']])
            else:
                svl=[]
                for i in range(len(dsearch['key'])):
                    svl.append(dsearch['object'][i][dsearch['key'][i]])
                vl.append(svl)
        return vl
        
    def get_resultarray(self,mkl):
        ra=[]
        for mk in mkl:
            ra.append(self.resultdict[str(mk)][0])
        return np.vstack(ra)

    def get_logstr(self,mkl):
        loglist=[]
        for mk in mkl:
            loglist.append(str(self.resultdict[str(mk)][0][-1])+' '+self.resultdict[str(mk)][1]+' '+str(mk)+' '+self.resultdict[str(mk)][2])
        return '\n'.join(loglist)
            
    def pick_best_par(self,rarank,toprankbest,cti):
        if 'fix' in self.dsearches[self.dsearchpos] and self.dsearches[self.dsearchpos]['fix']:
            return rarank
        else:
            return rarank[0:cti]
        
    def rank_model(self,ra):
        ras=ra[:,-1]
        if ras.min()==ras.max():
            return [range(len(ras)),False,len(ras)]
        rassort=ras.argsort()
        rascutoff=ras.min()+(ras.max()-ras.min())*self.currentcutoffperc
        if ras.max()>0:
            rascutoff=min(rascutoff,ras.max()*self.currentcutoffperc)
        else:
            rascutoff=min(rascutoff,ras.max()*(2-self.currentcutoffperc))
        cti=np.searchsorted(ras[rassort],rascutoff,side='right')
        cti=len(ras)-cti
        cti=max(np.rint((1.0-self.currentcutoffperc)*len(ras)),cti)
        rassort=rassort[::-1]
        ramax=ra[:,[0,2,4]].max(axis=0)
        if (ra[rassort[0],[0,2,4]]==ramax).all():
            toprankbest=True
        else:
            toprankbest=False
        return [rassort,toprankbest,cti]
    
    def update_dsearches(self,pickedindex,dsearch):
        nvl=[]
        for i in pickedindex:
            nvl.append(dsearch['valueset'][i])
        try:
            a=np.random.rand() # randomize the selection of the index picked...
            if a<0.04*0.00000001:
                j=2
            elif a<0.125*0.00000001:
                j=1
            else:
                j=0
            j=min(j, len(pickedindex)-1)
            self.set_model(dsearch,pickedindex[j])
        except:
            print "in spss.update_dsearches"
            pdb.set_trace()
        dsearch['valueset']=nvl
        self.write2log(self.get_valuesets())
        
    def oneloop_through_dsearches(self):
        self.numofevals=0
        for dsearch,i in zip(self.dsearches,range(len(self.dsearches))):
            if len(dsearch)==1:
                continue
            self.dsearchpos=i
            self.eval_single_parlist(dsearch)
            
    def count_searchspace(self):
        n=1
        for dsearch in self.dsearches:
            n=n*len(dsearch['valueset'])
        return n

    def find_best_par(self):
        k=0
        while k<self.maximumround and self.count_searchspace()>self.num2evalall:
            self.currentcutoffperc=self.currentcutoffpercinitial+k*self.currentcutoffpercratio           
            self.oneloop_through_dsearches()
            k=k+1
            #if self.numofevals==0:
            #    break
        if self.count_searchspace()<=self.num2evalall:
            self.eval_allpars()
        mkl=self.resultdict.keys()
        self.analyze_results(mkl,'FinalAll\n')
        
    def analyze_results(self,mkl,logprefix=''):
        ra=self.get_resultarray(mkl)
        rarank,toprankbest,cti=self.rank_model(ra)
        nmkl=[mkl[i] for i in rarank]
        self.bestpar=nmkl[0]
        logstr=logprefix+self.get_logstr(nmkl)
        self.write2log(logstr)
        return self.pick_best_par(rarank,toprankbest,cti)

    def write2log(self,logstr):
        os.chdir(self.logdir)
        with FileLock("log.shelve", timeout=100, delay=2) as lock:
            print("Lock acquired.")
            logfile=open(self.logdir+'log','a+')
            logstr=logstr.encode('utf-8')
            logfile.write(logstr+'\n')
            print logstr
            logfile.close()
            print os.system('cp '+self.logdir+'log '+self.logdir2)
            print("lock released")        
        
    def get_valuesets(self):
        dsstr=''
        for dsearch in self.dsearches:
            dsstr=dsstr+str(dsearch['valueset'])
        return dsstr
        
    def eval_allpars(self):
        ml3=[[],[],[]]#model, key,index
        self.get_all_models(ml3,0)
        self.eval_modellist(ml3[0],ml3[1],ml3[2],'EvalAll\n')
        
    def log(self):
        #not used anymore
        os.chdir(self.logdir)
        self.write2log(self.get_valuesets())
        self.write2log('\n'+str(self.resultdict[str(self.bestpar)][4]))
        self.write2log('\n'+str(self.resultdict[str(self.bestpar)][3]))
        #self.dump()
        fl=os.listdir(self.baselogdir+'runs/')
        bbd=[f for f in fl if f.startswith(self.resultdict[str(self.bestpar)][2])][0]
        bbd='runs/'+bbd
        print os.system('cp '+self.baselogdir+bbd+'/bestmodel.pickle '+self.logdir)
        print os.system('touch '+self.logdir+'BMF-'+bbd)
        print os.system('cp '+sys.path[1]+'/'+sys.argv[0]+' '+self.logdir)
        print os.system('cp '+self.logdir+'log '+self.logdir2)# keep the log in dropbox
        print os.system('cp '+self.logdir+'*py '+self.logdir2)# keep the log in dropbox
        #print "Waiting for pool to finish"
        #p.close()
        #print "Pool close"
        #p.join()
        #print "Pool finished"
        print "The pool will just be ignored......(can't get it to work properly), now plotting the best one"
        cvplot(self.baselogdir+bbd)
        os.chdir(self.logdir)
        for key in self.resultdict:
            bd=[f for f in fl if f.startswith(self.resultdict[key][2])][0]
            print os.system('rm '+self.baselogdir+'runs/'+bd+'/figinput*')
        print os.system('cp '+self.baselogdir+bbd+'/*eps '+self.logdir2)
        print os.system('mv '+self.logdir+' '+self.logdir[:-1]+'-'+bd[5:])
        print os.system('mv '+self.logdir2+' '+self.logdir2[:-1]+'-'+bd[5:])
        self.logdir=self.logdir[:-1]+'-'+bd[5:]+'/'
        self.logdir2=self.logdir2[:-1]+'-'+bd[5:]+'/'
        with FileLock("spsslog.shelve", timeout=100, delay=2) as lock:
            print("Lock acquired.")
            resultdictlog=shelve.open(self.baselogdir+'spsslog.shelve')
            resultdictlog[self.modelstr]=self.logdir[:-1]+'-'+bd[5:]
            resultdictlog.close()
            print("lock released")
        return bbd
        
    def write_best_pot(self):
        print "only works for the case with single pot, single reference state"
        fh=open(self.logdir+'bestmodel.pickle')
        bm=pickle.load(fh)
        fh.close()
        #pdb.set_trace()
        #bm['scorers'][1]['sftype']='logbins'
        bm['scorers'][0]['refs']=[bm['scorers'][1]]
        so=scorer(model=bm)
        mre=so.assess_basescore()
        ssppath=so.write_potential(type=type)
        print mre
        print ssppath+'.opt.'+type
        
    def dump(self):
        os.chdir(self.logdir)
        fh=open('spss.pickle','w')
        pickle.dump(self,fh)
        fh.close()
        
    def get_all_models(self,ml3, dsearchindex):
        dsearch=self.dsearches[dsearchindex]
        if dsearchindex==(len(self.dsearches)-1):
            [ml,mkl,mil]=self.get_model_list(dsearch)
            ml3[2]=ml3[2]+list(np.array(mil)+len(ml3[1]))
            ml3[0]=ml3[0]+ml
            ml3[1]=ml3[1]+mkl
        else:
            for i in range(len(dsearch['valueset'])):
                self.set_model(dsearch,i)
                self.get_all_models(ml3,dsearchindex+1)

    def show_finished_runresults(self):
        ml3=[[],[],[]]#model, key,index
        self.get_all_models(ml3,0)
        for key in self.resultdict:
            print self.resultdict[key][1]+', '+str(key)
    
def pickle2shelf(prefix):
    fh=open(prefix+'.pickle')
    a=pickle.load(fh)
    fh.close()
    print os.system('rm '+prefix+'.pickle')
    nd=shelve.open(prefix+'.shelve')
    for key in a:
        nd[key]=a[key]
    nd.close()
        
def get_ppdock_plotdata(path2logdir,scoretype='model',bm=[],mindex=-1,bmset='',scorers=[],sname='',accuracy='rmsd10ORirmsd4FIRST'):
    if bm==[]:
        bm=get_best_model(path2logdir)
    #del bm['searches']
    if bmset:
        bm['bmtype']['dslist']=[bmset]
    if scorers:
        bm['scorers']=scorers
        scoretype='scorer'
    #cvmodel=pickle.load(open(os.path.join(path2logdir,'cvmodel.pickle')))
    #bestpar=cvmodel['originalresult'][0][0]['bestpar']
    #cvo=k2cvcluster(model=bm,initialize=True)
    #cvo.prepare_cross_validation_sample()
    #so=cvo.scorerlist[-1][mindex]
    #pdb.set_trace()
    so=scorer(model=bm)
    #print so.assess_model(slevel='top1000_rlrank__'+accuracy)
    print so.assess_model(slevel='top1000_rlrank__'+accuracy)
    #print so.assess_model()
    print so.assess_basescore()
    tnr=range(0,30)
    ra=np.zeros(30)
    idealperc=so.assess_ideal(slevel='top1000__nonempty_'+accuracy)
    for p in tnr:
        tcp=int(10**(p/10.0))
        if scoretype=='sascore':
            ra[p]=so.assess_sascore(slevel='top'+str(tcp)+'__nonempty_'+accuracy)
        elif scoretype.startswith('model'):
            ra[p]=so.assess_model(slevel='top'+str(tcp)+'__nonempty_'+accuracy)
        elif scoretype=='random':
            ra[p]=so.assess_randomized_score(slevel='top'+str(tcp)+'__nonempty_'+accuracy)
        elif scoretype=='scorer':
            ra[p]=so.assess_basescore(slevel='top'+str(tcp)+'__nonempty_'+accuracy)            
    ra[...]=ra[...]/idealperc
    if sname:
        np.save(os.path.join(path2logdir,bmset+sname+'.npy'),ra)
    elif scoretype=='sascore':
        np.save(os.path.join(path2logdir,bmset+'sascore.npy'),ra)
    elif scoretype.startswith('model'):
        np.save(os.path.join(path2logdir,bmset+scoretype+'.npy'),ra)
    elif scoretype=='random':
        np.save(os.path.join(path2logdir,bmset+'random.npy'),ra)
    elif scoretype=='scorer':
        np.save(os.path.join(path2logdir,bmset+'scorer.npy'),ra)        
    return ra

def get_bestm_result(path2logdir,scoretype='model',bm=[],mindex=-1,dslist=''
                     ,scorerlist=[],sname='',ctl=[], ct='', testperc=0,alldecoys=False):
    if not bm:
        bm=get_best_model(path2logdir)
    #del bm['searches']
    if ct:
        bm['bmtype']['criteria']=ct
    if dslist:
        bm['bmtype']['dslist']=dslist
    if scorerlist:
        bm['scorers']=scorerlist
        scoretype='scorer'
        bm['searches']=[]
    if len(ctl)==0:
        ctl=[bm['bmtype']['criteria']]
    if testperc!=0:
        bm['testperc']=testperc
        bm['cvk']=0
        spso=sps(modellist=[bm])
        spso.cv()
        print spso.cvlist[0].logpath
        return spso.cvlist[0].logpath
    #cvo=k2cvcluster(model=bm,initialize=True)
    #cvo.prepare_cross_validation_sample()
    #if alldecoys:
    #    so=cvo.spsfo.scorer
    #else:
    #    so=cvo.scorerlist[-1][mindex]
    bm={'scorers':bm['scorers'],'bmtype':bm['bmtype'],'searches':bm['searches']}
    try:
        so=scorer(model=bm)
    except:
        so=scorer(model=bm)
    #idealperc=so.assess_ideal(slevel='top100_rmsd10_irmsd4')
    print so.assess_model()
    ra=np.zeros(len(ctl))
    for p,pi in zip(ctl,range(len(ctl))):
        if scoretype=='sascore':
            ra[pi]=so.assess_sascore(slevel=p)
        elif scoretype.startswith('model'):
            ra[pi]=so.assess_model(slevel=p)
        elif scoretype=='random':
            ra[pi]=so.assess_randomized_score(slevel=p)
        elif scoretype=='scorer':
            ra[pi]=so.assess_basescore(slevel=p)
        elif scoretype=='ideal':
            ra[pi]=so.assess_ideal(slevel=p)
    
    if sname:
        np.save(os.path.join(path2logdir,bmset+sname+'.npy'),ra)       
    return ra

def get_best_model(path2logdir):
    print "only works for the case with single pot, single reference state"
    #fh=open(os.path.join(path2logdir,'bestmodel.pickle'))
    #bm=pickle.load(fh)
    #fh.close()
    fh=open(os.path.join(path2logdir,'cvmodel.pickle'))
    bm=pickle.load(fh)
    fh.close()
    so=scorer(model=bm['model'])
    bestpars=bm['originalresult'][0][0]['bestpar']
    print so.assess_model()
    return so.model

def plot_ppdock_result(figpath,fl,lt):
    plt.clf()
    xb=range(0,30)
    x=[int(10**(p/10.0)) for p in xb]
    for f in fl:
        na=np.load(f)
        plt.semilogx(x,na)
    plt.xlabel('Number of predictions')
    plt.ylabel('Success rate')
    plt.legend(lt,loc=0)
    plt.savefig(figpath)
    
def plot_reference_state(figpath,reflist,lt):
    xb=range(0,30)
    x=[int(10**(p/10.0)) for p in xb]
    for f in reflist:
        sfo=sf(model=f)
        csf=sfo.get_sf()
        plt.plot(csf)
    plt.xlabel('Number of predictions')
    plt.ylabel('Success rate')
    plt.legend(lt,loc=0)
    plt.savefig(figpath)
  
def generate_uniformbins(pd):
    fbs=pd['featurebins']#0,6,0.5
    bs=fbs[1]
    be=fbs[2]
    bv=fbs[3]
    bins=list(np.arange(bs+(bv/2.0),be,bv))
    n=pd['uniform']+1
    if (n+1)>=len(bins):
        return bins[1:-1]
    par=[]
    if 'firstpoint' in pd:
        fps=pd['firstpoint']
        il=range(len(bins))
        il.reverse()
        for i in il:
            if bins[i]<pd['firstpoint']:
                bins.pop(i)
        par.append(pd['firstpoint'])
        n-=1
    if (n+1)>=len(bins):
        return bins[:-1]
    t=len(bins)-1
    b=t/n
    x=(b+1)*n-t
    inds=[b*(i+1) for i in range(x)]+[(b+1)*(i+1)+b*x for i in range(n-x)]
    if pd['distribute']=='lownumberpriority':
        par=par+[bins[i] for i in inds][:-1]
    else:
        print traceback.print_exec()
        raise Exception('not implemented')
    return par

class parlist(list):
    pass


def preprocessing(nm):
    for sco in nm['scorers']:
        for key in sco:
            try:
                if isinstance(sco[key],list) and len(sco[key])>0  and sco[key][0]=='container':
                    sco[key]=sco[key][1]
            except:
                pdb.set_trace()


def convert2old(nm):
    if 'dsearches' in nm:
        del nm['dsearches']
    if isinstance(nm['scorers'],dict):
        ssd=nm['scorers']
        if 'scorers' in ssd:
            nm['scorers']=ssd['scorers']
            nm['searches']=ssd['searches']
        elif 'index' in ssd:
            nm['scorers']=ssd['valueset'][ssd['index']][0]
            nm['searches']=ssd['valueset'][ssd['index']][1]
    k=-1
    preprocessing(nm) 
    for search in nm['searches']:
        #k+=1
        #search=nm['searches'][k]
        scorer=search['object']
        if len(scorer[search['key']])==0:
            nm['searches'].pop(k)
            continue
        print scorer
        #if scorer['parvalue']==None:
        #    pdb.set_trace()
        #print search['key']=='parvalue' and search['object']['sftype']=='bins'
        if (search['key']=='parvalue') and ('parvalue' in scorer) and isinstance(scorer['parvalue'],dict):
            npl=[]
            sopar=scorer['parvalue']
            pardict=sopar['valueset'][sopar['key']]
            scorer['parvalue']=copy.deepcopy(pardict['parvalue'])
            scorer['sftype']=copy.deepcopy(pardict['sftype'])
            scorer['par']=copy.deepcopy(pardict['par'])
            if 'rawsp' in pardict:
                scorer['rawsp']=pardict['rawsp']
        elif (search['key']=='parvalue') and ('par' in scorer) and isinstance(scorer['par'],dict):
            if 'uniform' in scorer['par']:
                par=generate_uniformbins(scorer['par'])
                scorer['par']=par
                scorer['parvalue']=[1 for i in range(len(par)+2)]
                search['pos']=range(0,len(scorer['parvalue']))
                search['pos']=range(0,len(scorer['par']))
            else:
                npl=[]
                sopar=scorer['par']
                parrange=sopar['object'][sopar['index1']][sopar['index2']]
                partype=sopar['typelist'][sopar['parindex']]
                if partype[0]=='fixbin':
                    par=range(partype[1],parrange,partype[2])
                    par.append(parrange)
                elif partype[0]=='numafter':
                    par=partype[1]+list(np.linspace(partype[1][-1],parrange,partype[2]+1))[1:]
            scorer['par']=par
        if search['key']=='parvalue' and search['object']['sftype']=='bins' and 'd' in search['object']['features']:
            #pdb.set_trace()
            fbs=search['object']['features'][-1]
            bs=fbs[1]
            be=fbs[2]
            bv=fbs[3]
            bins=list(np.arange(bs+(bv/2.0),be,bv))
            search['object']['parvalue']=bins
        if 0 and (search['key']=='parvalue') and ('par' in scorer) and len(scorer['par'])>0:
            scorer['parvalue']=[1 for item in scorer['par']]
        try:
            if (search['key']=='parvalue'):
                search['pos']=range(0,len(scorer['parvalue']))
            elif search['key']=='par':
                search['pos']=range(0,len(scorer['par']))
        except Exception,e:
            traceback.print_exc()
            print scorer
            print search
            pdb.set_trace()
    if ('bm' in nm['bmtype']) and isinstance(nm['bmtype']['bm'],list):
        npl=[]
        for item in nm['bmtype']['bm']:
            if isinstance(item,list):
                npl.append(''.join(item))                      
            else:
                npl.append(item)
        nnpl=[item for item in npl if len(item)>0]
        nm['bmtype']['bm']=''.join(nnpl)    
    for scorer in nm['scorers']:
        if ('features' in scorer) and isinstance(scorer['features'],list):
            nfl=[]
            for item in scorer['features']:
                if isinstance(item,list):
                    try:
                        featurerange=[item[1],item[2]]
                        if item[1]==0:
                            nfl.append(item[0]+str(int(0.5+(item[2]-item[1])/item[3]))+'#'+str(item[2]))
                        else:
                            nfl.append(item[0]+str(int(0.5+(item[2]-item[1])/item[3]))+'#'+str(item[2])+'#'+str(item[1]))                        
                    except:
                        pdb.set_trace()
                else:
                    nfl.append(item)
            scorer['features']=''.join(nfl)
        if ('pm' in scorer) and isinstance(scorer['pm'],list):
            npl=[]
            for item in scorer['pm']:
                if isinstance(item,list):
                    npl.append(''.join(item))                      
                else:
                    npl.append(item)
            nnpl=[item for item in npl if len(item)>0]
            scorer['pm']='+'.join(nnpl)
        for pn in ['bm','genmethod','pdbset']:
            if (pn in scorer) and isinstance(scorer[pn],list):
                npl=[]
                for item in scorer[pn]:
                    if isinstance(item,list):
                        npl.append(''.join(item))                      
                    else:
                        npl.append(item)
                nnpl=[item for item in npl if len(item)>0]
                scorer[pn]=''.join(nnpl)

        if ('sftype' in scorer) and isinstance(scorer['sftype'],dict):
            ssd=scorer['sftype']
            scorer['par']=copy.deepcopy(ssd['valueset'][ssd['key']]['par'])
            scorer['parvalue']=copy.deepcopy(ssd['valueset'][ssd['key']]['parvalue'])
            scorer['sftype']=copy.deepcopy(ssd['valueset'][ssd['key']]['sftype'])
            if 'rawsp' in ssd['valueset'][ssd['key']]:
                scorer['rawsp']=ssd['valueset'][ssd['key']]['rawsp']
        if ('parvalue' in scorer) and isinstance(scorer['parvalue'],dict):
            try:
                scorer['par']=copy.deepcopy(scorer['par']['par'])
                scorer['parvalue']=copy.deepcopy(scorer['parvalue']['parvalue'])
            except:
                pdb.set_trace()        

        if ('par' in scorer) and isinstance(scorer['par'],dict):
            if 'uniform' in scorer['par']:
                par=generate_uniformbins(scorer['par'])
                scorer['par']=par
                scorer['parvalue']=[1 for i in range(len(par)+2)]
            else:
                npl=[]
                sopar=scorer['par']
                parrange=sopar['object'][sopar['index1']][sopar['index2']]
                partype=sopar['typelist'][sopar['parindex']]
                if partype[0]=='fixbin':
                    par=range(partype[1],parrange,partype[2])
                    par.append(parrange)
                elif partype[0]=='numafter':
                    par=partype[1]+list(np.linspace(partype[1][-1],parrange,partype[2]+1))[1:]
                scorer['par']=par
                scorer['parvalue']=[1 for item in par]
        for search in nm['searches']:
            if scorer==search['object']:
                if search['key']=='parvalue':
                    search['pos']=range(0,len(scorer['parvalue']))
                elif search['key']=='par':
                    search['pos']=range(0,len(scorer['par']))
                    #search['object']['parvalue']=list(np.ones(len(par)))
                elif isinstance(search['key'],list) and 'parvalue' in search['key']:
                    ind=search['key'].index('parvalue')
                    search['pos'][ind]=range(0,len(scorer['parvalue']))
    print nm
    return nm

def modelinlog(nm):
    bdir=runenv.basedir+'results/'
    tdir=bdir+'_'.join(nm['bmtype']['dslist'])
    tdir=tdir+'/'+nm['bmtype']['criteria']
    logdir=tdir+'/'
    os.chdir(logdir)
    with FileLock("log.shelve", timeout=100, delay=2) as lock:
        print("Lock acquired.")
        if not os.path.isfile(logdir+'log.shelve'):
            nmr=[]
        else:
            resultdictlog=shelve.open(logdir+'log.shelve')
            if not 'str' in nm:
                nm['str']=any2str(nm)
            if nm['str'] in resultdictlog:
                nmr=resultdictlog[nm['str']]
            else:
                nmr=[]
            resultdictlog.close()
        print("lock released")
    return nmr



        
        
    