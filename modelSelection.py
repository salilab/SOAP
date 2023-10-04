"""
   SOAP Model Selection module

"""
from __future__ import print_function
from env import *
from crossValidate import *
import filelock
import shutil
import glob
import os
import subprocess

debug=True


def cvplot(cl):
    print('plotting started')
    os.chdir(cl)
    if not os.path.isfile('figinput.pickle'):
        return 0
    cvcluster=mypickle().load('figinput')
    cvcluster.plot()
    for g in glob.glob(os.path.join(cl, 'figinput*')):
        os.unlink(g)
    print('ploting finished')
    del cvcluster
    gc.collect()
    return 0



class sps(object):
    """
    Evaluate a list of models on sge cluster
    """
    def __init__(self,modellist=[],evalPotFunc=None):
        self.env=runenv
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
        if 'refineProtocal' in model:
            nm=copy.deepcopy(model)
            del nm['refineProtocal']
        else:
            nm=model
        cv=k2cvcluster(model=nm, initialize=False)
        self.cvlist.append(cv)
        t=cv.get_task()
        if t!=0:
            tcl.append(t)
        if 'refineProtocal' in model:#use the eval function to further evaluate the devired potential
            lfo=sprefine(refineProtocal=model['refineProtocal'])
            rt=lfo.get_task(cv)
            if rt!=0:
                tcl.append(rt)
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
            print(scorer.assess_basescore())

    def write_optimal_potential(self,type='lib'):
        for i in range(0,len(self.scorerlist)):
            scorer=self.scorerlist[i]
            scorer.build_model(self.cvresult[i][-1][0][0][0])
            scorer.write_potential(type)

    def show_summary(self):
        for cv in self.cvlist:
            print(cv.statsstr)
        for cv in self.cvlist:
            print(cv.bestmodelresult)
        for cv in self.cvlist:
            print(cv.resultstr)

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
        with open(runenv.basedir+'log.pickle','rb') as fh:
            self.resultdict=pickle.load(fh)

    def get_cvrundir(self):
        rdl=[]
        for cv in self.cvlist:
            rdl.append(cv.rundir)
        return rdl

class spss(object):
    """

    Select the best model.

    :Parameters:
      - `model`: the model space to search the best model within

    Common usage::
        spl=spss(model=model1)
        spl.find_best_par()
        spl.log()

    .. note::
        The model and model space are defined using python lists and dictionaries.

        Define the recovery function features, distance only, from 0 to 20A with bin size 0.05, see the :mod:`features`. ::

           rfeature=[['d',0,20,0.05]] # 'd' - distance, 0 - start position, 20 - end position, 0.05 - bin size

        Features can be defined using either list or string, list defintion will be converted into string defintion on the fly.

        Define the spliens used for generating recoery functions. see :mod:`recoveryFunction`::

            par={'uniform':4,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
            slo={'key':'sacnfflex',
            'valueset':{
                            'sacnfflex':{'sftype':'sacnf','par':par,'parvalue':[1,1,1,1]},
                            'sacnf52':{'sftype':'sacnf','par':[2.75,3.75],'parvalue':[1,1,1,1]},
            }}

        The most important parameter to vary is 'uniform', which defines the number of anchor points to use.

        Define the recovery function, see :mod:`recoveryFunction`::

            ref1={'type':'sf','features':rfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0]}

        Define the features for the statistics calculation from a set of structures, a158 represent residue dependent atom type::

            sfeature=[rfeature[0],'a158','as158']

        Define the processing method of the probabilistic table, 'npend' means normalized by the last bin value, see :class:`statsTable.scaledsp`::

            pm=['','npend']

        Processing method can be defined using either list or string, list defintion will be converted into string definition on the fly.


        Define the probabilistic table using for scoring::

            scaledsp1={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree','features':sfeature,'genmethod':'cs1','pm':pm,'refs':[ref1],'ratio':[1.0]}


        The benchmark criteria, defines how the statistical potentials are benchmarked, see :mod:`benchmark`::

            bmtype={'type':'dsscore','dslist':['fpdb','fpdu'],'criteria':'3xtop1_rmsdallif_mean_+2xtop3_rmsdallif_mean_','combine':'scoresum','bm':'cs1'}

        * 'type' can be 'dsscore' or 'refinescore' (use SOAP for refine).
        * 'combine' can only be 'scoresum' at the moment.
        * 'bm' has the same meaning as :ref:`genmethod <genmethod>`.

        Parameters to optimize::

            search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}
            search2={'object':ref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}

        * 'Object': the object to optimize, needs to be a dict
        * 'key': the key of the object to optimize, must be a list
        * 'pos': the positions of the values in the list to optimize.
        * 'InitialGenerator': ways to generate initial value

        Discrete search dictionaries, defining the model space - the model variables/options to vary::

            dsearch4={'object':par,'key':'uniform','valueset':[7,6,5,4,3]}

            dsearch9={'object':[scaledsp1,scaledsp1],'key':['genmethod','pdbset'],'valueset':[['cs1','X2_2.2A_0.25rfree']]}#,['cs1','X2_2.2A_0.25rfree_30'],['cs1','X2_2.2A_0.25rfree_60'],['cs1','X2_2.2A_0.25rfree_95'],['bs20dsp','X_2.2A_0.25rfree'],['bs20dsp','X_2.2A_0.25rfree_30'],['bs20dsp','X_2.2A_0.25rfree_60'],['bs15dsp','X_2.2A_0.25rfree'],['bs10dsp','X_2.2A_0.25rfree']]}#,,'cs1'

        * 'Object': the object to vary, needs to be a dict
        * 'key': the key of the object to vary, often a string
        * 'valueset': the set of values to search within.

        Parameters controling sampling and optimization, please check the code in the sampling module for detailed meanings::

            ni=40

            initvalues=list(np.arange(0,ni)/10.0+1)

            initvalues=np.array(initvalues)

            inner={'improved':2,'maxloop':100,'minloop':2}
            outer={'improved':4,'maxloop':5023}

            td=list(np.sqrt(np.linspace(1,float(10)**2,ni)))
            td=np.array(td)
            tune_interval=200

            sampleschedule={'inner':inner,'outer':outer}
            ssm={'sm':'mcp','reset_betweenruns':2,'blockupdate':True, 'exchange_method':1,
                  'sample_schedule':sampleschedule,
                  'stepmethod':'mxmp2.0','tune_interval':200,'add_bias':False,'temperature_distribution':td}

            ssm0={'sm':'mcs','reset_betweenruns':2,'blockupdate':False,'using_globalbest':True,
                  'sample_schedule':sampleschedule,
                  'stepmethod':'mxmp2.0','tune_interval':201,'add_bias':False,'temperature_distribution':td}


            ssm2={'sm':'mcp','reset_betweenruns':2,'blockupdate':False, 'exchange_method':1,
                  'sample_schedule':sampleschedule,
                  'stepmethod':'mxmp2.0','tune_interval':200,'add_bias':False,'temperature_distribution':td}

            ssm20={'sm':'mcs','reset_betweenruns':2,'blockupdate':True,'using_globalbest':True,
                  'sample_schedule':sampleschedule,
                  'stepmethod':'mxmp2.0','tune_interval':201,'add_bias':False,'temperature_distribution':td}


            ssm1={'sm':'mca','reset_betweenruns':2,'blockupdate':True, 'using_globalbest':True,
                  'sample_schedule':sampleschedule,
                  'stepmethod':'mxmp2.0','tune_interval':201,'add_bias':False,'temperature_distribution':np.zeros(ni)+2}

            ssm3={'sm':'powell','blockupdate':False}

            sml=[ssm20,ssm2,ssm0,ssm,ssm0,ssm,ssm0,ssm, ssm1]


        Define the final model::

            model1={'scorers':[scaledsp1,ref1],'bmtype':bmtype,'searches':[search1,search2], 'runsperscorer':ni,
                'dsearches':[dsearch2,dsearch5,dsearch4,dsearch7,dsearch8,dsearch9],'sml':sml,'cvk':2,'repeat':1,'fold':3}

        * 'scorers': scoring terms
        * 'bmtype' : benchmark criteria for juding a statistical potential
        * 'searches' : parameters to optimzie
        * 'dsearches' : model options/values for model seleciton
        * 'sml' : optimization method
        * 'runsperscorer' : how many runs per replica exchange
        * 'cvk' : k times cross validation
        * 'repeat' : how many replica exchanges to carry out for each optimization run
        * 'fold' : n fold cross validation
        * 'testperc': if included, this percentage of decoys will be left out for final validation.


    """
    def __init__(self, model=[],logfoldername='',evalPotFunc=None):
        self.dsearches=model['dsearches']
        self.env=runenv
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
            print(path+'/spss.pickle')
            print('skipping')
            with open(path+'/spss.pickle', 'rb') as fh:
                a=pickle.load(fh)
            for key in self.__dict__.keys():
                self.__dict__[key]=a.__dict__[key]
            self.logdir=a.logdir
            #self.logdir2=a.logdir2
            return None
        self.write2log(self.get_valuesets())
        self.currentcutoffpercinitial=0.09
        self.currentcutoffperc=0.09
        self.currentcutoffpercratio=0.1
        self.maximumround=10
        self.numofevals=0
        self.num2evalall=20
        if 'refineProtocal' in model:
            self.refineProtocal=model['refineProtocal']
            #del model['refineProtocal']
        else:
            self.refineProtocal={}
        self.currentOptimal=-np.inf
        self.bestpar=[]

    def modelinlog(self,nm):
        os.chdir(self.baselogdir)
        with filelock.FileLock("log.shelve", timeout=100, delay=2) as lock:
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
        with filelock.FileLock("spsslog.shelve", timeout=100, delay=2) as lock:
            print("Lock acquired.")
            if os.path.isfile(self.baselogdir+'spsslog.shelve'):
                resultdictlog=shelve.open(self.baselogdir+'spsslog.shelve')
                if self.modelstr in resultdictlog:
                    path=resultdictlog[self.modelstr]
                resultdictlog.close()
            print("lock released")
        return path

    def setup_logdir(self):
        bdir=os.path.join(runenv.basedir, 'results')
        tdir=bdir+'_'.join(self.model['bmtype']['dslist'])
        if not os.path.isdir(tdir):
            os.mkdir(tdir)
        tdir=os.path.join(tdir, self.model['bmtype']['criteria'])
        if not os.path.isdir(tdir):
            os.mkdir(tdir)
        self.baselogdir=tdir+'/'
        path=''#self.modelexist()
        if len(path)>0:
            return path
        #bdir2='~/Dropbox/spresults/'
        #tdir=bdir2+'_'.join(self.model['bmtype']['dslist'])
        #if not os.path.isdir(tdir):
        #    os.mkdir(tdir)
        #tdir=tdir+'/'+self.model['bmtype']['criteria']
        #if not os.path.isdir(tdir):
        #    os.mkdir(tdir)
        logfoldername=sys.argv[0][:-3]
        os.chdir(self.baselogdir)
        fl=os.listdir('./')
        existlogfoldernum=len([item for item in fl if item.startswith(logfoldername)])
        logfoldername=logfoldername+str(existlogfoldernum+1)
        self.logdir=self.baselogdir+logfoldername+'/'
        #self.logdir2=tdir+'/'+logfoldername+'/'
        #os.mkdir(self.logdir2)
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
            #nmresult=self.modelinlog(nm)
            #if len(nmresult)>0:
            #    self.resultdict[any2str(modelkey)]=nmresult
            #    continue
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
        #if len(nml)==0:
        #    print("no model to evaluate skip this round")
        #    return 0
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
            #print('some scorer object is still in memory')
            #pdb.set_trace()
        return pickedindex


    def show_current_result(self,spso,mkl):
        for cv,mk in zip(spso.cvlist,mkl):
            try:
                print(cv.resultstr+' '+str(mk))
            except:
                pass

    def get_current_keyvalues(self):
        vl=[]
        for dsearch in self.dsearches:
            try:
                if not isinstance(dsearch['key'],list):
                    vl.append(dsearch['object'][dsearch['key']])
                else:
                    svl=[]
                    for i in range(len(dsearch['key'])):
                        svl.append(dsearch['object'][i][dsearch['key'][i]])
                    vl.append(svl)
            except:
                pdb.set_trace()
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
        if ras.max()>=self.currentOptimal:
            self.currentOptimal=ras.max()
        else:
            pdb.set_trace()
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
            print("in spss.update_dsearches")
            pdb.set_trace()
        dsearch['valueset']=nvl
        self.write2log(self.get_valuesets())

    def oneloop_through_dsearches(self):
        self.numofevals=0
        for i,dsearch in enumerate(self.dsearches):
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
            oldbestpar=self.bestpar
            self.currentcutoffperc=self.currentcutoffpercinitial+k*self.currentcutoffpercratio
            self.oneloop_through_dsearches()
            k=k+1
            if self.currentcutoffpercratio==0 and self.bestpar==oldbestpar:
                break
            #if self.numofevals==0:
            #    break
            if self.count_searchspace()<=self.num2evalall:
                self.eval_allpars()
                break
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
        with filelock.FileLock("log.shelve", timeout=100, delay=2) as lock:
            print("Lock acquired.")
            with open(self.logdir+'log','a+') as logfile:
                logstr=logstr.encode('utf-8')
                logfile.write(logstr+'\n')
                print(logstr)
            #print(os.system('cp '+self.logdir+'log '+self.logdir2))
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
        try:
            os.chdir(self.logdir)
            self.write2log(self.get_valuesets())
            self.write2log('\n'+str(self.resultdict[str(self.bestpar)][4]))
            self.write2log('\n'+str(self.resultdict[str(self.bestpar)][3]))
            #self.dump()
            fl=os.listdir(os.path.join(self.baselogdir, 'runs'))
            bbd=[f for f in fl if f.startswith(self.resultdict[str(self.bestpar)][2])][0]
            bbd = os.path.join('runs', bbd)
            shutil.copy(os.path.join(self.baselogdir, bbd, 'bestmodel.pickle'),
                        self.logdir)
            open(os.path.join(self.logdir, 'BMF-'+bbd), 'w') # touch
            shutil.copy(os.path.join(sys.path[1], sys.argv[0]), self.logdir)
            #print(os.system('cp '+self.logdir+'log '+self.logdir2)# keep the log in dropbox)
            #print(os.system('cp '+self.logdir+'*py '+self.logdir2)# keep the log in dropbox)
            #print("Waiting for pool to finish")
            #p.close()
            #print("Pool close")
            #p.join()
            #print("Pool finished")
            print("The pool will just be ignored......(can't get it to work properly), now plotting the best one")
            cvplot(self.baselogdir+bbd)
            os.chdir(self.logdir)
            for key in self.resultdict:
                bd=[f for f in fl if f.startswith(self.resultdict[key][2])][0]
                for g in glob.glob(os.path.join(self.baselogdir, 'runs',
                                                bd, 'figinput*')):
                    os.unlink(g)
            #print(os.system('cp '+self.baselogdir+bbd+'/*eps '+self.logdir2))
            shutil.move(self.logdir, self.logdir[:-1]+'-'+bd[5:])
            #print(os.system('mv '+self.logdir2+' '+self.logdir2[:-1]+'-'+bd[5:]))
            self.logdir=self.logdir[:-1]+'-'+bd[5:]+'/'
            #self.logdir2=self.logdir2[:-1]+'-'+bd[5:]+'/'
            with filelock.FileLock("spsslog.shelve", timeout=100,
                                   delay=2) as lock:
                print("Lock acquired.")
                resultdictlog=shelve.open(os.path.join(self.baselogdir,
                                                       'spsslog.shelve'))
                resultdictlog[self.modelstr]=self.logdir[:-1]+'-'+bd[5:]
                resultdictlog.close()
                print("lock released")
            return bbd
        except:
            pass

    def write_best_pot(self):
        print("only works for the case with single pot, single reference state")
        with open(self.logdir+'bestmodel.pickle', 'rb') as fh:
            bm=pickle.load(fh)
        #pdb.set_trace()
        #bm['scorers'][1]['sftype']='logbins'
        bm['scorers'][0]['refs']=[bm['scorers'][1]]
        so=scorer(model=bm)
        mre=so.assess_basescore()
        ssppath=so.write_potential(type=type)
        print(mre)
        print(ssppath+'.opt.'+type)

    def dump(self):
        os.chdir(self.logdir)
        with open('spss.pickle','wb') as fh:
            pickle.dump(self,fh)

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
            print(self.resultdict[key][1]+', '+str(key))

def pickle2shelf(prefix):
    with open(prefix+'.pickle', 'rb') as fh:
        a=pickle.load(fh)
    os.unlink(prefix+'.pickle')
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
    #print(so.assess_model(slevel='top1000_rlrank__'+accuracy))
    print(so.assess_model(slevel='top1000_rlrank__'+accuracy))
    #print(so.assess_model())
    print(so.assess_basescore())
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
        print(spso.cvlist[0].logpath)
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
    print(so.assess_model())
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
    print("only works for the case with single pot, single reference state")
    #fh=open(os.path.join(path2logdir,'bestmodel.pickle'))
    #bm=pickle.load(fh)
    #fh.close()
    with open(os.path.join(path2logdir,'cvmodel.pickle'), 'rb') as fh:
        bm=pickle.load(fh)
    so=scorer(model=bm['model'])
    bestpars=bm['originalresult'][0][0]['bestpar']
    print(so.assess_model())
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
        print(traceback.print_exec())
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
        print(scorer)
        #if scorer['parvalue']==None:
        #    pdb.set_trace()
        #print(search['key']=='parvalue' and search['object']['sftype']=='bins')
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
            if search['key']=='parvalue':
                search['pos']=range(0,len(scorer['parvalue']))
            elif search['key']=='par':
                search['pos']=range(0,len(scorer['par']))
        except Exception as e:
            traceback.print_exc()
            #print(scorer)
            #print(search)
            pdb.set_trace()
    print('converting model')
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
                        traceback.print_exc()
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
    print(nm)
    return nm
