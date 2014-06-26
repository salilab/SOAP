"""
   SOAP cross vailiation and bayesian predictive densitites module.

"""

from scorer import *
from sampling import *


class k2cv(object):
    """
    Cross validate base class
    """
    def __init__(self,spsfo=None,model=[]):
        if model:
            so=scorer(model=model)
            spsfo=optimizer(scorer=so)

        self.spsfo=copy.deepcopy(spsfo)
        self.ospsfo=spsfo
        self.model=spsfo.scorer.model
        if spsfo:
            self.sample=spsfo.scorer.get_sample()
        else:
            self.sample=[]
        self.nots=2


    def partfunc(self,osample,k,i):
        sample=copy.copy(osample)
        sl=len(sample)/k
        re=len(sample)%k
        random.seed(i*1234)
        random.shuffle(sample)
        ci=0
        samplelist=[]
        for i in range(0,k):
            if i<re:
                ll=sl+1
            else:
                ll=sl
            samplelist.append(sample[ci:ci+ll])
            ci=ci+ll
        return samplelist

    def evalfunc(self,sample,parvalue):
        #model should be a dictionary defines the type of the model,the parameters, the initial value of parameters, the value range
        #model={type:;par:;parvalue:;parrange:;parsearchbins:;}
        if not sample:
            return 0
        self.spsfo.scorer=self.testscorer
        validationresult=self.spsfo.calc_score_point(parvalue)
        return validationresult

    def trainfunc(self,sample):
        self.spsfo.scorer=copy.deepcopy(self.ospsfo.scorer)
        [self.spsfo.scorer,self.testscorer]=self.spsfo.scorer.get_ds_part(sample)
        trainresulttask=self.spsfo.multi_optimizer_cluster()#nspo.optimize()[0:2]#
        return trainresulttask

    def process_cv_result(self,resultlist):
        aresult=[]
#        pdb.set_trace()
        for result in resultlist:
            tperc=0
            perc=0
            for item in result:
                try:
                    tperc=tperc+item[0][1]
                    perc=perc+item[1]
                except Exception,e:
                    traceback.print_exc()
                    pdb.set_trace()
            tperc=tperc/len(result)
            perc=perc/len(result)
            aresult.append([tperc,perc])
        return aresult

    def pickfunc(self,cvresult):
        pass

    def cross_validation_single_run(self,input):
        if len(input[0])==0:
            return []
        trainedresulttask=self.trainfunc(input[0])
        return tasklist(tasklist=[trainedresulttask],afterprocessing=self.single_run_after_processing,other=copy.deepcopy(input[1]))

    def single_run_after_processing(self,trainedresult,tasklist=[],other=[]):
        testresult=self.evalfunc(other,trainedresult[0])
        return [trainedresult,testresult]

    def cross_validation_single_model(self,input):
        k=input['k']
        sample=input['sample']
        testsample=input['testsample']
        modellist=[]
        resultlist=[]
        trainresultlist=[]
        inputlist=[]
        for i in range(0,k):
            samplelist=self.partfunc(sample,2,i)
            inputlist.append([samplelist[0],samplelist[1]])
            inputlist.append([samplelist[1],samplelist[0]])
        inputlist.append([sample,testsample])
        #if self.spsfo.scorer.model['cvm'].startswith('parallel'):
        #    srtask=task(self.env).parallel_local_withreturn(self.cross_validation_single_run,inputlist,len(inputlist))
        #else:
        srtask=task().serial_local(self.cross_validation_single_run,inputlist)#parallel_local_withreturn,10
        return tasklist(tasklist=srtask)

    def get_test_sample(self,sample,testperc=0,testsample=[]):
        if testperc:
            return self.spsfo.scorer.get_sample(perc=testperc)
        else:
            trainsample=[]
            for item in sample:
                if not item in testsample:
                    trainsample.append(item)
            return [testsample,trainsample]

    def cross_validation_with_test(self,testsample=[],testperc=0):
        sample=self.sample
        #should call this function to do cross validation under every circustance
        testsample,trainsample=self.get_test_sample(sample,testperc,testsample)
        input={'k':self.model['cvk'],'sample':sample,'testsample':testsample}
        srtask=self.cross_validation_single_model(input)
        self.testsample=testsample
        return tasklist(tasklist=[srtask],afterprocessing=self.analyze_cv)

    def analyze_cv(self,rresultlist,tasklist=[],other=[]):
        resultlist=[]
        testresult=[]
        for result in rresultlist:
            if 1:
                testresult.append(result[-1])
                trainresult=result[:-1]
            resultlist.append(trainresult)
        return [resultlist,self.process_cv_result(resultlist), testresult]

    def load_optarraylist(self,path):
        cr=scipy.io.loadmat(path+'cvresult.mat')
        self.optarraylist=[]
        rk=[item for item in cr if item.startswith('a')]
        for i in range(0,len(rk)):
            key='a'+str(i)
            na=cr[key]
            nref=[]
            for i in range(0,len(na)):
                nref.append(list(na[i]))
            nref=squeeze_list(nref)
            na=list2array(nref)
            self.optarraylist.append(na)

    def combine_optimization_results(self,optresultlist):
        #analyze the best result from different runs
        na=np.zeros([len(optresultlist),len(optresultlist[0])])
        optarraylist=[]
        for i in range(len(optresultlist)):
            item=optresultlist[i]
            for j in range(len(item)):
                na[i,j]=item[j][:,-1].max()
            optarraylist.append(np.vstack(item))
        self.optarraylist=optarraylist
        self.optstats=na
        print self.optstats
        self.statsstr='('+str(self.optstats[-1].min())+' '+str(self.optstats[-1].max())+') '+str(self.optstats[-1].mean())+u"\u00B1"+str(self.optstats[-1].std())

    def get_scorer_list(self):
        if 'filter' in self.spsfo.scorer.model:
            fl=self.spsfo.scorer.model['filter']
            cvmodel=pickle.load(open(fl['filterpath']))
            if len(cvmodel['allresult'])!=len(self.inputlist):
                raise Exception('The model filter can not be applied to this  cv case, number of tries different')
        self.scorerlist=[]
        for i in range(0,len(self.inputlist)):
            scorer1,scorer2=self.spsfo.scorer.get_ds_part(self.inputlist[i][0])
            if 'filter' in self.spsfo.scorer.model:
                scorer1.loadfilter=[{'criteria':fl['criteria'],'parvalue':cvmodel['allresult'][i]['bestpar']}]
            self.scorerlist.append([scorer1,scorer2])

    def clustering_optimization_results(self):
        parlist=[]
        resultlist=[]
        for optarray in self.optarraylist:
            ra,nay=self.spsfo.get_rrf(optarray)
            parlist.append(ra)
            resultlist.append(nay)
        self.clusters=cvclustering(parlist,resultlist,self.scorerlist,clustermethod=self.model['clustermethod'])
        self.clusters.analyze()

        sd={'lastone':lastone,'parlist':parlist,'resultlist':resultlist,'testresultlist':testresultlist,
            'bestpar':bestpar, 'bestresult':resultlist[0],'testresult':testscorer.assess(bestpar),
            'fulltestresult':testscorer.assess(bestpar,report='full'),
            'pd':pd/k,'rlpd':rlpd}

    def analyze_clustering_results(self):
        if self.withlastone:
            self.finaloptresult=self.rdlist[-1]['bestresult']
            self.finaloptpar=self.rdlist[-1]['bestpar']
            if 'finalbestresult' in self.rdlist[-1]:
                self.finalfinaloptresult=self.rdlist[-1]['finalbestresult']
                self.finalfinaltestresult=self.rdlist[-1]['finaltestresult']
            else:
                self.finalfinaloptresult=0
                self.finalfinaltestresult=0
            tn=len(self.rdlist)-1
            self.lasttestresult=self.rdlist[-1]['testresult']
            self.bestmodelresult=self.rdlist[-1]['bestresult']
        else:
            self.finaloptresult=0
            self.finaloptpar=[]
            tn=len(self.rdlist)
            self.lasttestresult=0
            self.bestmodelresult=[]
            self.finalfinaloptresult=0
            self.finalfinaltestresult=0
        optlist=np.zeros(tn)
        testlist=np.zeros(tn)
        finaloptlist=np.zeros(tn)
        finaltestlist=np.zeros(tn)
        pdlist=np.zeros(tn)
        rlpdlist=np.zeros(tn)
        for i in range(0,tn):
            optlist[i]=self.rdlist[i]['bestresult']
            testlist[i]=self.rdlist[i]['testresult']
            if 'finalbestresult' in self.rdlist[-1]:
                finaloptlist[i]=self.rdlist[i]['finalbestresult']
                finaltestlist[i]=self.rdlist[i]['finaltestresult']
            else:
                finaloptlist[i]=0
                finaltestlist[i]=0
            pdlist[i]=self.rdlist[i]['pd']
            rlpdlist[i]=self.rdlist[i]['rlpd']
        self.optresults=optlist
        self.testresults=testlist
        self.optmean=optlist.mean()
        self.optstd=optlist.std()/np.sqrt(tn)
        self.finaloptmean=finaloptlist.mean()
        self.finaloptstd=finaloptlist.std()/np.sqrt(tn)
        self.pdsum=pdlist.mean()
        self.rlpdmean=rlpdlist.mean()
        self.testmean=testlist.mean()
        self.teststd=testlist.std()/np.sqrt(tn)
        self.finaltestmean=finaltestlist.mean()
        self.finalteststd=finaltestlist.std()/np.sqrt(tn)
        self.resultstr=u"{5:.5f},{6:.3f}, {0:2.3f}\u00B1{1:2.3f}, {2:2.3f}\u00B1{3:2.3f}, {4:2.3f}, {7:2.3f}__{8:2.3f}\u00B1{9:2.3f}, {10:2.3f}\u00B1{11:2.3f}, {12:2.3f}, {13:2.3f} ".format(
            self.optmean,self.optstd,self.testmean,self.teststd,self.finaloptresult,self.pdsum,self.rlpdmean,self.lasttestresult, self.finaloptmean,self.finaloptstd,self.finaltestmean,self.finalteststd,self.finalfinaloptresult,self.finalfinaltestresult)
        self.resultarray=np.array([self.optmean,self.optstd,self.testmean,self.teststd,self.finaloptresult,self.pdsum,self.lasttestresult, self.finaloptmean,self.finaloptstd,self.finaltestmean,self.finalteststd,self.finalfinaloptresult,self.finalfinaltestresult,self.rlpdmean])
        print self.resultstr

    def save_optimization_results(self):
        bdir=runenv.basedir+'results/'
        os.chdir(bdir)
        tdir=bdir+'_'.join(self.spsfo.scorer.bmtype['dslist'])
        if not os.path.isdir(tdir):
            os.mkdir(tdir)
        tdir=tdir+'/'+self.spsfo.scorer.bmtype['criteria']
        if not os.path.isdir(tdir):
            os.mkdir(tdir)
        self.logdir=tdir+'/'
        tdir=tdir+'/runs/'
        if not os.path.isdir(tdir):
            os.mkdir(tdir)
        cdp=tdir+'/'+self.rundir+'-'+self.resultstr.replace(u"\u00B1","+").replace(', ','_')+'/'
        os.mkdir(cdp)
        os.chdir(cdp)
        #self.clusters.figpath=cdp
        self.logpath=cdp
        #mypickle().dump(self.clusters, 'figinput')
        fh=open('cvmodel.pickle','w')
        pickle.dump({'model':self.spsfo.scorer.originalmodel, 'runtime':self.task.runduration,
                     'inputlist':self.inputlist,'allresult':self.rdlist,'originalresult':self.originalresult},fh)
        fh.close()
        #self.clusters.dump('cvcluster.pickle')
        print os.system('cp '+sys.path[1]+'/'+sys.argv[0]+' ./')
        if self.withlastone:
            self.bestmodel=self.spsfo.scorer.build_model(self.finaloptpar)
            fh=open('bestmodel.pickle','w')
            pickle.dump(self.bestmodel,fh)
            fh.close()
            fh=open('bestmodel','w')
            fh.write(any2str(self.bestmodel))
            fh.close()
            fh=open('bestmodelresult','w')
            fh.write(str(self.bestmodelresult)+'\n'+str(self.rdlist[-1]['fulltestresult']))
            fh.close()
        else:
            self.bestmodel=self.spsfo.scorer.originalmodel
            fh=open('bestmodel.pickle','w')
            pickle.dump(self.bestmodel,fh)
            fh.close()
            fh=open('bestmodel','w')
            fh.write(any2str(self.bestmodel))
            fh.close()
        self.write2logshelve()
        self.write2db()
        return self.resultstr

    def plot_bestpar(self):
        par=self.finaloptpar
        bestmodel=self.bestmodel
        so=scorer(model=bestmodel)
        so.assess(par)
        plt.clf()
        for i in len(so.model['scorers']):
            if so.model['scorers']['type']!='sf':
                continue
            sfo=so.scorerlist[i]
        plt.plot(sfo.sfv)


    def write2logshelve(self,model=None):
        if model==None:
            model=self.spsfo.scorer.originalmodel
        os.chdir(self.logdir)
        with FileLock("log.shelve", timeout=100, delay=2) as lock:
            print("Lock acquired.")
            resultdictlog=shelve.open(self.logdir+'log.shelve')
            del model['str']
            model['str']=any2str(model)
            resultdictlog[model['str']]=[self.resultarray,self.resultstr,self.rundir,self.bestmodel,self.bestmodelresult,self.logpath]
            resultdictlog.close()
            print("lock released")


    def modelinlog(self,model=None):
        if model==None:
            model=self.spsfo.scorer.originalmodel
        os.chdir(self.logdir)
        if 'str' in model:
            del model['str']
        model['str']=any2str(model)
        with FileLock("log.shelve", timeout=100, delay=2) as lock:
            print("Lock acquired.")
            if not os.path.isfile(self.logdir+'log.shelve'):
                nmr=[]
            else:
                resultdictlog=shelve.open(self.logdir+'log.shelve')
                if model['str'] in resultdictlog:
                    nmr=resultdictlog[model['str']]
                    #ko=k2cvcluster()
                    #rundirname=[item for item in os.listdir(os.path.join(self.baselogdir,'runs')) if item.startswith(nmr[2])][0]
                    #ko.load_fromlogdir(os.path.join(self.baselogdir,'runs',rundirname))
                    #nmr[0]=ko.resultarray
                else:
                    nmr=[]
                resultdictlog.close()
            print("lock released")
        if len(nmr)==0:
            return False
        else:
            self.resultarray,self.resultstr,self.rundir,self.bestmodel,self.bestmodelresult,self.logpath=nmr
            return True

    def write2db(self):
        #create the database if non exist
        import sqlite3
        class integerList(list):
            pass
        class textList(list):
            pass
        class realList(list):
            pass
        class anyList(list):
            pass

        def adapt_slist(s):
            return '; '.join([str(item) for item in s])


        def adapt_anylist(textList):
            return '; '.join([pickle.dumps(item) for item in textList])

        def convert_integerList(s):
            return map(int, s.split("; "))

        def convert_realList(s):
            return map(float, s.split("; "))

        def convert_textList(s):
            return s.split("; ")

        def convert_anyList(s):
            return map(pickle.loads, s.split("; "))
        sqlite3.register_adapter(integerList, adapt_slist)
        sqlite3.register_converter("integerlist", convert_integerList)
        sqlite3.register_adapter(realList, adapt_slist)
        sqlite3.register_converter("reallist", convert_realList)
        sqlite3.register_adapter(textList, adapt_slist)
        sqlite3.register_converter("textlist", convert_textList)
        sqlite3.register_adapter(anyList, adapt_anylist)
        sqlite3.register_converter("anylist", convert_anyList)
        con = sqlite3.connect(runenv.resultdbpath,detect_types=sqlite3.PARSE_DECLTYPES)
        conn= con.cursor()
        #self.optmean,self.optstd,self.testmean,self.teststd,
        #self.finaloptresult,self.lasttestresult,
        #self.finaloptmean,self.finaloptstd,self.finaltestmean,self.finalteststd,
        #self.finalfinaloptresult,self.finalfinaltestresult,
        #self.rlpdmean,self.pdsum,
        conn.execute("""create table if not exists Result (result TEXT,
                     dslist textlist,type TEXT,criteria TEXT, bm TEXT,combine TEXT,
                     filters TEXT, criteria2 TEXT, scorers integerlist, optmethod INTEGER,
                     optmean REAL, optstd REAL,testmean REAL, teststd REAL, finaloptresult REAL,finaltestresult REAL,
                     optmean2 REAL, optstd2 REAL,testmean2 REAL, teststd2 REAL, finaloptresult2 REAL,finaltestresult2 REAL,
                     pdresultmean REAL, pd REAL,bestpar reallist,runnum INTEGER primary key);""")
        conn.execute("""create table if not exists Scorers (scoind INTEGER, runnum INTEGER)""")
        conn.execute("""create table if not exists OptMethod (optind INTEGER unique, optmethod TEXT primary key, sml TEXT, cvk INTEGER,
                     repeat INTEGER, testperc REAL,fold INTEGER);""")
        conn.execute("""create table if not exists Scorer (scoind INTEGER unique, scorer any primary key, type TEXT, features TEXT,
                     sftype TEXT, par TEXT, searchpar TEXT, parvalue TEXT,  searchparvalue TEXT, ratio TEXT,searchratio TEXT,
                     bm TEXT, pdbset TEXT, pm TEXT, genmethod TEXT);""")
        conn.execute("""create table if not exists Counter (numoptm INTEGER, numsco INTEGER, numres INTEGER)""")
        numofcount=conn.execute("""select count(*) from Counter""").fetchall()[0][0]
        if numofcount==0:
            conn.execute("insert into Counter values (1000000,2000000,3000000)")
        elif numofcount>1:
            print "more than 1 line of count"
            pdb.set_trace()

        if self.spsfo==None:
            m=self.originalmodel
        else:
            m=self.spsfo.scorer.originalmodel
        optList=[]
        for key in ['sml','cvk','repeat','testperc','fold']:
            if key in m:
                optList.append(m[key])
            else:
                optList.append(None)
        optList[0]=any2str(optList[0])
        optList.insert(0,any2str(optList))
        conn.execute("""SELECT *
                        FROM OptMethod
                        WHERE optmethod=?;""",(optList[0],))
        rows=conn.fetchall()
        if len(rows)>1:
            print "duplicate in database"
            pdb.set_trace()
        elif len(rows)==1:
            if rows[0][1:]!=tuple(optList):
                print "duplicate but not matching optmethod"
                pdb.set_trace()
            optList.insert(0, rows[0][0])
        elif len(rows)==0:
            conn.execute('update Counter set numoptm=numoptm+1')
            optind=conn.execute('Select numoptm from Counter').fetchall()[0][0]
            optList.insert(0,optind)
            conn.execute('INSERT INTO OptMethod values (?,?,?,?,?,?,?)',tuple(optList))
        for search in m['searches']:
            search['object'][search['key']]=['search',len(search['object'][search['key']]),any2str(search['InitialGenerator'])]
        allScorers=[]
        for scorer in m['scorers']:
            scoList=[]
            scoList.append(any2str(scorer))
            for key in ['type','features','sftype','par','parvalue','ratio','bm','pdbset','pm','genmethod']:
                if key in scorer:
                    if isinstance(scorer[key],list) and len(scorer[key])>0 and scorer[key][0]=='search':
                        scoList.append(any2str(scorer[key][1]))
                        scoList.append(scorer[key][2])
                    elif key in ['par','parvalue','ratio']:
                        scoList.append(any2str(scorer[key]))
                        scoList.append(None)
                    else:
                        scoList.append(any2str(scorer[key]))
                elif key in ['par','parvalue','ratio']:
                    scoList.append(None)
                    scoList.append(None)
                else:
                    scoList.append(None)
            conn.execute("""SELECT *
                            FROM Scorer
                            WHERE scorer=?""",(scoList[0],))
            rows=conn.fetchall()
            if len(rows)>1:
                print "duplicate in database"
                pdb.set_trace()
            elif len(rows)==1:
                if rows[0][1:]!=tuple(scoList):
                    print "duplicate but not matching scorer"
                    pdb.set_trace()
                scoList.insert(0, rows[0][0])
            elif len(rows)==0:
                conn.execute('update Counter set numsco=numsco+1')
                scoind=conn.execute('Select numsco from Counter').fetchall()[0][0]
                scoList.insert(0,scoind)
                conn.execute('INSERT INTO Scorer values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',tuple(scoList))
            allScorers.append(scoList[0])
        resList=[]
        resList.append(any2str(m))
        resList.append(textList(m['bmtype']['dslist']))
        for key in ['type','criteria','bm','combine','filters','finalcriteria']:
            if key in m['bmtype']:
                resList.append(any2str(m['bmtype'][key]))
            else:
                resList.append(None)
        resList.append(textList(allScorers))
        resList.append(optList[0])
        resList.extend([self.optmean,self.optstd,self.testmean,self.teststd,self.finaloptresult,self.lasttestresult,
            self.finaloptmean,self.finaloptstd,self.finaltestmean,self.finalteststd,self.finalfinaloptresult,self.finalfinaltestresult,
            self.rlpdmean,self.pdsum, self.rdlist[-1]['bestpar'],int(self.rundir)])
        if  isinstance(resList[-2],np.float):
            resList[-2]=[resList[-2]]
        resList[-2]=realList(resList[-2])
        if isinstance(resList[-2][0],np.ndarray):
            nl=realList()
            for i in range(len(resList[-2][0])):
                nl.append(resList[-2][0][i])
            resList[-2]=nl
        print resList[-2]
        conn.execute("""SELECT *
                        FROM Result
                        WHERE runnum=?""",(resList[-1],))
        rows=conn.fetchall()
        if len(rows)>1:
            print "duplicate in database"
        elif len(rows)==1:
            for i in range(len(resList)):
                try:
                    if np.isnan(resList[i]):
                        resList[i]=None
                except:
                    continue
            if rows[0][:-2]!=tuple(resList)[:-2] or rows[0][-1]!=tuple(resList)[-1] or (not np.allclose(rows[0][-2],resList[-2])):
                print "duplicate but not matching result"
                pdb.set_trace()
        elif len(rows)==0:
            for sc in allScorers:
                conn.execute('INSERT INTO Scorers values (?,?)',(sc,int(self.rundir)))
            conn.execute('INSERT INTO Result values (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)',tuple(resList))
        con.commit()
        con.close()

class k2cvlocal(k2cv):
    """
    Cross validate a model locally // may not work on its own anymore
    """
    def trainfunc(self,sample):
        #pdb.set_trace()
        self.spsfo.scorer=copy.deepcopy(self.ospsfo.scorer)
        [self.spsfo.scorer,self.testscorer]=self.spsfo.scorer.get_ds_part(sample)
        trainresult=self.spsfo.multi_optimizer_local()#nspo.optimize()[0:2]#
        return trainresult

    def cross_validation_single_model(self,input):
        k=input['k']
        sample=input['sample']
        testsample=input['testsample']
        modellist=[]
        resultlist=[]
        trainresultlist=[]
        inputlist=[]
        for i in range(0,k):
            samplelist=self.partfunc(sample,2,i)
            inputlist.append([samplelist[0],samplelist[1]])
            inputlist.append([samplelist[1],samplelist[0]])
        inputlist.append([sample,testsample])
        #if self.spsfo.scorer.model['cvm'].startswith('parallel'):
        #    srtask=task(runenv).parallel_local_withreturn(self.cross_validation_single_run,inputlist,len(inputlist))
        #else:
        srtask=task().serial_local(self.cross_validation_single_run,inputlist)#parallel_local_withreturn,10
        return srtask

    def cross_validation_single_run(self,input):
        if len(input[0])==0:
            return []
        trainedresult=self.trainfunc(input[0])
        return self.single_run_after_processing(trainedresult)

    def cross_validation_with_test(self,testsample=[],testperc=0):
        sample=self.sample
        #should call this function to do cross validation under every circustance
        testsample,trainsample=self.get_test_sample(sample,testperc,testsample)
        input={'k':self.model['cvk'],'sample':sample,'testsample':testsample}
        self.testsample=testsample
        srtask=self.cross_validation_single_model(input)
        return self.analyze_cv([srtask])

class k2cvcluster(k2cvlocal):
    """
    Cross validate a model on SGE cluster
    """
    def __init__(self,spsfo=None,model=[], initialize=False,logpath=''):
        self.model=model
        self.spsfo=spsfo
        self.runspernode=1
        self.runpath='./'
        self.statsstr=''
        self.testperc=-1
        self.testsample=[]
        self.rid=0
        self.rsn=0
        self.repeat=1
        if logpath:
            self.load_fromlogdir(logpath)
        if initialize:
            self.initialize_model()
        if runenv.hostn==0 and len(model)>0:
            try:
                bdir=runenv.basedir+'results/'
                tdir=bdir+'_'.join(self.model['bmtype']['dslist'])
                if not os.path.isdir(tdir):
                    os.mkdir(tdir)
                tdir=tdir+'/'+self.model['bmtype']['criteria']
                if not os.path.isdir(tdir):
                    os.mkdir(tdir)
                self.logdir=tdir+'/'
            except:
                traceback.print_exc()
                pdb.set_trace()

    def initialize_scorer(self):
        so=scorer(model=self.originalmodel)
        spsfo=optimizer(scorer=so)
        self.spsfo=spsfo
        self.model=self.originalmodel
        self.initialize_model()
        self.prepare_cross_validation_sample()
        self.get_scorer_list()

    def load_fromlogdir(self,lp=''): #note that repeat is ignored here
        cvmodel=pickle.load(open(os.path.join(lp,'cvmodel.pickle')))
        self.rdlist=cvmodel['allresult']
        if 'testperc' in cvmodel['model']:
            self.withlastone=True
        else:
            self.withlastone=False
        self.originalmodel=cvmodel['model']
        self.spsfo=None
        if 0:
            bm=pickle.load(open(os.path.join(lp,'bestmodel.pickle')))
            if 'repeat' in bm:
                bm['repeat']=1
            so=scorer(model=bm)
            spsfo=optimizer(scorer=so)
            self.spsfo=spsfo
            self.model=bm
            self.initialize_model()
            self.prepare_cross_validation_sample()
        self.analyze_clustering_results()
        #self.rundir=os.path.split(lp)[-1].split('-')[0]

    def initialize_model(self):
        if self.model!=[]:
            so=scorer(model=self.model)
            self.spsfo=optimizer(scorer=so)
            if 'testperc' in self.model:
                self.testperc=self.model['testperc']
            if 'testsample' in self.model:
                self.testsample=self.model['testsample']
        if self.spsfo:
            #self.ospsfo=optimizer(scorer=self.spsfo.scorer)
            self.sample=self.spsfo.scorer.get_sample()
            self.model=self.spsfo.scorer.model
        else:
            self.sample=[]

    def cross_validation_single_run(self,input):
        pass

    def cross_validation_single_model(self,input):
        k=input['k']
        sample=input['sample']
        modellist=[]
        resultlist=[]
        trainresultlist=[]
        inputlist=[]
        parnum=2
        if 'fold' in self.model:
            parnum=self.model['fold']
        for i in range(0,k):
            samplelist=self.partfunc(sample,parnum,i)
            for j in range(parnum):
                trainlist=[]
                for m in range(parnum):
                    if m!=j:
                        trainlist=trainlist+samplelist[m]
                inputlist.append([trainlist,samplelist[j]])
        if 'testsample' in input:
            inputlist.append([sample,input['testsample']])
        if self.repeat>1:
            self.inputlist=[]
            for i in range(self.repeat):
                self.inputlist=self.inputlist+inputlist
        else:
            self.inputlist=inputlist
        self.get_scorer_list()
        self.distribute_runs()
        if "initialmodelpath" in self.model:
            self.get_initial_values_from_path(self.model['initialmodelpath'])

    def get_initial_values_from_path(self,path):
        if isinstance(path,str):
            path=[path]
        csp=0
        for currentpath in path:
            cvm=pickle.load(open(os.path.join(currentpath,'cvmodel.pickle')))
            om=cvm['model']
            #if om['cvk']!=self.model['cvk'] or ('fold' in self.model and om['fold']!=self.model['fold']) or ('testperc' in self.model and om['testperc']!=self.model['testperc']):
            #    raise Bugs('Models does not match, can not load the model as initial conditions')
            #match models How to???
            if 'initialpars' in self.model and self.model['initialpars']=='best':
                parname='bestrepna'
            else:
                parname='repna'
            if len(cvm['allresult'])<(len(self.scorerlist)/self.repeat):
                rdlist=[]
                for i in range(len(self.scorerlist)/self.repeat):
                    k=i%len(cvm['allresult'])
                    crd=cvm['allresult'][k]
                    a=range(crd[parname].shape[0])
                    random.shuffle(a)
                    crd[parname]=crd[parname][a,:]
                    rdlist.append(crd)
            else:
                rdlist=cvm['allresult'][:len(self.scorerlist)/self.repeat]
            oso=scorer(model=om)
            for i in range(len(self.scorerlist)/self.repeat):
                repna=rdlist[i][parname]
                sp=0
                trn=self.runsperscorer*self.repeat/repna.shape[0]+1
                repna=np.vstack([repna for kk in range(trn)])[0:self.runsperscorer*self.repeat,:]
                for j in range(self.repeat):
                    so=self.scorerlist[j*len(rdlist)+i][0]
                    ep=sp+self.runsperscorer
                    if len(so.initialvalues)==0:
                        so.initialvalues=np.zeros((self.runsperscorer,len(self.spsfo.scorer.parvalues)))-9999
                    csptemp=self.get_parvalue_from_rf(so.initialvalues,om,repna[sp:ep,:-1],oso,sp=csp)
                    sp=sp+self.runsperscorer
            csp=csptemp
        print "please check the initialvalue performance"
        self.spsfo.scorer.initialvalues=so.initialvalues
        self.spsfo.get_initial_value()
        print self.spsfo.scorer.assess(self.spsfo.bestpar)

    def get_parvalue_from_rf(self,newpar,originalmodel,originalpar,oso,sp=0):
        #originalpar is the list containing all the search pars, we need to somehow convert the original pars to the pars for the new model
        #pdb.set_trace()
        parshape=originalpar.shape
        #the scorers from the original model and current model are matched 1 by 1.
        ospi=0 # original model  position
        for spi in range(sp,len(self.spsfo.scorer.model['scorers'])): #current model posiiotn
            print spi
            #if spi==(len(self.spsfo.scorer.model['scorers'])-2):
            #    pdb.set_trace()
            if len(self.spsfo.scorer.scorersearchlist[spi])==0:
                continue
            bins=[]
            if 'par' in self.model['scorers'][spi]:
                for j in range(parshape[0]):
                    bins.append(self.model['scorers'][spi]['par'])
            #skip original model to the first search position
            while ospi<len(originalmodel['scorers']):
                if len(oso.scorersearchlist[ospi])>0:
                    break
                else:
                    ospi+=1
            if ospi==len(originalmodel['scorers']):
                spi-=1
                break
            setvalues=False
            for si in self.spsfo.scorer.scorersearchlist[spi]:
                currentsearch=self.model['searches'][si]
                if currentsearch['key']=='ratio':
                    for si2 in oso.scorersearchlist[ospi]:
                        if originalmodel['searches'][si2]['key']=='ratio':
                            if np.any(newpar[:,self.spsfo.scorer.searchlistspos[si]]!=-9999):
                                print "resetting the values twice"
                            for j in range(parshape[0]):
                                newpar[j,self.spsfo.scorer.searchlistspos[si]]=originalpar[j,oso.searchlistspos[si2]]
                            setvalues=True
                elif currentsearch['key']=='par':
                    setbins=False
                    for si2 in oso.scorersearchlist[ospi]:
                        if originalmodel['searches'][si2]['key']=='par' and ((self.spsfo.scorer.searchlistspos[si+1]-self.spsfo.scorer.searchlistspos[si])==(oso.searchlistspos[si2+1]-oso.searchlistspos[si2])):
                            setbins=True
                            bins=[]
                            if np.any(newpar[:,self.spsfo.scorer.searchlistspos[si]:self.spsfo.scorer.searchlistspos[si+1]]!=-9999):
                                print "resetting the values twice"
                                pdb.set_trace()
                            for j in range(parshape[0]):
                                newpar[j,self.spsfo.scorer.searchlistspos[si]:self.spsfo.scorer.searchlistspos[si+1]]=originalpar[j,oso.searchlistspos[si2]:oso.searchlistspos[si2+1]]
                                bins.append(originalpar[j,oso.searchlistspos[si2]:oso.searchlistspos[si2+1]])
                            setvalues=True
                    if setbins==False and len(originalmodel['scorers'][ospi]['par'])==len(self.model['scorers'][spi]['par']):
                        bins=[]
                        if np.any(newpar[:,self.spsfo.scorer.searchlistspos[si]:self.spsfo.scorer.searchlistspos[si+1]]!=-9999):
                            print "resetting the values twice"
                            pdb.set_trace()
                        for j in range(parshape[0]):
                            newpar[j,self.spsfo.scorer.searchlistspos[si]:self.spsfo.scorer.searchlistspos[si+1]]=originalmodel['scorers'][ospi]['par']
                            bins.append(originalmodel['scorers'][ospi]['par'])
                        setvalues=True
            for si in self.spsfo.scorer.scorersearchlist[spi]:
                currentsearch=self.model['searches'][si]
                if currentsearch['key']=='parvalue':
                    for si2 in oso.scorersearchlist[ospi]:
                        if originalmodel['searches'][si2]['key']=='parvalue':
                            if np.any(newpar[:,self.spsfo.scorer.searchlistspos[si]:self.spsfo.scorer.searchlistspos[si+1]]!=-9999):
                                print "resetting the values twice"
                                pdb.set_trace()
                            for j in range(parshape[0]):
                                oso.parvalues=originalpar[j,:]
                                oso.assign_values2model()
                                sfmodel=oso.model['searches'][si2]['object']
                                sfo=sf(**sfmodel)
                                if self.model['scorers'][spi]['sftype']=='bins':
                                    sfo.bins=sf(**self.model['scorers'][spi]).bins
                                    sfv=np.exp(sfo.get_sf(returnreft=True))
                                    newpar[j,self.spsfo.scorer.searchlistspos[si]:self.spsfo.scorer.searchlistspos[si+1]]=sfv
                                else:
                                    sfo.bins=[list(bins[j])+[sfo.bins[0][-1]]]
                                    sfv=np.exp(sfo.get_sf(returnreft=True))
                                    newpar[j,self.spsfo.scorer.searchlistspos[si]:self.spsfo.scorer.searchlistspos[si+1]]=sfv[:-1]
                            setvalues=True
            if setvalues:
                ospi+=1
            if ospi==(len(originalmodel['scorers'])):
                break
            #pdb.set_trace()
        if np.any(newpar==-9999):
            print 'Some Initial values are undefined '
        print spi
        spi+=1
        return spi

    def distribute_runs(self):
        #get runs per node
        self.optn=len(self.inputlist) #total number of trials, number of different sets
        perruntime=len(self.spsfo.scorer.dsss.ds.sa)/250000.0 #approximate run time in minutes
        optin=self.spsfo.scorer.model['runsperscorer']
        self.runsperscorer=optin
        rsnf=optimizerjobtime/perruntime
        if (rsnf+1)>=optin:
            self.runspernode=optin
        elif rsnf>1:
            for i in range(int(np.rint(rsnf+1)),1,-1):
                if optin%i==0:
                    self.runspernode=i
                    break
        self.runspernode=1
        if optin%self.runspernode!=0:
            raise Exception('Make sure the total number of optimization tries is equal to interger * '+str(self.runspernode))
        rpert=optin/self.runspernode #number of runs per trail
        k=0
        rnl=[]
        rnd={}
        for i in range(self.optn):
            if self.withlastone and i==self.optn-1:
                rpn=np.rint(self.runspernode/2.0)#
                if rpn<1:
                    rpn=1
            else:
                rpn=self.runspernode
            rnl.append([])
            for j,rj in zip(range(optin),range(1,optin+1)):
                nn=j%rpn
                if nn==0:
                    k=k+1 # run number, sort of series number
                    rnl[i].append(k)
                    rnd[k]=[i,[rj]]
                else:
                    rnd[k][1].append(rj)
        self.nors=k#number of runs
        self.runnumlist=rnl
        self.runnumdict=rnd

    def get_task(self):
        if self.modelinlog(self.model):
            return 0
        else:
            self.task=task('','',afterprocessing=self,preparing=self)
            return self.task

    def prepare_cross_validation_sample(self):
        sample=self.sample
        #should call this function to do cross validation under every circustance
        try:
            if 'repeat' in self.model:
                self.repeat=self.model['repeat']
            if self.testperc>=0:
                testsample,trainsample=self.get_test_sample(sample,self.testperc,self.testsample)
                self.testsample=testsample
                self.trainsample=trainsample
                input={'k':self.model['cvk'],'sample':trainsample,'testsample':testsample}
                self.withlastone=True
            else:
                input={'k':self.model['cvk'],'sample':sample}
                self.withlastone=False
        except:
            pdb.set_trace()
        self.cross_validation_single_model(input)
        #pdb.set_trace()
        #preparing task object

    def cross_validation_with_test(self):
        to=self.get_task()
        tasklist(tasklist=[to]).monitor2end()

    def runtask_cluster(self,inputcode):
        #the individual nodes on the cluster will run this function, if other funtions needs to be run here, chage the code here.
        self.rid=int(inputcode)
        self.path='./'
        self.runpath=runenv.basedir
        self.inputcode=str(inputcode)
        inputindex=self.runnumdict[int(self.inputcode)][0]
        self.runstatus=np.zeros(len(self.runnumlist[inputindex]))
        self.spsfo.scorer=self.scorerlist[inputindex][0]
        self.testscorer=self.scorerlist[inputindex][1]
        self.spsfo.scorerid=inputindex
        self.spsfo.jobid=inputcode
        self.spsfo.runsperscorer=len(self.runnumlist[inputindex])
        res=[]
        #initialize the daemon thread array if specified...
        if 'thread' in self.spsfo.scorer.model and self.spsfo.scorer.model['thread']>0:
            runenv.numofthread=self.spsfo.scorer.model['thread']
            runenv.setup_threadpool()
        #thread pool setup finished
        self.trainscorer=self.scorerlist[inputindex][0]
        self.testscorer=self.scorerlist[inputindex][1]
        if not os.path.isfile(self.runpath+self.inputcode+'.pickle'):
            self.spsfo.starttime=time.time()
            self.spsfo.quit=False
            for i in self.runnumdict[int(self.inputcode)][1]:
                self.spsfo.scorer=load_scorer(self.scorerlist[inputindex][0],self.spsfo.indexlen)
                res.append(self.spsfo.runtask_local(i,self.testscorer))
                if self.spsfo.quit:
                    break
            print os.system('touch Finished'+str(self.spsfo.scorerid)+'#'+str(self.spsfo.jobid)+'#'+str(int(time.time()-self.spsfo.starttime)))
            print os.system('touch '+self.runpath+self.inputcode+'.pickle')
            fh=open(self.runpath+self.inputcode+'.pickle','w')
            pickle.dump(res,fh)
            fh.close()
        try:
            if self.runnumlist[inputindex][-1]==int(self.inputcode):#last one will collect all the jobs
                self.wait_analyze_single_try(inputindex)
                #oc.dump(self.inputcode+'.optimizer.pickle')
                for i in self.runnumlist[inputindex]:
                    print 'removing '+str(i)
                    print os.system('rm '+str(i)+'.pickle')
            else:
                fh=open(self.inputcode+'.optimizer.pickle','wb')
                pickle.dump('Empty',fh)
                fh.close()
            runsuccess=True
        except Exception,e:
            print e
            print "FatalError : can not write the pickle file"
            traceback.print_exc()
            runsuccess=False
        report_job_runstatus(self.runpath, runsuccess, inputcode, '.optimizer',inputname='runme.py')
        return 0

    def wait_analyze_single_try(self,inputindex):
        print 'waiting for other runs to finish '
        #rn=self.check_status_single_try(inputindex)
        #while rn>0:
        #    time.sleep(10)
        #    rn=self.check_status_single_try(inputindex)
        #    print os.system('touch '+self.runpath+str(inputindex)+'-waiting-'+str(rn))
        rn=self.check_status_single_try(inputindex,type='jobfinished')
        while rn>1: #change to rn>1
            time.sleep(30)
            rn=self.check_status_single_try(inputindex,type='jobfinished')
            print 'waiting for others to finish'
            #print os.system('touch '+self.runpath+str(inputindex)+'-waitingjob-'+str(rn))
        #print os.system('rm '+self.runpath+str(inputindex)+'-waiting*')
        return self.analyze_results_for_single_try(inputindex)

    def check_status_single_try(self,inputindex,type='pickle'):
        k=0
        fl=os.listdir(self.runpath)
        for i,ii in zip(self.runnumlist[inputindex],range(len(self.runnumlist[inputindex]))):
            if type=='pickle':
                cfn=str(i)+'.pickle'
            elif type=='jobfinished':
                cfn='jobfinished.'+str(i)+'.'+str(i)+'.optimizer.tar.gz'
            if self.runstatus[ii] or cfn in fl:
                self.runstatus[ii]=1
                k=k+1
            else:
                pass
                #print os.system('touch '+str(inputindex)+'-waitingnum-'+str(i))
        print "Finished run num ("+type+") :"+str(k)
        return len(self.runnumlist[inputindex])-k

    def analyze_results_for_single_try(self,inputindex):
        dl=[]
        brl=[]
        for i in self.runnumlist[inputindex]:
            tnt=3
            success=False
            while tnt>0 and not success:
                tnt=tnt-1
                try:
                    fh=open(self.runpath+str(i)+'.pickle')
                    ndl=pickle.load(fh)
                    fh.close()
                    for nd in ndl:
                        dl.append(nd) #result dictionary list
                        brl.append(np.array(nd['trainresultlist']).max())
                    success=True
                except:
                    traceback.print_exc()
                    success=False
                    time.sleep(1)
        if inputindex==len(self.runnumlist)-1 and self.withlastone:
            lastone=True # whether this is training on the whole set...
        else:
            lastone=False
        bbr=max(brl)
        parlist=[]
        resultlist=[]
        testresult=0
        pd=0
        k=0
        testresultlist=[]
        allparlist=[]
        alltrainresult=[]
        alltestresult=[]
        for i in range(len(brl)):
            if brl[i]==bbr: # only take runs that has reach "the global minimum"
                parlist=parlist+dl[i]['parlist']
                resultlist=resultlist+dl[i]['trainresultlist']
                testresultlist=testresultlist+dl[i]['testresultlist']
                pd=pd+dl[i]['PredictiveDensity']
                k=k+1
            allparlist=allparlist+dl[i]['parlist']
            alltrainresult=alltrainresult+dl[i]['trainresultlist']
            #alltestresult=alltestresult+dl[i]['testresultlist']
        #get representative pars
        allna=np.hstack([np.array(allparlist),np.array(alltrainresult).reshape(len(alltrainresult),1)])
        #bestrepna=get_representative_pars_forbest(allna,maxcluster=400,cutoffdistance=0.0000001,clusteringperc=0.999, othermaxnumber=0)
        #repna=get_representative_pars_forbest(allna,maxcluster=400,cutoffdistance=0.0000001,clusteringperc=0.5, othermaxnumber=0)
        plna=np.squeeze(np.array(parlist))
        if len(plna.shape)==1 and len(resultlist)==1:
            bestpar=plna
            bestresult=resultlist[0]
        else:
            ra=np.array(resultlist)
            plnam=plna.mean(axis=0)
            plnas=plna[ra==ra.max()]
            if len(plna.shape)==1:
                plnad=(((plnas-plnam))**2)
            else:
                plnad=(((plnas-plnam))**2).sum(axis=1)
            mindex=np.argmin(plnad)
            bestpar=plnas[mindex]
            bestresult=ra.max()
        ta=np.array(testresultlist)
        ta=-(ta-dl[0]['idealresult'])**2/(2*dl[0]['std']**2)
        tamax=ta.max()
        rlpd=tamax+np.log((np.exp(ta-tamax)).sum())-np.log(len(testresultlist))
        finalbestresult=0
        if 'finalcriteria' in self.model['bmtype']:
            self.spsfo.scorer=load_scorer(self.trainscorer,self.spsfo.indexlen)
            finalbestresult=self.spsfo.scorer.assess(bestpar,slevel=self.model['bmtype']['finalcriteria'])
            del self.spsfo.scorer
        testscorer=load_scorer(self.testscorer,self.spsfo.indexlen)
        finaltestresult=0
        if 'finalcriteria' in self.model['bmtype']:
            finaltestresult=testscorer.assess(bestpar,slevel=self.model['bmtype']['finalcriteria'])
        sd={'lastone':lastone,#'parlist':parlist,'resultlist':resultlist,'testresultlist':testresultlist,
            'bestpar':bestpar, 'bestresult':resultlist[0],'testresult':testscorer.assess(bestpar),
            'fulltestresult':testscorer.assess(bestpar,report='full'), 'allna':allna
            ,'pd':pd/k,'rlpd':rlpd,'finalbestresult':finalbestresult,'finaltestresult':finaltestresult}
        fh=open(self.inputcode+'.optimizer.pickle','w')
        pickle.dump(sd,fh)
        fh.close()

    def prepare_task_input(self):
        self.initialize_model()
        self.prepare_cross_validation_sample()
        self.dir=runenv.runs
        self.task.get_sn()
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
        self.save()
        freememory=0.5+self.spsfo.arraysize*self.get_runmemory_percentage()/1000000000.0
        self.freememory=freememory
        #if self.spsfo.arraysize>50000000:
        #    self.task.queues='-q lab.q'
        nors=self.nors
        #write the input list for different runs
        inputlist=open('inputlist','w')
        inputlist.write(','.join(['runme.py' for i in range(1,nors+1)]))
        inputlist.close()
        makemdt=open('runme.py','w')
        makemdt.write('keeptry=5\nwhile keeptry:\n    try:\n        from SOAP import *\n        keeptry=0\n    except:\n        keeptry-=1\n'
                      +'import sys \n \nspopt=k2cvcluster()\n'
                      +'spopt.freememory='+str(freememory)+'\n'
                      +'spopt.rsn='+str(self.task.rsn)+'\n'
                      +'spopt=spopt.load()\n'
                      +'spopt.runtask_cluster(sys.argv[1])')
        makemdt.flush()
        if self.model['bmtype']['type']=='sprefine':
            parallel=8
        else:
            parallel=False
        #if self.freememory>20:#this is only necessary when lots of people request memory but not used
            #generate_job_submit_script(0.2,self.rundir,runtime,nors,mem_total=10)
        #else:
        additional_resources=''
        if 'thread' in self.spsfo.scorer.model and self.spsfo.scorer.model['thread']>0:
            parallel=self.spsfo.scorer.model['thread']
            if self.freememory>10:
                additional_resources='\n#$ -l hostname=ih*|xe5520=true\n#$ -q lab.q\n'
        if not parallel:
            dp=1.0
        else:
            dp=float(parallel)
        generate_job_submit_script(freememory/dp,self.rundir,runtime,nors,parallel=parallel,additional_resources=additional_resources)
        return []

    def get_runmemory_percentage(self):
        ml=0
        for sp in self.scorerlist:
            for sc in sp:
                nindexlist=squeeze_indexlist(sc.dsss.ds.indexlist)
                if len(nindexlist)>0 and nindexlist[-1][-1]>ml:
                    ml=nindexlist[-1][-1]
        return ml/float(self.spsfo.indexlen)

    def afterprocessing(self):
        os.chdir(self.dir+self.rundir)
        self.rdlist=[[] for i in range(0,self.optn)]
        for i in range(self.optn):
            for inputn in self.runnumlist[i]:
                if os.path.isfile(str(inputn)+'.optimizer.pickle'):
                    fh=open(str(inputn)+'.optimizer.pickle','rb')
                    res=pickle.load(fh)
                    fh.close()
                    if res=='Empty':
                        continue
                    else:
                        self.rdlist[i]=res
        if self.repeat>1:
            rdll=[]
            nrdlist=[]
            srn=self.optn/self.repeat
            for i in range(self.repeat):
                rdll.append(self.rdlist[srn*i:srn*(i+1)])
            self.originalresult=rdll
            for i in range(srn):
                print "num"+str(i)
                bestresult=-np.inf #choose the best train restult
                print bestresult
                bestrd=[]
                for j in range(self.repeat):
                    if rdll[j][i]['bestresult']>bestresult:
                        bestresult=rdll[j][i]['bestresult']
                        print bestresult
                        bestrd=rdll[j][i]
                nrdlist.append(bestrd)
            self.rdlist=nrdlist
        else:
            self.originalresult=[self.rdlist]
        for i in range(len(self.rdlist)):
            rd=self.rdlist[i]
            allnalist=[]
            for arl in self.originalresult:
                allnalist.append(arl[i]['allna'])
            if 0:
                bestrepna=get_representative_pars_forbest(np.vstack(allnalist),maxcluster=400,cutoffdistance=0.0000001,clusteringperc=0.999, othermaxnumber=0)
                repna=get_representative_pars_forbest(np.vstack(allnalist),maxcluster=800,cutoffdistance=0.0000001,clusteringperc=0.5, othermaxnumber=0)
                rd['bestrepna']=bestrepna
                rd['repna']=repna
        #pdb.set_trace()
        #self.clusters=cvclustering([],[],[],clustermethod=self.model['clustermethod'])
        #self.clusters.optclusterlist=optclusterlist
        #self.clusters.optclusterlist[-1].optscorer=self.scorerlist[-1][0]
        self.analyze_clustering_results()
        cvresult=self.save_optimization_results()
        print os.system('rm -r '+self.dir+self.rundir)
        #del self.clusters
        del self.task
        gc.collect()
        return cvresult

    def save(self):
        self.spsfo.arraysize=mypickle().dump([self.spsfo.scorer,self.scorerlist],'input')
        fh=open('input.pickle','w')
        pickle.dump(self,fh)
        fh.close()
        # the scorer object will be useless unlesed the savednumpy is loaded from numpy arrays.

    def load(self):
        os.chdir(self.runpath)
        #fm=get_system_freemem()
        #if fm/1000.0<self.freememory:
        #    sys.exit(0)
        if self.freememory<20 or runenv.hostn==0:
            return mypickle().loadpickleonly('input')
        else:
            if runenv.hostn==1:
                time.sleep(self.rid*0.2)
                fl=get_starting_jobs()
                while len(fl)>0:
                    sleep(5)
                    print "jobs are starting on  "+hcn
                fm=get_system_freemem()
                if fm/1000.0>(self.freememory+1.5):
                    return mypickle().loadpickleonly('input')
                else:
                    print "not enough free memory on node"
                    sys.exit(0)

    def get_ppdock_plotdata(self):
        ara=np.zeros((30,len(self.scorerlist),2))
        for i in range(len(self.scorerlist)):
            ts=self.scorerlist[i][0]
            vs=self.scorerlist[i][1]
            k=0
            for so in [ts,vs]:
                if so.loadfilter!=None:
                    so.filter(ts.loadfilter)
                so.loadfilter=None
                par=self.rdlist[i]['bestpar']
                tnr=range(0,30)
                print so.assess(par)
                for p in tnr:
                    tcp=int(10**(p/10.0))
                    ara[p,i,k]=so.assess(par,slevel='top'+str(tcp)+'__nonempty_rmsd10ORirmsd4FIRST')
                k+=1
        return ara

    def generate_potential(self):
        so=scorer(model=self.model)
        if not hasattr(self,'finaloptpar'):
            self.load_fromlogdir(self.logpath)
        print so.assess(self.finaloptpar)
        print so.assess_model()
        rp=so.write_potential(filetype='lib')
        print rp
        libpath=rp[0]+'.opt.lib'
        print libpath
        print os.system('scp '+libpath+' '+runenv.jobserver+':'+os.path.join(runenv.serverUserPath+'lib',self.rundir+'.lib'))
