"""
   SOAP Optimizer for optimizing parameters 

"""

from env import *

from scorer import scorer
from recoveryFunction import sf
from mcmc import *

class optimizer(object):
    """
    Sampling and optimizing parameters for a scorer.
    """
    
    def __init__(self,scorer=[],seed=0):
        print 'Optimizer'
        self.inputcode=''
        self.scorer=scorer
        self.seed=seed
        self.searchpos=[]
        self.bestpar=[]
        self.currentpar=[]#parameter for next optimization step
        self.currentresult=-100000000
        self.bestresult=-100000000000000
        self.oldbestresult=-100000000000000
        self.improved=2
        self.testscorer=[]
        self.scorerid=0
        self.jobid=0
        self.runsperscorer=20
        self.quit=False
        self.indexlen=len(self.scorer.dsss.ds.sa)
        self.singlepar=False
        self.starttime=time.time()
        #self.partypelist=[] #the type of the par
        self.parindex=[0] # the index of the specific par in self.initialvalue/parvalues
        self.parproposemethod=[] # the proposemethod for the par
        self.parproposepars=[]
        self.parsearchindex=[] # the index of the corresponding search object
        self.globalresultlist=[]
        self.optstep=2
        self.stochastics=[]
        self.seed=0
        #self.optmethod=self.scorer['optmethod']
        self.currentmcmcm=None
        self.usesamemcmodel=False#whether build new model or use the same model
        self.mc_inner_converge={'improved':2,'maxloop':99999,'minloop':1}
        self.mc_outer_converge={'improved':4,'maxloop':99999}
        self.globalbestresult=-9999999
        self.oldglobalbestresult=-9999999
        self.globalbestpar=[]
        self.optm=None
         
    def optimize(self):
        self.optmind=0
        self.optm=None
        sp=self.start_from_existing_runs()# start from where it failed.
        print "starting from "+str(sp)
        for self.optmind in range(sp,len(self.scorer.model['sml'])):
            optm=self.scorer.model['sml'][self.optmind ]
            if self.optm!=None and optm['blockupdate']!=self.optm['blockupdate']:
                print "reset-model"
                self.currentmcmcm=None
            self.optm=optm
            if not os.path.isdir(str(self.scorerid)):
                try:
                    os.mkdir(str(self.scorerid))
                except:
                    pass
            if self.seed==0:
                print '///running'+str(self.scorerid)+'-'+str(self.optmind)
                print os.system('rm -f ./'+str(self.scorerid)+'/running-'+str(self.scorerid)+'-*')
                print os.system('touch ./'+str(self.scorerid)+'/running-'+str(self.scorerid)+'-'+str(self.optmind))
            func=optm['sm']
            if func.startswith('pbp'): #partbypart
                result=self.optimize_partbypart(func[4:])
            else:
                result=self.optimize_value(self.initialvalue,func)
            if func.startswith('p') or func.startswith('m'):
                self.initialvalue[...]=self.bestpar
            else:
                result=list(result)
                self.initialvalue[...]=result[0]
        return self.analyze_optimization_results(result)

    def start_from_existing_runs(self):
        fl=os.listdir('./')
        if not str(self.scorerid) in fl:
            return 0
        afl=os.listdir('./'+str(self.scorerid)+'/')
        runprogreesfile=[f for f in afl if f.startswith('running-'+str(self.scorerid)+'-')]
        if len(runprogreesfile)==0:
            return 0
        else:
            rpfl=[[f, int(f.split('-')[-1])] for f in runprogreesfile]
            rpfl.sort(key=itemgetter(1))
            return rpfl[-1][-1]

    def optimize_partbypart(self,method):
        #break the whole parameter set into small ones
        pl=len(self.scorer.searchlistspos)-1
        parvalues=np.copy(self.initialvalue)
        self.scorer.parvalues[...]=parvalues
        if pl==1:
            return self.optimize_value(self.initialvalue,method)
        for i in range(pl):
            pos=range(self.scorer.searchlistspos[i],self.scorer.searchlistspos[i+1])
            res=self.optimize_value(parvalues,method,pos)
            parvalues[self.scorer.searchlistspos[i]:self.scorer.searchlistspos[i+1]]=res[0]
            self.scorer.parvalues[...]=parvalues
        return [parvalues,res[1]]
    
    def optimize_value(self,initialvalue,method='powell',pos=[]):
        if len(pos)>0:
            self.searchpos=pos
            initialvalue=convert2ndarray(initialvalue)
            initialvalue=initialvalue[pos]
        else:
            self.searchpos=range(self.scorer.searchlistspos[-1])
        if method=='powell':
            res=scipy.optimize.fmin_powell(self.calc_score_point,initialvalue,full_output=True,maxfun=1000)
        elif method=='fmin':
            res=scipy.optimize.fmin(self.calc_score_point,initialvalue,full_output=True)
        elif method=='fmins':
            res=self.fmins(initialvalue)
        elif method=='fmins2':
            res=self.fmins(initialvalue,allpowell=True)
        elif method=='cg':
            res=scipy.optimize.fmin_cg(self.calc_score_point,initialvalue,full_output=True)
        elif method=='bfgs':
            res=scipy.optimize.fmin_bfgs(self.calc_score_point,initialvalue,full_output=True)
        elif method=='mcsa':
            res=self.mcsa(initialvalue)
        elif method.startswith('pd'):
            if len(method)==2:
                self.std=1
            elif method.startswith('pdre'):
                if self.seed>=len(runenv.mcmc_temperature_ratios1):
                    self.std=float(method[4:])
                else:
                    self.std=float(method[4:])*runenv.mcmc_temperature_ratios1[self.seed]
            else:
                self.std=float(method[2:])
            res=self.mcmc()
        elif method.startswith('pc'):#no exchange at first
            if len(method)==2:
                self.std=1
            elif method.startswith('pcre'):
                if self.seed>=len(runenv.mcmc_temperature_ratios1):
                    self.std=float(method[4:])
                else:
                    self.std=float(method[4:])*runenv.mcmc_temperature_ratios1[self.seed]
            else:
                self.std=float(method[2:])
            res=self.mcmc_noexchange_beforeconverge()
        elif method.startswith('pb'):#exchange based on probabilites
            if len(method)==2:
                self.std=1
            elif method.startswith('pbre'):
                if self.seed>=len(runenv.mcmc_temperature_ratios1):
                    self.std=float(method[4:])
                else:
                    self.std=float(method[4:])*runenv.mcmc_temperature_ratios1[self.seed]
            else:
                self.std=float(method[2:])
            res=self.mcmc_re()
        elif method.startswith('mccdsaf'):
            if len(method)==3:
                self.std=1
            else:
                self.std=float(method[3:])
            res=self.mcmc_continuous()
        elif method.startswith('mcae'):
            if len(method)==4:
                self.std=1
            else:
                self.std=float(method[4:])*runenv.mcmc_temperature_ratios1[self.seed]
            res=self.mcmc_tillconverge()
        elif method.startswith('mca'):
            if len(method)>3:
                self.std=float(method[3:])
            self.setup_mcpars()
            res=self.mcmc_tillconverge()
        elif method.startswith('mcs'):
            if len(method)>3:
                self.std=float(method[3:])
            self.setup_mcpars()
            res=self.mcmc_tillconverge_sametime()            
        elif method.startswith('mcp'):
            if len(method)>3:
                self.std=float(method[3:])
            self.setup_mcpars()
            res=self.mcmc_tillconverge_pt()           
        return res

    def fmins(self,initialvalue, allpowell=False):
        scorefunc=self.calc_score_point
        if allpowell:
            searchfunc=[scipy.optimize.fmin_powell,scipy.optimize.fmin_powell,scipy.optimize.fmin_powell,
                        scipy.optimize.fmin_powell,scipy.optimize.fmin_powell]
        else:
            searchfunc=[scipy.optimize.fmin,scipy.optimize.fmin_powell,scipy.optimize.fmin,
                        scipy.optimize.fmin_powell,scipy.optimize.fmin_powell]
        xtollist=[0.1,0.01,0.001,0.001,0.0001]
        k=0
        for func in searchfunc:
            res=func(scorefunc,initialvalue,full_output=True,xtol=xtollist[k])
            k=k+1
            initialvalue=res[0]
        if len(res[0].shape)==0:
            res=list(res)
            res[0]=[res[0]]
        return res
        
    def get_initial_value(self):
        seed=self.seed
        self.initialvalue=[]
        for search in self.scorer.model['searches']:
            id=search['InitialGenerator']
            pos=search['pos']
            #default propose method
            self.parindex.append(self.parindex[-1]+len(pos))
            if 'sftype' in search['object'] and search['object']['sftype'].startswith(('pl','sinb','normb','lognb','expb','rexpb','svdb','sacf')):
                self.parproposemethod.append('normal')
            else:
                self.parproposemethod.append('pnormal')
            self.parproposepars.append([])
            if 'sftype' in search['object'] and search['object']['sftype'].startswith('s') and search['key']=='par':
                sfo=sf(**search['object'])
                par=search['object']['par']
                self.initialvalue=self.initialvalue+list(par)
                self.parproposemethod[-1]='splineparn'
                self.parproposepars[-1]=([sfo.bins[0][0],sfo.bins[0][-1], len(par),(sfo.bins[0][-1]-sfo.bins[0][0])/(len(par)+1),sfo.bins[0][1]-sfo.bins[0][0]])
                continue              
            if isinstance(id,list):
                if isinstance(id[0],list):
                    self.initialvalue=self.initialvalue+id[seed]
                else:
                    self.initialvalue=self.initialvalue+id
            elif id['type']=='dfire':
                if search['object']['sftype'] in ['dfire','ig']:
                    self.initialvalue=self.initialvalue+[id['values'][seed]]
                elif search['object']['sftype'].startswith('s') and search['key']=='parvalue':
                    sfo=sf(**search['object'])
                    par=search['object']['par']
                    par=[sfo.bins[0][0]]+par+[sfo.bins[0][-1]]
                    dfirevalue=((np.array(par)/float(par[-1]))**(id['values'][seed]))
                    #how to generate the splines
                    ptype=search['object']['sftype']
                    if ptype[4] in ['l','s']:
                        dfirevalue=list(np.log(dfirevalue))
                        self.parproposemethod[-1]='normal'
                    #how to convert parvalue to control point values
                    if ptype[3]=='c':
                        dfirevaluediff=np.array([dfirevalue[0]]+list(dfirevalue[1:]-dfirevalue[0:-1]))
                        dfirevalue=list((dfirevaluediff/dfirevaluediff[-1])[:-1])
                    if ptype[1]=='a':
                        self.initialvalue=self.initialvalue+list(dfirevalue)
                    else:
                        self.initialvalue=self.initialvalue+list(dfirevalue[:-1])
                  
                elif search['object']['sftype'] in 'psn4' and search['key']=='parvalue':
                    sfo=sf(**search['object'])
                    par=search['object']['par']
                    par=list(np.linspace(sfo.bins[0][0], sfo.bins[0][-1],len(par)+1))
                    dfirevalue=((np.array(par)/float(par[-1]))**(id['values'][seed]))
                    dfirevaluediff=[dfirevalue[0]]+list(dfirevalue[1:]-dfirevalue[0:-1])
                    parvalue=dfirevaluediff
                    self.initialvalue=self.initialvalue+parvalue 
                elif search['object']['sftype'] in 'psn3' and search['key']=='parvalue':
                    sfo=sf(model=search['object'])
                    par=search['object']['par']
                    par=list(np.linspace(sfo.bins[0][0], sfo.bins[0][-1],len(par)+1))
                    dfirevalue=((np.array(par)/float(par[-1]))**(id['values'][seed]))
                    dfirevaluediff=np.array([dfirevalue[0]]+list(dfirevalue[1:]-dfirevalue[0:-1]))
                    parvalue=list((dfirevaluediff/dfirevaluediff[-1])[:-1])
                    self.initialvalue=self.initialvalue+parvalue 
                elif search['object']['sftype'] in 'psn2' and search['key']=='parvalue':
                    sfo=sf(model=search['object'])
                    par=search['object']['par']
                    par=list(np.linspace(sfo.bins[0][0], sfo.bins[0][-1],len(par)+1))
                    parvalue=list(((np.array(par)/float(par[-1]))**(id['values'][seed])))
                    self.initialvalue=self.initialvalue+parvalue                  
                elif search['object']['sftype']=='psn' and search['key']=='parvalue':
                    sfo=sf(model=search['object'])
                    par=search['object']['par']
                    par=list(np.linspace(sfo.bins[0][0], sfo.bins[0][-1],len(par)+1))
                    parvalue=list(((np.array(par)/float(par[-1]))**(id['values'][seed])))
                    self.initialvalue=self.initialvalue+parvalue[:-1]
                elif search['object']['sftype'].startswith('psn') and search['key']=='par':
                    sfo=sf(model=search['object'])
                    par=search['object']['par']
                    rpar=np.linspace(sfo.bins[0][0], sfo.bins[0][-1],len(par)+1)
                    self.initialvalue=self.initialvalue+list(rpar[:-1])
                    self.parproposemethod[-1]=('splineparn')
                    self.parproposepars[-1]=([sfo.bins[0][0],sfo.bins[0][-1], len(par),rpar[-1]-rpar[-2],sfo.bins[0][1]-sfo.bins[0][0]])
                else:
                    par=np.array(search['object']['par'])
                    print str(seed)+'dfire'+str(id['values'][seed])
                    self.initialvalue=self.initialvalue+list(((par/float(par[-1]))**(id['values'][seed]))[pos])
            elif id['type']=='sin':
                par=np.array(search['object']['par'])
                print str(seed)+'sin'+str(id['values'][seed])
                self.initialvalue=self.initialvalue+list(((np.abs(np.sin(par*3.1415926/180)+0.000000001))**(id['values'][seed]))[pos])            
            elif id['type']=='values':
                if not isinstance(id['values'][0],list):
                    id['values']=[[item] for item in id['values']]
                self.initialvalue=self.initialvalue+id['values'][seed]
            elif id['type']=='random':
                if search['key']=='par':
                    sfo=sf(**search['object'])
                    par=search['object']['par']
                    self.initialvalue=self.initialvalue+list(par)
                    self.parproposemethod[-1]=('splineparn')
                    self.parproposepars[-1]=([sfo.bins[0][0],sfo.bins[0][-1], len(par),(sfo.bins[0][-1]-sfo.bins[0][0])/(len(par)+1),sfo.bins[0][1]-sfo.bins[0][0]])                    
                else:
                    parlen=len(search['object'][search['key']])
                    print parlen
                    if 1:#stupid bug in scipy on the shared cluster, because the version is too old
                        rvl=[]
                        for i in range(parlen):
                            srv=scipy.stats.uniform.rvs(0,2,1)
                            try:
                                rvl.append(srv[0])
                            except:
                                rvl.append(srv)                        
                        self.initialvalue=self.initialvalue+rvl
                    #self.initialvalue=self.initialvalue+scipy.stats.uniform.rvs(0,2,parlen).tolist()
            elif id['type'].startswith('sf') and search['object']['type']=='sf':
                initialtype=id['type'][2:]
                sfo=sf(**search['object'])
                sfo.gen_initialvalue(id['values'][seed])
            elif id['type']=='ig':
                xi=np.arange(0.5,30,1)
                aprpr=[   144.434703,   1412.820011,   3834.829326,   7202.871092,
                 11341.975029,  16096.371468,  21329.303897,  26919.548441,
                 32742.984165,  38692.779078,  44664.067503,  50566.74941 ,
                 56312.288375,  61810.089186,  66974.927745,  71737.414983,
                 76029.089773,  79803.212568,  83024.246634,  85671.647628,
                 87740.847002,  89233.801466,  90173.44354 ,  90561.097026,
                 90429.358327,  89795.74556 ,  88703.604215,  87182.984942,
                 85264.625698,  83004.011034]
                us=scipy.interpolate.InterpolatedUnivariateSpline(xi*id['values'][seed],aprpr)
                ref=us(par)[pos]
                self.initialvalue=self.initialvalue+list(ref)
        self.initialvalue=np.array(self.initialvalue)
        if len(self.scorer.initialvalues)>0:
            boolfilter=self.scorer.initialvalues[-1,:]!=-9999
            self.initialvalue[boolfilter]=self.scorer.initialvalues[self.seed,:][boolfilter]
        if len(self.initialvalue)==1:
            self.singlepar=True
        if self.initialvalue.min()<0.0000000001:
            self.initialvalue+=0.0000000001 #pymc has problems with values at zero, which I think mess up the step method...
        self.bestpar=np.copy(self.initialvalue)
        self.currentpar=np.copy(self.bestpar)
        print "initialvalue performance "+str(self.calc_score_point(self.bestpar)) 
                
    def calc_score_point(self,parvalues):
        #print parvalues
        if self.singlepar and not (isinstance(parvalues,np.ndarray) and len(parvalues.shape)==1):
            parvalues=[parvalues]
        perc=self.scorer.assess(parvalues)#,self.searchpos
        #print 'Assessing: '+str(values)
        print perc
        return -perc

    def initialize_bias_array(self):
        self.biasnbin=10000
        self.biasna=np.zeros(self.biasnbin)
        self.biasbin=(self.idealresult-self.worstresult)/self.biasnbin
        self.biasw0=self.biasbin*2
        self.biasgamma=1.5
        self.biassigma=self.biasbin*1

    def get_bias(self,score):
        index=np.floor((score-self.worstresult)/self.biasbin)
        index=max(0,min(self.biasnbin-1,index))
        return self.biasna[index]        
                
    def update_bias(self, score):
        if(score < self.worstresult or score >self.idealresult):
            return 0
        vbias=self.get_bias(score)
        ww=self.biasw0*np.exp(-vbias/((2*self.likehoodstd**2)*(self.biasgamma-1.0)))
        i0=int(np.floor((score-4.0*self.biassigma-self.worstresult)/self.biasbin))
        i1=int(np.floor((score+4.0*self.biassigma-self.worstresult)/self.biasbin)+1)
        for i in range(max(0,i0),min(i1,self.biasnbin-1)):
            xx=self.worstresult + float(i)*self.biasbin
            dp=(xx-score)/self.biassigma
            self.biasna[i] = self.biasna[i]+ww*np.exp(-0.5*dp*dp)
        return vbias

    def calc_score_point_block(self,parvalues,blockind):
        #print parvalues
        if self.singlepar and not (isinstance(parvalues,np.ndarray) and len(parvalues.shape)==1):
            parvalues=[parvalues]
        perc=self.scorer.assess_block(parvalues,blockind)#,self.searchpos
        #print 'Assessing: '+str(values)
        print perc
        return -perc

    def runtask_cluster(self,inputcode):
        #the individual nodes on the cluster will run this function, if other funtions needs to be run here, chage the code here.
        self.runtask_local(inputcode)
        res.append(self.optlist2array(self.scoredict.values()))
        fh=open('optimizer'+self.inputcode+'.pickle','wb')
        pickle.dump(res,fh)
        fh.close()
        return 0
    
    def runtask_local(self,inputcode=5,testscorer=[]):
        self.testscorer=testscorer
        #pdb.set_trace()
        self.idealresult=self.scorer.assess_ideal()
        self.worstresult=self.scorer.assess_worst()
        print "inputcode"
        print inputcode
        print os.system('pwd')
        self.path='./'
        self.inputcode=inputcode
        self.seed=int(inputcode)-1
        self.get_initial_value()
        print self.initialvalue
        res=self.optimize()
        return res  

    def get_MCMC_single_stochastic_block(self,initialvalue,i):
        @stochastic(observed=False,name='a'+str(i))
        def a(value=initialvalue):
            if (value>=-1000).all() and (value<1000).all():
                return -1
            else:
                #print "Not accepatable value"
                return -np.Inf
        return a

    def get_MCMC_model_withblocks(self,initialvalue=[]):
        parlist=[]
        for i in range(len(self.scorer.searchblocklist)):
            block=self.scorer.searchblocklist[i][-1]
            a=self.get_MCMC_single_stochastic_block(initialvalue[block[0]:block[1]],i)
            parlist.append(a)
        
        print parlist
        print runenv.stepblockind
        print parlist[runenv.stepblockind]
        self.stochasticlist=parlist
        @deterministic()
        def p(v=parlist):
            return -self.calc_score_point_block(v[runenv.stepblockind],runenv.stepblockind)
            
        @stochastic(observed=True)
        def D(value=[],par=p):
            perc=par
            if self.optm['add_bias']:
                bias=self.update_bias(perc)
                print bias
                perc=perc-bias
            rvalue=-(perc-self.idealresult)**2/(2*self.likehoodstd**2)            
            return rvalue
        m=MyMCMC(parlist+[p,D])
        return m

    def get_MCMC_model_singleblock(self,initialvalue,blockind): #not implemented properly yet
        @stochastic(observed=False)
        def a(value=initialvalue):
            if (value>=-1000).all() and (value<1000).all():
                return -1
            else:
                #print "Not accepatable value"
                return -np.Inf
        
        self.stochasticlist.append(a)
        @deterministic()
        def p(v=a):            
            return -self.calc_score_point_block(v,blockind)
            
        @stochastic(observed=True)
        def D(value=[],par=p):
            perc=par
            rvalue=-(perc-self.idealresult)**2/(2*self.likehoodstd**2)            
            return rvalue
        
        m=MyMCMC(parlist+[p,D])
        return m
                
    def get_MCMC_model(self,initialvalue=[]):
        @stochastic(observed=False)
        def a(value=initialvalue):
            if (value>=-1000).all() and (value<1000).all():
                return -1
            else:
                #print "Not accepatable value"
                return -np.Inf
        
        self.stochasticlist=[a]
        @deterministic()
        def p(v=a):            
            return -self.calc_score_point(v)
    
        @stochastic(observed=True)
        def D(value=[],par=p):
            perc=par
            if self.optm['add_bias']:
                bias=self.update_bias(perc)
                print bias
                perc=perc-bias
            rvalue=-(perc-self.idealresult)**2/(2*self.likehoodstd**2)            
            return rvalue
        
        m=MyMCMC([a,p,D])
        return m
    
    def save_currentbest(self):
        os.system('rm -f ./'+str(self.scorerid)+'/'+'B'+str(self.optmind)+'est'+str(self.scorerid)+'#'+str(self.jobid)+'#*')
        fh=open('./'+str(self.scorerid)+'/B'+str(self.optmind)+'est'+str(self.scorerid)+'#'+str(self.jobid)+'#'+str(self.bestresult)+'#'+str(self.likehoodstd)+'#.pickle','w')
        pickle.dump([self.currentpar,self.currentresult,self.bestpar,self.bestresult],fh)
        fh.flush()
        fh.close()
        print 'B'+str(self.optmind)+'est'+str(self.scorerid)+'#'+str(self.jobid)+'#'+str(self.bestresult)+'#'+str(self.likehoodstd)+'#.pickle'
        print "saving result finished"

    def get_filelist(self): # collect result from otherruns, update globalbestresult/par
        print "-----------get_filelist--------------"
        self.allfilelist=os.listdir('./'+str(self.scorerid)+'/')
        #print os.getcwd()
        #get the global best result
        fl=[item for item in self.allfilelist if item.startswith('B'+str(self.optmind)+'est'+str(self.scorerid)+'#')]
        random.shuffle(fl)
        self.globalresultlist=[]
        bestresult=-9999999999999
        if len(fl)>0:
            for f in fl:
                fss=f.split('#')
                cr=float(fss[2])
                self.globalresultlist.append([float(fss[3]),cr,fss[0],int(fss[1]),f])
                if cr>bestresult:
                    bestresult=cr
                    bestf=f
            self.globalresultlist.sort(key=itemgetter(3))
            ngrl=[]
            for sgr in self.globalresultlist:
                if len(ngrl)==0:
                    ngrl.append(sgr)
                else:
                    if sgr[3]==ngrl[-1][3]:#from the same run
                        if sgr[1]>ngrl[-1][1]:
                            ngrl[-1]=sgr
                    else:
                        ngrl.append(sgr)
            self.globalresultlist=ngrl
            print 'Global Bestresult '+bestf+'   '+str(bestresult)
            time.sleep(0.1)
            try:
                fh=open('./'+str(self.scorerid)+'/'+bestf)
                bp=pickle.load(fh)
                fh.close()
            except:
                self.allfilelist=os.listdir('./'+str(self.scorerid)+'/')
                bestf=[f for f in self.allfilelist if f.startswith(bestf.split('#')[0]+'#')][0]
                fh=open('./'+str(self.scorerid)+'/'+bestf)
                bp=pickle.load(fh)
                fh.close()                
            self.globalbestresult=bp[-1]
            self.globalbestpar=bp[-2]
            print 'Global Bestresult '+bestf+'   '+str(self.globalbestresult)
        print "-----------get_filelist  END----------"

    def check_globalresult(self):
        print "////////////check_globalresult/////////"
        if self.globalbestresult>self.bestresult:
            self.improved=self.mc_outer_converge['improved']
            print "Using global best result "+str(self.globalbestresult)
            self.bestpar=self.globalbestpar
            self.bestresult=self.globalbestresult
        print "////////////check_globalresult END/////////"
            
    def set_globalbest2current(self):
        if 'using_globalbest' in self.optm and not self.optm['using_globalbest']:
            return 0
        self.get_filelist()
        self.check_globalresult()
        self.currentpar=self.bestpar
        self.currentresult=self.bestresult
                  
    def run_takewaytoolong(self):
        bestresult=self.globalbestresult
        if self.bestresult>bestresult+abs(bestresult)*0.001:
            print "Will conitune because current result is much better than global best"+str(self.bestresult)+' '+str(bestresult)
            return False
        tl=[float(item.split('#')[-1]) for item in self.allfilelist if item.startswith('Finished'+str(self.scorerid)+'#')]
        if len(tl) < self.runsperscorer*0.7999:
            print "Will continue because less than 80% runs has finished"
            return False
        elif len(tl) > self.runsperscorer:
            raise Exception("Fatal error, bugs in code, the runsperscorer variable is wrong")
        ta=np.array(tl)
        tastd=np.std(ta)
        tamean=ta.mean()
        crt=(time.time()-self.starttime)
        if crt>tamean*2 and crt>tamean+3*tastd:
            print "Runs will quit because it taks way too long "+str(crt)+" compared to "+str(tamean)
            self.quit=True
            return True
        else:
            print "will continue because it is not too long yet "+str(crt)+' mean '+str(tamean)+' std '+str(tastd)
            return False
        
    def mcmc(self,steps=2000): #repeated mcmc
        std=self.std
        if 'steps' in self.scorer.model:
            steps=self.scorer.model['steps']
        self.idealresult=self.scorer.assess_ideal()
        self.worstresult=self.scorer.assess_worst()
        self.likehoodstd=std*(self.idealresult)/(2*len(self.scorer.dsss.ds.indexlist))
        pl=self.mcsa_singleT(self.initialvalue, self.likehoodstd, steps,nonadaptive=True)
        for i in range(1):
            if self.quit:
                break
            pl=self.mcsa_singleT(self.bestpar, self.likehoodstd, steps,nonadaptive=True)
            if self.quit:
                break
        self.check_globalresult()
        self.save_currentbest()
        while self.improved>0 and (not self.quit):
            pl=self.mcsa_singleT(self.bestpar, self.likehoodstd, steps,nonadaptive=True)
            self.check_globalresult()
            self.save_currentbest()
        parl=[]
        rl=[]
        for i in range(0,steps):
            parl.append(list(m.trace('a')[i][0]))
            rl.append(pl[i])
        return [parl,rl]

    def setup_mcpars(self):
        if 'tune_interval' in self.optm:
            self.steps=self.optm['tune_interval']+1
        else:
            self.steps=201
        print "tune_interval "+str(self.steps)
        self.idealresult=self.scorer.assess_ideal()
        self.worstresult=self.scorer.assess_worst()
        if 'temperature_distribution' in self.optm:
            self.std=self.optm['temperature_distribution'][self.seed]
            self.exchange_times=np.zeros(self.runsperscorer)
            self.exchange_totaltimes=0
        self.likehoodstd=self.std*(self.idealresult)/(2*len(self.scorer.dsss.ds.indexlist))
        if ('sample_schedule' in self.optm) and isinstance(self.optm['sample_schedule'], dict):
            if 'inner' in self.optm['sample_schedule']:
                self.mc_inner_converge=self.optm['sample_schedule']['inner']
            if 'outer' in self.optm['sample_schedule']:
                self.mc_outer_converge=self.optm['sample_schedule']['outer']
        if 'add_bias' in self.optm and self.optm['add_bias']:
            self.initialize_bias_array()
        else:
            self.optm['add_bias']=False
        
    def setup_mcmcm(self):
        if self.currentmcmcm==None:
            if self.optm['blockupdate']:
                runenv.stepblockind=0
                m=self.get_MCMC_model_withblocks(self.currentpar)
                self.setup_stepmethod_allblocks(m)
            else:
                m=self.get_MCMC_model(self.currentpar)
                self.setup_stepmethod(m)
            self.currentmcmcm=m
            return m
        elif self.optm['stepmethod'].startswith('amp'):#adaptive
            if self.optm['blockupdate']:
                for i in range(len(self.stochasticlist)):
                    block=self.scorer.searchblocklist[i][-1]
                    sm=self.currentmcmcm.step_method_dict[self.stochasticlist[i]][0]
                    sorted(sm.stochastics)[0].value=np.copy(self.currentpar[block[0]:block[1]])
                    sm.accepted=0.0
                    sm.rejected=0.0
            else:
                sm=self.currentmcmcm.step_method_dict.values()[0][0]
                sorted(sm.stochastics)[0].value=np.copy(self.currentpar)
                sm.accepted=0.0
                sm.rejected=0.0
        elif self.optm['reset_betweenruns']==0:
            if self.optm['blockupdate']:
                for i in range(len(self.stochasticlist)):
                    block=self.scorer.searchblocklist[i][-1]
                    sm=self.currentmcmcm.step_method_dict[self.stochasticlist[i]][0]
                    sm.propose_sd=np.copy(self.currentpar[block[0]:block[1]])
                    sm.adaptive_scale_factor=1
                    sm.stochastics.value=np.copy(self.currentpar[block[0]:block[1]])
                    sm.accepted=0.0
                    sm.rejected=0.0
            else:
                sm=self.currentmcmcm.step_method_dict.values()[0][0]
                sm.propose_sd=np.copy(self.currentpar)
                sm.adaptive_scale_factor=1
                print "reset adaptive_scale_factor "+str(sm.adaptive_scale_factor)
                sm.stochastics.value=np.copy(self.currentpar)
                sm.accepted=0.0
                sm.rejected=0.0
        else:
            #reset stochastic.value, propose_sd, but not adaptive_scale_factor
            if self.optm['blockupdate']:
                for i in range(len(self.stochasticlist)):
                    block=self.scorer.searchblocklist[i][-1]
                    sm=self.currentmcmcm.step_method_dict[self.stochasticlist[i]][0]
                    sm.propose_sd=np.copy(self.currentpar[block[0]:block[1]])
                    sm.adaptive_scale_factor=sm.adaptive_scale_factor*self.optm['reset_betweenruns']
                    sm.stochastic.value=np.copy(self.currentpar[block[0]:block[1]])
                    sm.accepted=0.0
                    sm.rejected=0.0
            else:
                sm=self.currentmcmcm.step_method_dict.values()[0][0]
                sm.propose_sd=np.copy(self.currentpar)
                sm.adaptive_scale_factor=sm.adaptive_scale_factor*self.optm['reset_betweenruns']
                sm.stochastic.value=np.copy(self.currentpar)
                sm.accepted=0.0
                sm.rejected=0.0
        return self.currentmcmcm
            
    def mcmc_tillconverge(self):
        print "------> mcmc_tillconverge"
        #start sampling
        self.improved=self.mc_outer_converge['improved']
        pl=self.mcmc_converge_single()
        self.save_currentbest()
        k=0
        self.quit=False
        while True:
            self.set_globalbest2current()
            pl=self.mcmc_converge_single()
            self.save_currentbest()
            k=k+1
            if self.improved<=0 or self.quit or k>=self.mc_outer_converge['maxloop']:
                self.currentpar=self.bestpar
                self.likehoodstd=self.likehoodstd/self.std
                pl=self.mcmc_converge_single()
                break
        return pl

    def mcmc_converge_single(self):
        m=self.setup_mcmcm()     
        pl=self.mcmc_sample(m, self.steps,update_status=False)
        cbr=pl[1].max()
        improved=self.mc_inner_converge['improved']
        k=0
        while k<self.mc_inner_converge['minloop'] or (improved>0 and (not self.quit) and k<self.mc_inner_converge['maxloop']):
            #m.restore_sampler_state()
            pl=self.mcmc_sample(m, self.steps, update_status=False)
            k=k+1
            if pl[1].max()>cbr:
                cbr=pl[1].max()
                improved=self.mc_inner_converge['improved']
            else:
                improved-=1
            print "##########INNER########### improved "+str(improved)
            fl=os.listdir('./'+str(self.scorerid)+'/')
            if 'r'+str(self.optmind)+'quit' in fl:
                self.quit=True
                print "instructed to quit"
                break             
        self.update_runstatus(pl[0],pl[1])
        print '###################Single MCMC now converged'
        return pl

    def mcmc_tillconverge_sametime(self):
        print "------> mcmc_tillconverge_sametime"
        #start sampling
        #pl=self.mcmc_converge_single()
        #self.save_currentbest()
        self.improved=self.mc_outer_converge['improved']
        k=0
        self.quit=False
        while self.improved>0 and not (self.quit) and k<self.mc_outer_converge['maxloop']:
            if k>0:
                self.set_globalbest2current()
            pl=self.mcmc_converge_single_sametime()
            self.save_currentbest()
            k=k+1
        print os.system('touch ./'+str(self.scorerid)+'/'+'Finished'+str(self.optmind)+'-'+str(self.scorerid))
        return pl

    def mcmc_converge_single_sametime(self):
        m=self.setup_mcmcm()     
        #pl=self.mcmc_sample(m, self.steps,update_status=False, quittype=2)
        cbr=-999999
        improved=self.mc_inner_converge['improved']
        k=0
        while k<self.mc_inner_converge['minloop'] or (improved>0 and (not self.quit) and k<self.mc_inner_converge['maxloop']):
            #m.restore_sampler_state()
            pl=self.mcmc_sample(m, self.steps, update_status=False, quittype=2)
            k=k+1
            if pl[1].max()>cbr:
                cbr=pl[1].max()
                improved=self.mc_inner_converge['improved']
            else:
                improved-=1
            print "##########INNER########### improved "+str(improved)
            if k>=self.mc_inner_converge['minloop'] and ('Finished'+str(self.optmind)+'-'+str(self.scorerid) in self.allfilelist):
                print "breaking because other has finished"
                self.quit=True
                break            
        self.update_runstatus(pl[0],pl[1])
        print '###################Single MCMC now converged'
        return pl

    def mcmc_tillconverge_pt(self):
        print "------> mcmc_tillconverge_pt"
        self.improved=self.mc_outer_converge['improved']
        #self.mcmc_converge_single_pt(save=False)
        k=0
        self.quit=False
        while True:
            pl=self.mcmc_converge_single_pt()
            k=k+1
            self.wait_pt()
            if self.quit  or k>=self.mc_outer_converge['maxloop']:
                #self.currentpar=self.bestpar
                #self.likehoodstd=self.likehoodstd/self.std
                #pl=self.mcmc_converge_single()
                break
        print self.optm['temperature_distribution']
        print self.exchange_times/self.exchange_totaltimes
        return pl
    
    def mcmc_converge_single_pt(self,save=True):
        m=self.setup_mcmcm()             
        #pl=self.mcmc_sample(m, self.steps,update_status=False,quittype=1)
        cbr=-99999
        improved=self.mc_inner_converge['improved']
        k=0
        while k<self.mc_inner_converge['minloop'] or (improved>0 and (not self.quit) and k<self.mc_inner_converge['maxloop']):
            #m.restore_sampler_state()
            pl=self.mcmc_sample(m, self.steps, update_status=False,quittype=1)
            k=k+1
            if pl[1].max()>cbr:
                cbr=pl[1].max()
                improved=self.mc_inner_converge['improved']
            else:
                improved-=1
            if pl[1].max()>self.bestresult:
                mi=pl[1].argmax()
                self.bestpar=pl[0][mi]
                self.bestresult=pl[1].max()
            print "##########INNER########### improved "+str(improved)
            fl=[item for item in self.allfilelist if item.startswith('B'+str(self.optmind)+'est'+str(self.scorerid)+'#')]
            if len(fl)>0:
                break
        print "starting replica exchange"
        self.currentpar=pl[0][-1]
        self.currentresult=pl[1][-1]
        if save:
            self.save_currentbest()
        return True

    def wait_pt(self):
        print "**********pt********"
        if self.seed==0:
            while True:
                self.allfilelist=os.listdir('./'+str(self.scorerid)+'/')
                fl=[item for item in self.allfilelist if item.startswith('B'+str(self.optmind)+'est'+str(self.scorerid)+'#')]
                if len(fl)>=self.runsperscorer:
                    break
                else:
                    print "waiting for others to report status on disk: "+str(fl)+' '+str(fl)
                    time.sleep(5)
                #enabling quit with mark
                if 's'+str(self.scorerid)+'#'+str(self.optmind)+'quitwithfinalround' in fl:
                    self.quit=True
                    return 
            time.sleep(1)
            self.get_filelist()
            print self.globalresultlist
            print "global best result "+str(self.oldglobalbestresult)+' > '+str(self.globalbestresult)
            if self.globalbestresult>self.oldglobalbestresult:
                self.improved=self.mc_outer_converge['improved']
            else:
                self.improved-=1
            self.oldglobalbestresult=self.globalbestresult
            print "*************** pt improved "+str(self.improved)
            if self.improved<=0:
                print os.system('touch ./'+str(self.scorerid)+'/s'+str(self.scorerid)+'#'+str(self.optmind)+'quitwithfinalround')
            else:
                #doing exchange
                self.pt_exhange()
        #loading result
        while True:
            lfn='p'+str(self.optmind)+'t'+str(self.scorerid)+'#'+str(self.jobid)+'.pickle'
            fl=os.listdir('./'+str(self.scorerid)+'/')
            if 's'+str(self.scorerid)+'#'+str(self.optmind)+'quitwithfinalround' in fl:
                self.quit=True
                break
            elif lfn in fl:
                print "currentresult "+str(self.currentresult)
                time.sleep(0.5)
                self.currentpar, self.currentresult=pickle.load(open('./'+str(self.scorerid)+'/'+lfn))
                print "exchange to result "+str(self.currentresult)
                print os.system('rm ./'+str(self.scorerid)+'/'+lfn)
                break
            else:
                print "waiting for the replica exchange status"
                if not 'B'+str(self.optmind)+'est'+str(self.scorerid)+'#'+str(self.jobid)+'#'+str(self.bestresult)+'#'+str(self.likehoodstd)+'#.pickle' in fl:
                    self.save_currentbest()
                print lfn
                time.sleep(5)
        #delete files
        if self.seed==0:
            print os.system('rm -f ./'+str(self.scorerid)+'/B'+str(self.optmind)+'est'+str(self.scorerid)+'#'+'*')
            print os.listdir('./'+str(self.scorerid)+'/')
        print "*********ptEND*******"
    
    def pt_exhange(self):
        self.exchange_totaltimes+=1.0
        grl=self.globalresultlist
        #loading actual result
        for sr in grl:
            sr.append(pickle.load(open('./'+str(self.scorerid)+'/'+sr[-1])))
        #starting exhange
        if self.optm['exchange_method']==0:
            il=range(len(self.globalresultlist)-1)
            il.reverse()
        elif self.optm['exchange_method']==1:
            il=range(1, len(self.globalresultlist)-1,2)
            self.optm['exchange_method']=2
        elif self.optm['exchange_method']==2:
            il=range(0, len(self.globalresultlist)-1,2)
            self.optm['exchange_method']=1
        for i in il:
            ri=self.globalresultlist[i]
            rj=self.globalresultlist[i+1]
            ti=2*(ri[0]**2)
            tj=2*(rj[0]**2)
            ei=(ri[-1][1]-self.idealresult)**2
            ej=(rj[-1][1]-self.idealresult)**2
            p=min(1,np.exp((ei-ej)*(1/ti-1/tj)))
            if np.random.rand()<p:
                print "doing exchange# "+str(i)+' '+str(ri[-1][1])+' '+str(rj[-1][1])+' '+str(p)
                self.exchange_times[i]+=1.0
                si=ri[-1]
                ri[-1]=rj[-1]
                rj[-1]=si
        #save exchange result
        for rl in self.globalresultlist:
            fh=open('./'+str(self.scorerid)+'/'+rl[2].replace('B'+str(self.optmind)+'est','p'+str(self.optmind)+'t')+'#'+str(rl[3])+'.pickle','w')
            pickle.dump(rl[-1][:2],fh)
            fh.flush()
            fh.close()
        #replica exhange finished...
    
    def mcmc_tillconverge_withblocks_singlebysingle(self):
        if 'tune_interval' in self.scorer.model:
            steps=self.scorer.model['tune_interval']+1
        else:
            steps=201
        self.idealresult=self.scorer.assess_ideal()
        self.worstresult=self.scorer.assess_worst()
        self.likehoodstd=self.std*(self.idealresult)/(2*len(self.scorer.dsss.ds.indexlist))
        #start sampling
        pl=self.mcmc_converge_single_allblocks(steps)
        pl=self.mcmc_converge_single_allblocks(steps)
        self.save_currentbest()
        self.check_globalresult()
        while self.improved>0 and (not self.quit):
            pl=self.mcmc_converge_single_allblocks(steps)
            self.save_currentbest()
            self.check_globalresult()
        return pl

    def mcmc_converge_single_loopthroughsingleblock(self,steps):
        rl=[]
        for blockind in range(len(self.scorer.searchblocklist)):
            pl=self.mcmc_converge_single_singleblock(steps,blockind)
            rl.append(pl)
        #check whether initial value is in the trace . Need to combine the traces
        self.update_runstatus(pl[0],pl[1])
        
    def mcmc_converge_single_singleblock(self,steps,blockind):
        print "####################Updating single block"
        self.optm['blockupdate']=False
        m=self.get_MCMC_model(self.currentpar,blockind)
        self.setup_stepmethod(m)        
        pl=self.mcmc_sample(m, steps,update_status=False,verbose=2)
        cbr=pl[1].max()
        improved=4
        while improved>0 and (not self.quit):
            m.restore_sampler_state()
            pl=self.mcmc_sample(m, steps, update_status=False)
            if pl[1].max()>cbr:
                cbr=pl[1].max()
                improved=4
            else:
                improved-=1
            print "##################### improved "+str(improved)
        print '###################Single MCMC now converged'
        return pl

    def global_exchange(self):
        grl=self.globalresultlist
        if len(grl)==0:
            return 0
        pa=np.zeros(len(grl))
        for i in range(len(grl)):
            pa[i]=np.exp(-(grl[i][1]-self.idealresult)**2/(2*self.likehoodstd**2))
        pa[...]=pa/pa.sum()
        pa=pa.cumsum()
        i=np.searchsorted(pa, np.random.rand())
        gb=pickle.load(open(grl[i][2]))
        gbp=gb[0]
        gbr=gb[1]
        if gbr>self.bestresult:
            if gbr>self.oldbestresult+abs(self.oldbestresult)*0.001:
                self.improved=2
                print "Global imporved more than 0.1% "+str(self.improved)
            else:
                if self.improved<2:
                    if self.stepimproved==0:
                        self.improved+=0.5
                    print "Global imporved less than 0.1% "+str(self.improved)
            print "Using global best result"
        self.bestpar=gbp
        self.bestresult=gbr

    def replica_exchange(self):
        grl=self.globalresultlist
        if len(grl)==0:
            return 0
        pa=np.zeros(len(grl))
        for i in range(len(grl)):
            pa[i]=np.exp(-(grl[i][1]+self.idealresult)**2/(2*self.likehoodstd**2))
        pa[...]=pa/pa.sum()
        pa=pa.cumsum()
        i=np.searchsorted(pa, np.random.rand())
        gb=pickle.load(open(grl[i][2]))
        gbp=gb[0]
        gbr=gb[1]
        if gbr>self.bestresult:
            if gbr>self.oldbestresult+abs(self.oldbestresult)*0.001:
                self.improved=2
                print "Global imporved more than 0.1% "+str(self.improved)
            else:
                if self.improved<2:
                    if self.stepimproved==0:
                        self.improved+=0.5
                    print "Global imporved less than 0.1% "+str(self.improved)
            print "Using global best result"
        self.bestpar=gbp
        self.bestresult=gbr

    def mcmc_re(self,steps=2000): #repeated mcmc
        std=self.std
        if 'steps' in self.scorer.model:
            steps=self.scorer.model['steps']
        self.idealresult=self.scorer.assess_ideal()
        self.worstresult=self.scorer.assess_worst()
        self.likehoodstd=std*(self.idealresult)/(2*len(self.scorer.dsss.ds.indexlist))
        pl=self.mcsa_singleT(self.initialvalue, self.likehoodstd, steps,nonadaptive=True)
        for i in range(1):
            if self.quit:
                break
            pl=self.mcsa_singleT(self.bestpar, self.likehoodstd, steps,nonadaptive=True)
            if self.quit:
                break
        self.global_exchange()
        self.save_currentbest()        
        while self.improved>0 and (not self.quit):
            pl=self.mcsa_singleT(self.bestpar, self.likehoodstd, steps,nonadaptive=True)
            self.global_exchange()
            self.save_currentbest()
        parl=[]
        rl=[]
        for i in range(0,steps):
            parl.append(list(m.trace('a')[i][0]))
            rl.append(pl[i])
        return [parl,rl]

    def mcmc_noexchange_beforeconverge(self,steps=2000): #repeated mcmc
        std=self.std
        if 'steps' in self.scorer.model:
            steps=self.scorer.model['steps']
        self.idealresult=self.scorer.assess_ideal()
        self.worstresult=self.scorer.assess_worst()
        self.likehoodstd=std*(self.idealresult)/(2*len(self.scorer.dsss.ds.indexlist))
        pl=self.mcsa_singleT(self.initialvalue, self.likehoodstd, steps,nonadaptive=True)
        while self.improved>0 and (not self.quit):
            pl=self.mcsa_singleT(self.bestpar, self.likehoodstd, steps,nonadaptive=True)
        #begin exchange
        self.check_globalresult()
        self.save_currentbest()
        while self.improved>0 and (not self.quit):
            pl=self.mcsa_singleT(self.bestpar, self.likehoodstd, steps,nonadaptive=True)
            self.check_globalresult()
            self.save_currentbest()
        parl=[]
        rl=[]
        for i in range(0,steps):
            parl.append(list(m.trace('a')[i][0]))
            rl.append(pl[i])
        return [parl,rl]

    def mcmc_continuous(self,steps=2000): #repeated mcmc with memory
        std=self.std
        if 'steps' in self.scorer.model:
            steps=self.scorer.model['steps']
        self.idealresult=self.scorer.assess_ideal()
        self.worstresult=self.scorer.assess_worst()
        self.likehoodstd=std*(self.idealresult)/(2*len(self.scorer.dsss.ds.indexlist))
        mdl=self.get_MCMC_model(self.initialvalue)
        self.setup_stepmethod(mdl)
        pl=self.mcmc_sample(mdl, steps)
        for i in range(1):
            if self.quit:
                break
            pl=self.mcmc_sample(mdl, steps)
            if self.quit:
                break
        self.save_currentbest()
        self.check_globalresult()
        while self.improved>0 and (not self.quit):
            pl=self.self.mcmc_sample(mdl, steps)
            self.save_currentbest()
            self.check_globalresult()
        parl=[]
        rl=[]
        for i in range(0,steps):
            parl.append(list(m.trace('a')[i][0]))
            rl.append(pl[i])
        return [parl,rl]
            
    def analyze_optimization_results(self, results):
        opl=results[0]
        orl=results[1]
        if not (isinstance(orl, list) or (isinstance(orl,np.ndarray) and len(orl.shape)>0)):
            orl=[orl]
            opl=[opl]
        parlist=[]
        pln=[]
        trl=[]
        rl=[]
        rrl=[]
        del self.scorer
        if isinstance(self.testscorer,scorer):
            testscorer=self.testscorer
            testscorer=load_scorer(testscorer,self.indexlen)
            testideal=testscorer.assess_ideal()
            testworst=testscorer.assess_worst()        
            if len(testscorer.dsss.ds.indexlist)>0:
                teststd=self.std*(testideal)/(2*len(testscorer.dsss.ds.indexlist))
            else:
                teststd=0        
            for i in range(len(opl)):
                ni=opl[i]
                cr=orl[i]
                ni=list(ni)
                nis=any2str(ni)
                if not nis in parlist:
                    parlist.append(nis)
                    pln.append(ni)
                    rl.append(cr)
                    trl.append(testscorer.assess(ni))
    #rrl.append(((trl[-1]-testworst)/(testideal-testworst))/((cr-self.worstresult)/(self.idealresult-self.worstresult)))
            ta=np.array(trl)
            print "Number of sample points: "+str(len(ta))
            pa=ta.mean() # np.exp(-(ta-self.idealresult)**2/(2*rstd**2))
            print 'pa: '+str(pa)
            #pdb.set_trace()
            return {'PredictiveDensity':pa
                ,'parlist':pln,'trainresultlist':rl,'testresultlist':trl,'std':teststd,
                'idealresult':testscorer.assess_ideal()}#,'generalizability':np.array(rrl).mean()}
        else:
            return {'PredictiveDensity':0,'parlist':opl,'trainresultlist':orl,'testresultlist':[],'std':0,'idealresult':0}

    def setup_stepmethod_allblocks(self,m):
        #if len(cov)>0:
        #    print 'using adaptive metropolis method'
        #    m.use_step_method(AdaptiveMetropolis,list(m.stochastics)[0],cov=cov,interval=200,delay=200,verbose=0)
        if 'stepmethod' in self.optm:
            sms=self.optm['stepmethod']
            if sms.startswith('mymp'):
                runenv.acc_rate_ratio=float(sms[4:])
                #m.use_step_method(MyMetropolis,list(m.stochastics)[0])
                for i in range(len(self.stochasticlist)):
                    m.use_step_method(MyMetropolis,self.stochasticlist[i],self.parindex,self.parproposemethod,self.parproposepars,blockind=i,searchind=self.scorer.searchblocklist[i][0])
            elif sms.startswith('mxmp'):
                runenv.acc_rate_ratio=float(sms[4:])
                #m.use_step_method(MyMetropolis,list(m.stochastics)[0])
                for i in range(len(self.stochasticlist)):
                    m.use_step_method(MxMetropolis,self.stochasticlist[i],self.parindex,self.parproposemethod,self.parproposepars,blockind=i,searchind=self.scorer.searchblocklist[i][0])
            elif sms.startswith('mzmp'):
                runenv.acc_rate_ratio=float(sms[4:])
                #m.use_step_method(MyMetropolis,list(m.stochastics)[0])
                for i in range(len(self.stochasticlist)):
                    m.use_step_method(MzMetropolis,self.stochasticlist[i],self.parindex,self.parproposemethod,self.parproposepars,blockind=i,searchind=self.scorer.searchblocklist[i][0])
            elif sms.startswith('mp'):
                for i in range(len(self.stochasticlist)):
                    m.use_step_method(MMetropolis,self.stochasticlist[i],blockind=i,verbose=1,searchind=self.scorer.searchblocklist[i][0])
            elif sms.startswith('amp'):
                for i in range(len(self.stochasticlist)):
                    m.use_step_method(AdaptiveMetropolis,self.stochasticlist[i],interval=200,delay=200,verbose=1)
        else:
            raise Exception('please specify step methos in you model')
            
    def setup_stepmethod(self,m,cov=[]):
        #if len(cov)>0:
        #    print 'using adaptive metropolis method'
        #    m.use_step_method(AdaptiveMetropolis,list(m.stochastics)[0],cov=cov,interval=200,delay=200,verbose=0)
        if 'stepmethod' in self.optm:
            sms=self.optm['stepmethod']
            if sms.startswith('mymp'):
                runenv.acc_rate_ratio=float(sms[4:])
                #m.use_step_method(MyMetropolis,list(m.stochastics)[0])
                m.use_step_method(MyMetropolis,list(m.stochastics)[0],self.parindex,self.parproposemethod,self.parproposepars)
            elif sms.startswith('mxmp'):
                runenv.acc_rate_ratio=float(sms[4:])
                #m.use_step_method(MyMetropolis,list(m.stochastics)[0])
                m.use_step_method(MxMetropolis,list(m.stochastics)[0],self.parindex,self.parproposemethod,self.parproposepars)
            elif sms.startswith('mzmp'):
                runenv.acc_rate_ratio=float(sms[4:])
                #m.use_step_method(MyMetropolis,list(m.stochastics)[0])
                m.use_step_method(MzMetropolis,list(m.stochastics)[0],self.parindex,self.parproposemethod,self.parproposepars)
            elif sms.startswith('mp'):
                m.use_step_method(Metropolis,list(m.stochastics)[0],verbose=1)
            elif sms.startswith('amp'):
                m.use_step_method(AdaptiveMetropolis,list(m.stochastics)[0],interval=200,delay=200,verbose=1)
        else:
            m.use_step_method(Metropolis,list(m.stochastics)[0])

    def mcmc_sample(self,m,step,update_status=False,quittype=0):
        # sample, save state, get result, check global result and whether runs take too long
        #pdb.set_trace()
        #m.restore_sampler_state()
        print time.asctime()
        tn=step/50
        if quittype==1:
            for i in range(tn):
                m.sample(iter=50,tune_interval=51,verbose=mcverbose)
                self.allfilelist=os.listdir('./'+str(self.scorerid)+'/')
                fl=[item for item in self.allfilelist if item.startswith('B'+str(self.optmind)+'est'+str(self.scorerid)+'#')]
                if len(fl)>0:
                    print "quit becasue other has finished for pt"
                    break
            m.tune()
        elif quittype==2:
            for i in range(tn):
                m.sample(iter=50,tune_interval=51,verbose=mcverbose)
                self.allfilelist=os.listdir('./'+str(self.scorerid)+'/')
                if ('Finished'+str(self.optmind)+'-'+str(self.scorerid) in self.allfilelist):
                    print "quit becasue other has finished"
                    break
            m.tune()
        else:
            m.sample(iter=step,tune_interval=step-1,verbose=mcverbose)
        #m.save_state()
        resultlist=self.get_perc_trace(m)
        if update_status:
            self.update_runstatus(resultlist[0],resultlist[1])
        print time.asctime()
        return resultlist        

    def mcsa_singleT(self,initialvalue,T,step,cov=[],nonadaptive=False):
        #build a new model and do ssampling
        self.likehoodstd=T
        m=self.get_MCMC_model(initialvalue)
        self.setup_stepmethod(m)
        return self.mcmc_sample(m,step)
    
    def get_perc_trace(self,m):
        if self.optm['blockupdate']:
            al=[]
            for s in self.stochasticlist:
                al.append(m.trace(s).gettrace())
            naa=np.hstack(al)
            nap=m.trace('p').gettrace()
        else:
            naa=m.trace(self.stochasticlist[0]).gettrace()
            nap=m.trace('p').gettrace()
        return [naa,nap]
        
    def update_runstatus(self,naa,nap):
        print "*********update_runstatus*********"
        mi=np.argmax(nap)
        bp=naa[mi]
        br=nap[mi]
        self.oldbestresult=self.bestresult
        if br>=self.bestresult:
            if br>self.bestresult:#*self.mc_outer_converge['improved_less'][0]:
                self.stepimproved=self.mc_outer_converge['improved']
                self.improved=self.mc_outer_converge['improved']
                print "Imporved more than 0.1% "+str(self.improved)
            #elif br>self.bestresult:
            #    self.stepimproved=1
            #    self.improved-=self.mc_outer_converge['improved_less'][1]
            #    print "Improved but not more than 0.1% "+str(self.improved)                
            else:
                self.stepimproved=0
                self.improved-=1
                print "Not improved "+str(self.improved)
            self.bestresult=br
            self.bestpar=bp
        else:
            self.stepimproved=0
            self.improved-=1
            print "Not improved "+str(self.improved)
        print 'Best result: '+str(self.bestresult)
        self.get_filelist()
        #self.roundruntime=time.time()-st
        self.run_takewaytoolong()
        self.currentpar=naa[-1]#self.bestpar
        self.currentresult=nap[-1]#self.bestresult
        print self.bestresult
        print self.bestpar
        print "*********update_runstatus END*********"
            
