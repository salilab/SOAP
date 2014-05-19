"""
   SOAP Scorer, caculating a score of a statistical potential(with model and parameter) for benchmark

"""
from env import *
from rankScores import *
from benchmark import dsscorestats
from refineScores import *
from loop import sprefine
        
class scorer(object):
    """
    Calculate benchmark score/value given a model
    """
    
    def __init__(self,model=[]):
        self.bmtype=model['bmtype']
        self.model=model
        self.originalmodel=copy.deepcopy(model)
        self.scorertype=[]#0 basescore, 1:search for ratio, 2: sf scorer
        self.scorerlist=[]
        self.basescore=[]
        self.score=[]
        self.searchlists=[]
        self.searchlistspos=[0]
        self.clean=True
        self.searchscorerlist=[[] for i in range(len(self.model['searches']))]# the index of the scorers affected by the searches        
        self.scorersearchlist=[[] for i in range(len(self.model['scorers']))]
        if len(self.bmtype['dslist'])>0:
            self.initialize_stats()        
            self.initialize_searches()        
            self.initialize_scores()
            self.loadfilter=None # filters to apply when loading data
            self.blockind=-999# no block update
            if 'filters' in self.bmtype:
                self.filter(self.bmtype['filters'])
            self.find_scorersearch_cluster()
            self.initialvalues=[]
            self.singlescorerlist=[]
        #pdb.set_trace()

    def filter(self,filters): #use only part of the decoys for benchmarking...
        fa=self.get_filters(filters)
        self.apply_filter(fa)

    def get_filters(self,filters):
        safilter=np.ones(len(self.dsss.ds.sa),dtype=np.bool)
        for filter in filters:
            if filter['criteria'].startswith('top') and 'model' in filter:
                tn=float(filter['criteria'][3:])
                cm={'scorers':filter['model']['scorers'],'searches':[],'bmtype':filter['model']['bmtype']}
                so=scorer(cm)
                currentscore=np.copy(so.basescore)
                del cm
                del so
                for spos in self.dsss.ds.indexlist:
                    posscore=currentscore[spos[0]:spos[1]]
                    if tn<1:
                        tn=(spos[1]-spos[0])*tn
                    tn=int(tn)
                    si=np.argsort(posscore)
                    cfa=np.zeros(spos[1]-spos[0],dtype=np.bool)
                    cfa[si[:tn]]=1
                    safilter[spos[0]:spos[1]]=np.logical_and(cfa, safilter[spos[0]:spos[1]])
            elif filter['criteria'].startswith('top') and 'parvalue' in filter:
                tn=float(filter['criteria'][3:])
                self.parvalues[...]=filter['parvalue']
                self.assign_values2model()
                self.get_model_dsscore()
                currentscore=self.score
                for spos in self.dsss.ds.indexlist:
                    posscore=currentscore[spos[0]:spos[1]]
                    if tn<1:
                        stn=(spos[1]-spos[0])*tn
                    stn=int(stn)
                    si=np.argsort(posscore)
                    cfa=np.zeros(spos[1]-spos[0],dtype=np.bool)
                    cfa[si[:stn]]=1
                    safilter[spos[0]:spos[1]]=np.logical_and(cfa, safilter[spos[0]:spos[1]])
            elif filter['criteria'].startswith('proplb'):
                tn=float(filter['criteria'][3:])
                propname=filter['propname']
                lb=filter['lb']
                for spos in self.dsss.ds.indexlist:
                    proparr=self.dsss.ds.sa[propname][spos[0]:spos[1]]
                    cfa=np.zeros(spos[1]-spos[0],dtype=np.bool)
                    cfa[proparr>lb]=1
                    safilter[spos[0]:spos[1]]=np.logical_and(cfa, safilter[spos[0]:spos[1]])                
            elif filter['criteria'].startswith('rep'):
                tn=int(filter['criteria'][3:])+1
                for spos in self.dsss.ds.indexlist:
                    ai=np.linspace(0,spos[1]-spos[0],tn).astype(np.int)[:-1]
                    cfa=np.zeros(spos[1]-spos[0],dtype=np.bool)
                    cfa[ai]=1
                    safilter[spos[0]:spos[1]]=np.logical_and(cfa, safilter[spos[0]:spos[1]])                
            elif filter['criteria']=='has':
                cd=filter['criteriadict']
                for spos in self.dsss.ds.indexlist:
                    sa=self.dsss.ds.sa[spos[0]:spos[1]]
                    keep=False
                    for key in cd:
                        if np.any(sa[key]<cd[key]):
                            keep=True
                            break
                    if not keep:
                        #print spos
                        safilter[spos[0]:spos[1]]=False
            elif filter['criteria']=='nocode':
                cl=filter['codelist']
                cl=[self.dsss.ds.dsname+c for c in cl]
                for spos,code in zip(self.dsss.ds.indexlist,self.dsss.ds.codelist):
                    if code in cl:
                        safilter[spos[0]:spos[1]]=False
                        print code
            else:
                raise Exception('filter not recoginzed'+str(filter))
        #pdb.set_trace()
        return safilter

    def apply_filter(self,fa):
        #use only part of the big array// this is the only function that changes the real size of the score object
        #******pdb.set_trace()
        self.dsss.ds.filter(fa)
        self.dsss.initialize_slevel(self.bmtype['criteria'])
        #pdb.set_trace()
        self.scorearray=np.array(self.scorearray[:,fa],order='C')
        self.score=self.score[fa]
        self.basescore=self.scorearray[-1,:]
        k=-1
        for i in range(len(self.scorerlist)):
            st=self.scorertype[i]            
            if st==0:
                continue
            else:
                k+=1
                if st==1:
                    self.scorerlist[i]=self.scorearray[k,:]
                elif st==2:
                    if self.model['scorers'][i]['type']=='sf':
                        self.scorerlist[i].filter(fa, self.scorearray[k,:])
                    elif self.model['scorers'][i]['type']=='scaledsp':
                        pass
                    
    def relink_scorearray(self):
        #relink the score arrays
        self.basescore=self.scorearray[-1,:]
        k=-1
        for i in range(len(self.scorerlist)):
            st=self.scorertype[i]            
            if st==0:
                continue
            else:
                k+=1
                if st==1:
                    self.scorerlist[i]=self.scorearray[k,:]
                elif st==2:
                    if self.model['scorers'][i]['type']=='sf':
                        self.scorerlist[i].score=self.scorearray[k,:]
                    elif self.model['scorers'][i]['type']=='scaledsp':
                        pass        

    def initialize_stats(self):
        #initialize the scorestats object that can be used to analyze score, and give the performance of the score
        if self.bmtype['type']=='dsscore':
            if 'dsss' in self.bmtype:
                self.dsss=self.bmtype['dsss']
            else:
                self.dsss=dsscorestats(dslist=self.bmtype['dslist'])
            self.dsss.initialize_slevel(self.bmtype['criteria'])
            self.basescore=np.zeros(len(self.dsss.ds.sa),dtype=np.float32)
            self.score=np.zeros(len(self.basescore),dtype=np.float32)
        elif self.bmtype['type']=='sprefine':
            bmt=copy.deepcopy(self.bmtype)
            del bmt['type']
            self.sprefiner=sprefine(**bmt)
            self.sprefiner.initialize_dslist()
            self.atomclasslib='$(LIB)/atmcls-mf.lib'
            self.potlib='$(LIB)/dist-mf.lib'
            self.sprefiner.refpot=['$(LIB)/atmcls-mf.lib','$(LIB)/dist-mf.lib']
            for k in range(0,len(self.model['scorers'])):
                if self.model['scorers'][k]['type']=='scaledsp':
                    self.atomclasslib=runenv.libdir+feature(self.model['scorers'][k]['features']).get_lib()
                    self.potlib=scaledsp(model=self.model['scorers'][k]).ssppath+'.opt.lib'
                    break
        elif self.bmtype['type']=='structurescore':
            pass

    def initialize_scores(self):
        #add the scores that is not in the search list, because those stores won't change during our search
        if self.bmtype['type']=='dsscore' and self.bmtype['combine']=='scoresum':
            k=0
            for i in range(len(self.model['scorers'])):
                scorer=self.model['scorers'][i]
                if not sameobjectIn2list(scorer.values(),self.searchlists):
                    if len(self.basescore)>0:
                        try:
                        #if 1:
                            #print 'adding '+str(k)+'\n'+str(scorer)
                            self.basescore+=scorer['ratio'][0]*self.get_dsscore(scorer)
                            k=k+1
                        except Exception,e:
                            traceback.print_exc()
                            pdb.set_trace()
                    else:
                        self.basescore=scorer['ratio'][0]*self.get_dsscore(scorer)
                    self.scorertype.append(0)
                    self.scorerlist.append([])
                else: # the scorer is varied during search process
                    kl,il=findsameobjectkeys(scorer,self.searchlists)
                    for j in il:
                        self.searchscorerlist[j].append(i)
                        self.scorersearchlist[i].append(j)
                    if len(kl)==1 and kl[0]=='ratio':
                        #we are only varying the ratio between different scores
                        self.scorerlist.append(self.get_dsscore(scorer))
                        self.scorertype.append(1)
                    elif scorer['type']=='sf':
                        sc=scorer
                        if 'bm' in sc:
                            bmt=sc['bm']
                            sfscorer=sfdsscore(dslist=self.bmtype['dslist'], **sc)
                        else:
                            bmt=self.bmtype['bm']
                            sfscorer=sfdsscore(dslist=self.bmtype['dslist'],bm=bmt, **sc)
                        self.scorerlist.append(sfscorer)
                        self.scorertype.append(2)
                    else:
                        raise Exception('Donot know how to hintialize scorers')
            td=0
            for st in self.scorertype:
                if st>0:
                    td+=1
            nna=np.zeros([td+1,len(self.basescore)],dtype=np.float32)#number of different score + base score for each decoy
            ci=-1
            for i in range(len(self.scorertype)):
                so=self.scorerlist[i]
                st=self.scorertype[i]
                if st==1:
                    ci+=1
                    nna[ci,:]=so[...]
                    self.scorerlist[i]=nna[ci,:]
                    del so
                elif st==2:
                    ci+=1
                    nna[ci,:]=so.score[...]
                    del so.score
                    so.score=nna[ci,:]
            nna[-1,:]=self.basescore[...]
            del self.basescore
            self.basescore=nna[-1,:]
            self.scorearray=nna
            self.ratioarray=np.ones(td+1,dtype=np.float32)
                         
    def initialize_searches(self):
        #Initialize the scores need to be searched...
        if not 'searches' in self.model:
            return 0
        #self.searchobjects=searches
        for search in self.model['searches']:
            sl=0
            if isinstance(search['key'],list):
                raise Exception('search[key] can not be list')
            else:
                key=search['key']
                pos=search['key']
                slist=search['object'][key]
                convert2float(slist)
                self.searchlists.append(slist) # reference to the original list object whose value we want to change in the search process
                sl+=len(slist)
                self.searchlistspos.append(sl+self.searchlistspos[-1]) # the index of the list values in self.parvalues
        self.parvalues=np.zeros(self.searchlistspos[-1])
 
    def find_scorersearch_cluster(self):
        self.searchblocklist=[]
        csi=0
        while csi<len(self.model['searches']):
            searches=[csi]
            scorers=self.searchscorerlist[csi]
            while True:
                changed=False
                for si in scorers:
                    ssil=self.scorersearchlist[si]
                    for searchind in ssil:
                        if not searchind in searches:
                            searches.append(searchind)
                            changed=True
                for si in searches:
                    ssil=self.searchscorerlist[si]
                    for scorerind in ssil:
                        if not scorerind in scorers:
                            scorers.append(scorerind)
                            changed=True
                if not changed:
                    break
            if not len(searches)==(max(searches)-min(searches)+1):
                pdb.set_trace()
                raise Exception('searching block is not continuous???????')
            self.searchblocklist.append([searches,scorers,[self.searchlistspos[min(searches)],self.searchlistspos[max(searches)+1]]])
            csi=max(searches)+1
                        
    def get_dsscore(self,scorer):
        if scorer['type']=='scaledsp':
            return self.get_spdsscore(scorer).get_score()
        elif scorer['type']=='sf':
            sc=copy.deepcopy(scorer)
            del sc['type']
            #pdb.set_trace()
            if 'ratio' in sc:
                del sc['ratio']
            if not 'bm' in sc:
                sc['bm']=self.bmtype['bm']
            score=sfdsscore(dslist=self.bmtype['dslist'],**sc).get_score()
            return score
            
    def get_spdsscore(self, scorer):
        sc=copy.deepcopy(scorer)        
        if 'bm' in sc:
            bmt=sc['bm']
        else:
            bmt=self.bmtype['bm']
        rsp=rawsp(model=sc)
        ssp=scaledsp(pm=scorer['pm'],rspo=rsp)
        sspscore=spdsscore(dslist=self.bmtype['dslist'],bm=bmt,ssp=ssp)
        return sspscore
       
    def get_sf(self,parlist=[]):
        raise Exception('scorer.get_sf Not working yet')
        sf=[]
        if not (isinstance(parlist,list) or isinstance(parlist,tuple)  or isinstance(parlist,np.ndarray)):
            parlist=[parlist]
        if len(parlist)==0:
            return []
        for i in range(0,len(self.searches)):
            search=self.searches[i]
            if search['object']['type']=='sf':
                scorer=search['scorer']
                if search['key']=='parvalue':
                    scorer.parvalue[search['pos']]=parlist[self.parpos[i]:self.parpos[i+1]]
                    sf.append(list(scorer.get_sf()))
        return sf

    def get_model_dsscore(self):
        if self.blockind>=0:
            sil=self.searchblocklist[self.blockind][1]
        k=-1
        poslist=self.dsss.ds.indexlist
        #for index in poslist:
        #    self.score[index[0]:index[1]]=self.basescore[index[0]:index[1]]
        for sm,st,so in zip(self.model['scorers'],self.scorertype,self.scorerlist):
            if st==0:#already included in the basescore
                continue
            else:
                k=k+1
                #if   st==1: #just need to apply the ratio to the score
                self.ratioarray[k]=sm['ratio'][0]
                if st==2: #need to recalculate the score
                    if sm['type']=='sf':
                        if self.blockind<0 or (k in sil): #no block or within the specific block
                            so.parvalue[:]=sm['parvalue']
                            if so.par!=None:
                                so.par[:]=sm['par']
                            so.get_score(indexlist=self.dsss.ds.indexlist) #only calculate scores within indexlist
                    elif sm['type']=='scaledsp':
                        pass
        for index in poslist:
            np.dot(self.ratioarray,self.scorearray[:,index[0]:index[1]],self.score[index[0]:index[1]])

    def assess_rrf(self,rrf=[],report='single'):
        raise Exception('Function not working scorer.assess_rrf')
        if self.bmtype['type']=='dsscore':
            self.score[...]=self.basescore[...]
            sp=0
            for i in range(0,len(self.searches)):
                search=self.searches[i]
                if search['object']['type']=='sf':
                    scorer=search['scorer']
                    ep=sp+scorer.idist.shape[-1]
                    self.score+=scorer.get_score(indexlist=self.dsss.ds.indexlist,sfv=rrf[sp:ep])
                    sp=ep
            res=self.dsss.analyze_score(slevel=self.bmtype['criteria'],score=self.score,report=report)
            print res
            return res

    def assess_model(self,report='single',slevel=''):
        #if self.loadfilter!=None:
        #    self.filter(ts.loadfilter)
        #    self.loadfilter=None
        if self.bmtype['type']=='dsscore':
            if len(self.dsss.ds.indexlist)==0:
                if report=='single':
                    return 0
                elif report=='full':
                    return [0,'empty']

            try:
                self.get_model_dsscore()
            except OutOfRange,e:
                return 999.999999
            except NanInScore, e:
                return 998.999999
            #except Exception,e:
            #    raise e
            if not slevel:
                slevel=self.bmtype['criteria']
            res=self.dsss.analyze_score(slevel=slevel,score=self.score,report=report)
            #print "assessing model result: "+str(self.dsss.ds.indexlist)+'  '+str(res)
            self.result=res
            return res
        elif self.bmtype['type']=='sprefine':
            if not self.sprefiner.envinitialized:
                self.sprefiner.initialize_runenv()
            #self.write_potential(type='lib')
            #self.sprefiner.refpot=[self.atomclasslib,self.potlib]
            return -self.sprefiner.assess()
    
    def assess(self,values=[],pos=[],report='single',slevel=''):
        if values==[]:
            print "parameter values are needed for scorer.assess"
            return 0
        if len(pos)==0:
            try:
                self.parvalues[...]=values
            except Exception,e:
                print e
                pdb.set_trace()
        else:
            if not self.parvalues.flags['WRITEABLE']: 
                self.parvalues=np.copy(self.parvalues)
            self.parvalues[pos]=values
        self.assign_values2model()
        # do parallel jobs if specified, otherwise do single job
        #pdb.set_trace()        
        if report=='single' and runenv.numofthread>0:
            for singlescorer in self.singlescorerlist:
                runenv.queue.put([singlescorer.assess_model,[report, slevel]])
            runenv.queue.join()
            scoresum=0
            for singlescorer in self.singlescorerlist:
                scoresum+=singlescorer.result
            if len(self.singlescorerlist)==0:
                return 0
            else:
                return scoresum/len(self.singlescorerlist)        
        else:
            return self.assess_model(report,slevel)

    def assess_block(self,values=[],blockind=0,report='single',slevel=''):
        if values==[]:
            print "parameter values are needed for scorer.assess_block"
            return 0
        try:
            self.parvalues[self.searchblocklist[blockind][-1][0]:self.searchblocklist[blockind][-1][1]]=values
        except:
            pdb.set_trace()
        self.blockind=blockind
        self.assign_values2model()
        if report=='single' and runenv.numofthread>0:
            for singlescorer in self.singlescorerlist:
                runenv.queue.put([singlescorer.assess_model,[report, slevel]])
            runenv.queue.join()
            scoresum=0
            for singlescorer in self.singlescorerlist:
                scoresum+=singlescorer.result
            result=scoresum/len(self.singlescorerlist)        
        else:
            result=self.assess_model(report,slevel)        
        self.blockind=-999
        return result
        
    def assign_values2model(self):
        if self.blockind>=0:
            print "block# "+str(self.blockind)
            il=self.searchblocklist[self.blockind][0]
        else:
            il=range(len(self.searchlists))
        for i in il:
            l=self.searchlists[i]
            si=self.searchlistspos[i]
            ei=self.searchlistspos[i+1]
            l[0:]=self.parvalues[si:ei]

    def get_ds_part(self,sample):
        if self.bmtype['type']=='dsscore':
            #get part of the basescore, the self.sfscore.idst, and the corresponding dsss.
            nop=copy.copy(self)
            nop2=copy.copy(self)
            nop2.dsss=copy.copy(self.dsss)
            nop.dsss=copy.copy(self.dsss)
            [nop.dsss.ds,nop2.dsss.ds]=nop.dsss.ds.get_ds_part(sample)
            nop.dsss.initialize_slevel(nop.dsss.slevel)
            if len(nop2.dsss.ds.indexlist)>0:
                nop2.dsss.initialize_slevel(nop2.dsss.slevel)
            return [nop,nop2]
        elif self.bmtype['type']=='sprefine':
            nop=copy.deepcopy(self)
            nop.sprefiner=nop.sprefiner.get_ds_part(sample)
            return nop

    def get_ds_list(self):
        if self.bmtype['type']=='dsscore':
            #get part of the basescore, the self.sfscore.idst, and the corresponding dsss.
            noplist=[copy.copy(self) for code in self.dsss.ds.codelist]
            nopdslist=self.dsss.ds.get_ds_list()
            for nop,ds in zip(noplist,nopdslist):
                nop.dsss=copy.copy(self.dsss)
                nop.dsss.ds=ds
                nop.dsss.initialize_slevel(nop.dsss.slevel)
            return noplist
        elif self.bmtype['type']=='sprefine':
            raise Exception('Not implemented yet, scorer.get_ds_list')
    
    def assess_basescore(self,report='single',slevel=''):
        if self.bmtype['type']=='dsscore':
            if not slevel:
                slevel=self.bmtype['criteria']
            return self.dsss.analyze_score(slevel=slevel,score=self.basescore,report=report)

    def build_model(self,par):
        self.parvalues[...]=par
        self.assign_values2model()
        #self.assess(par)
        nm=copy.deepcopy(self.model)
        return nm
    
    def write_potential(self,filetype='hdf5',affix='opt',permute=False):
        atompairclustering=False
        decomposition=False
        ratiolist=[]
        svdilist=[]
        pathlist=[]
        refdict=dict([(id(r),r) for r in self.model['scorers'] if r['type']=='sf'])
        for spmodel in self.model['scorers']:
            ratiolist=[]
            if spmodel['type']!='scaledsp':
                continue
            if not atompairclustering and 'bm' in spmodel and re.search('(ap[0-9]+)n',spmodel['bm']):
                rer=re.search('(ap[0-9]+)n',spmodel['bm'])
                apname=rer.group(1)
                atompairclustering=True
                apca=np.load(runenv.libdir+apname+'.npy')
                reflist=[]
            if 'bm' in spmodel and re.search('(svd[0-9]+)',spmodel['bm']):
                decomposition=True
                reflist=[]
                rer=re.search('svd([0-9]',spmodel['bm'])
                svdilist.append(rer.group(1))                
            ref=np.array(0)
            if len(spmodel['refs'])==0:
                pdb.set_trace()
                raise Exception('NO REFS')
            #pdb.set_trace()
            for sfmodel in spmodel['refs']:
                sfv=sf(**sfmodel).get_sf()
                ref=np.squeeze(extendsum(ref,sfv))
                del refdict[id(sfmodel)]
            if atompairclustering:
                reflist.append(ref)
            ratiolist.append(spmodel['ratio'])
            print ref
            ssp=scaledsp(model=spmodel)
            if atompairclustering:
                ssp.add_ref(reflist,affix,filetype,apca=apca,ratio=ratiolist,permute=permute)
            elif decomposition:
                ssp.add_ref(ref,affix,filetype,ratio=ratiolist,svdilist=svdilist,permute=permute)
            else:
                ssp.add_ref(ref,affix,filetype,ratio=ratiolist,permute=permute)
            pathlist.append(ssp.ssppath)
        for i, sfmodel in refdict.items():
            sfo=sf(**sfmodel)
            sfo.write_potential(affix,self.model['bmtype']['bm'],ratio=sfmodel['ratio'])
        print pathlist
        return pathlist

    def get_sample(self, perc=0):
        if self.bmtype['type']=='dsscore':
            if perc==0:
                return copy.deepcopy(self.dsss.ds.codelist)
            else:
                testsample=[]
                trainsample=[]
                for i in range(len(self.dsss.ds.setpos)-1):
                    wcl=copy.deepcopy(self.dsss.ds.codelist[self.dsss.ds.setpos[i]:self.dsss.ds.setpos[i+1]])
                    sl=int(len(wcl)*perc)
                    random.seed(93487923)#347
                    random.shuffle(wcl)
                    testsample=testsample+copy.copy(wcl[0:sl+1])
                    trainsample=trainsample+copy.copy(wcl[sl+1:])
                return [testsample,trainsample]
        elif self.bmtype['type']=='sprefine':
            return self.sprefiner.loopdict.keys()

    def assess_sascore(self,report='single',slevel=''):
        if self.bmtype['type']=='dsscore':
            if not slevel:
                slevel=self.bmtype['criteria']
            res=self.dsss.analyze_score(slevel=slevel,score=self.dsss.ds.sa['score'],report=report)
            print res
            return res
        
    def randomize_score(self):
        poslist=self.dsss.ds.indexlist
        for i in range(0,len(poslist)):
            np.random.shuffle(self.score[poslist[i][0]:poslist[i][1]])

    def assess_randomized_score(self,report='single',slevel=''):
        if self.bmtype['type']=='dsscore':
            if not slevel:
                slevel=self.bmtype['criteria']            
            res=0
            for i in range(0,100):
                self.randomize_score()
                res=res+self.dsss.analyze_score(slevel=slevel,score=self.basescore,report=report)
            print res/100.0
            return res/100.0

    def assess_ideal(self,report='single',slevel=''):
        if self.bmtype['type']=='dsscore':
            if len(self.dsss.ds.indexlist)==0:
                if report=='single':
                    return 0
                elif report=='full':
                    return [0,'empty']
            if not slevel:
                slevel=self.bmtype['criteria']    
            self.dsss.getidealvalue=True
            res=self.dsss.analyze_score(slevel=slevel,score=self.dsss.ds.sa['score'],report=report)
            self.dsss.getidealvalue=False
            print res
            return res
        
    def assess_worst(self,report='single',slevel=''):
        if self.bmtype['type']=='dsscore':
            if len(self.dsss.ds.indexlist)==0:
                if report=='single':
                    return 0
                elif report=='full':
                    return [0,'empty']
            if not slevel:
                slevel=self.bmtype['criteria']    
            self.dsss.getworstvalue=True
            res=self.dsss.analyze_score(slevel=slevel,score=self.dsss.ds.sa['score'],report=report)
            self.dsss.getworstvalue=False
            print res
            return res

    def calc_structurelistscore(self,structurepath,structurelist):
        #Loop through the scorers and calculated score using the routinue for calculating scores from a single alignment file
        os.chdir(structurepath)
        self.apcidist=[]
        self.apcscore=[]
        score=np.zeros(len(structurelist))
        if self.bmtype['combine']=='scoresum':
            print os.getcwd()
            pir(path=os.path.join('pdbs0')).make_pir_fromlist(structurelist)
            for scorer in self.model['scorers']:
                score+=self.calc_structurescore_singlescorer(scorer,structurepath)
        else:
            raise Exception('this type of model is not supported yet.')
        scorestr=[f+', '+str(ss) for f,ss in zip(structurelist,list(score))]
        fh=open('score','w')
        fh.write('\n'.join(scorestr))
        fh.close()
        return score
        
    def calc_structurescore(self,structurepath):
        #Loop through the scorers and calculated score using the routinue for calculating scores from a single alignment file
        os.chdir(structurepath)
        fl=os.listdir(structurepath)
        fl=[f for f in fl if not (f.endswith('hdf5') or f.endswith('pickle') or f=='pdbs0')]
        fl=[f for f in fl if f.endswith('pdb')]
        self.calc_structurelistscore(structurepath,fl)
  
    def calc_structurescore_singlescorer(self,scorer,structurepath):
        #Use the routinue for calculating scores from a single alignment file
        if 'bm' in scorer and re.search('ap[0-9]{1,3}n[0-9]{1,3}',scorer['bm']):
            apc=True
            rer=re.search('ap[0-9]{1,3}n([0-9]{1,3})',scorer['bm'])
            n=int(rer.group(1))
        else:
            apc=False
        if scorer['type']=='scaledsp':
            if apc==True:
                if len(self.apcscore)==0:
                    if not 'bm' in scorer:
                        scorer['bm']=self.bmtype['bm']
                    scorer['dslist']=self.bmtype['dslist']
                    spdss=spdsscore(model=scorer)
                    os.chdir(structurepath)
                    spdss.calc_score(0,reportstatus=False)
                    self.apcscore=np.load('0.score.npy')
                finalscore=self.apcscore[:,n]
            else:
                if not 'bm' in scorer:
                    scorer['bm']=self.bmtype['bm']
                scorer['dslist']=self.bmtype['dslist']
                spdss=spdsscore(model=scorer)
                os.chdir(structurepath)
                spdss.calc_score(0,reportstatus=False)                
                finalscore=spdss.read_singlescorefile('0.score')
        elif scorer['type']=='sf':
            sc=copy.deepcopy(scorer)
            del sc['type']
            #pdb.set_trace()
            if not 'bm' in sc:
                sc['bm']=self.bmtype['bm']
            sfo=sf(**sc)
            if apc==True:
                if len(self.apcidist)==0:
                    rsp=rawsp(pdbset='_'.join(self.bmtype['dslist']),features=sc['features'],genmethod=sc['bm'],decoy=True,routine='calc_sp_i')
                    os.chdir(structurepath)
                    idist=rsp.calc_sp_i(0,reportstatus=False)
                    self.apcidist=idist
                idist=self.apcidist[:,n]
            else:
                rsp=rawsp(pdbset='_'.join(self.bmtype['dslist']),features=sc['features'],genmethod=sc['bm'],decoy=True,routine='calc_sp_i')
                os.chdir(structurepath)
                idist=rsp.calc_sp_i(0,reportstatus=False)
            sfv=np.squeeze(sfo.get_sf())
            finalscore=np.dot(idist,sfv)
        if 'ratio' in scorer:
            finalscore=finalscore*scorer['ratio']
        return finalscore

    def plot_rfs(self,r=5,c=6):
        k=0
        for sm,st,so in zip(self.model['scorers'],self.scorertype,self.scorerlist):
            if st==0:
                continue
            elif st==1:
                pass
            elif st==2:
                plt.subplot(r,c,k)
                k=k+1
                if sm['type']=='sf':
                    so.parvalue[:]=sm['parvalue']
                    plt.plot(np.exp(so.get_sf()))
                elif search['object']['type']=='scaledsp':
                    pass
              