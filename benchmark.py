"""
   SOAP benchmark module, calculating statistics on the scores.

"""
from env import *
from decoys import DecoySets
from rankScores import dsscore

class dsscorestats(object):
    """
    Summarize the ranking results based on RMSDs and scores.

    :param list dslist: the names of decoy sets used for benchmark, for retrieving the RMSD values
    :param dsscore dsscore: the score object :class:`rankScores.dsscore`, with the scores

    """
    def __init__(self,dslist=[],dsscore=dsscore()):
        self.dsscore=dsscore
        self.score=dsscore.score
        self.dslist=dslist
        self.ds=[]
        self.get_rmsdlist()
        self.slevel=''
        self.setcombine='raw'
        self.getidealvalue=False
        self.getworstvalue=False
        #self.sa_fields=[]
        #self.sa={}# reorder sa to save it in continuous array
        #self.sal=[]

    def get_rmsdlist(self):
        self.ds=DecoySets(self.dslist)
        #clean up the DecoySets object to eliminate unnecessary objects
        self.ds=self.ds.load(mode='dss')
        if not self.ds.sa.flags['OWNDATA']:
            self.ds.sa=copy.deepcopy(self.ds.sa)

    def analyze_score(self,slevel='NativeSelection',dsscoreo=[],score=[], report='single'):
        """
        Calculate the performance measure.

        :param str slevel: benchmark criteria/performance measure
        :param dsscore dsscoreo: the score object :class:`rankScores.dsscore`, with the scores
        :param list score: the scores
        :param str report: single|full, report single number summarizing the performance or detailed measures

        """
        #reported result are the large the better
        if len(score)>0:
            self.score=score
        elif dsscoreo:
            self.dsscore=dsscoreo
            self.score=dsscoreo.score
        if slevel!=self.slevel:
            self.initialize_slevel(slevel)
        if self.slevel=='dgcc' or slevel=='dgcc':
            return np.corrcoef(score,self.dga)[0,1]
        self.statsarray=np.zeros([len(self.ds.indexlist),self.numofc])
        for key in self.criteriadict:
            if key =='NativeSelection':
                self.analyze_score_0()
            elif key=='top':
                self.analyze_topmodels()
            elif key=='cc':
                self.analyze_cc()
            elif key=='bestrank':
                self.analyze_bestrank()
            elif key=='lbslope':
                self.analyze_lbslope()
            elif key=='dcg':
                self.analyze_dcg()
            elif key=='enrichment':
                self.analyze_enrichment()
            else:
                raise Exception('Do not know how to analyze score with criteria: '+slevel)
        return self.report_stats(report)

    def report_stats(self,report):
        csa=np.copy(self.statsarray)
        #print csa
        #print csa.mean()
        rv=0
        #report a single number indicating the performance of the potential
        for i in range(self.numofc):
            if self.setcombine=='raw':
                rv+=(self.criterialist[i]['ratio']*csa[:,i]).mean()
            elif self.setcombine=='equal':
                srv=0
                for j in range(len(self.ds.setpos)-1):
                    srv+=csa[self.ds.setpos[j]:self.ds.setpos[j+1]].mean()
                srv=srv/(len(self.ds.setpos)-1)
                rv+=self.criterialist[i]['ratio']*srv
        if report=='single':
            return rv
        ra=[]
        #report details on each criteria
        for i in range(self.numofc):
            ra.append((self.statsarray[:,i].mean(),self.statsarray[:,i].std()/np.sqrt(self.statsarray.shape[0])))
        if report=='detail': #different criteria
            return [rv, ra]
        fra={}
        nod=[]
        #report full details on each criteria on each set
        for i in range(len(self.ds.setpos)-1):
            nod.append(self.ds.setpos[i+1]-self.ds.setpos[i])
            sa=self.statsarray[self.ds.setpos[i]:self.ds.setpos[i+1]]
            sra=[]
            for j in range(self.numofc):
                sra.append((sa[:,j].mean(),sa[:,j].std()/np.sqrt(sa.shape[0]),sa[:,j].min(),sa[:,j].max()))
            fra[self.ds.dslist[i]]=sra
        if report=='full': #different criteria, different set,number of data(code)
            return [rv, ra,fra,nod]
        else:
            raise Exception('Type '+report+' is not defined in dsscorestats.report_stats()')

    def initialize_slevel(self,slevel):
        if slevel=='dgcc': #whole sets correlation?
            poslist=self.ds.indexlist
            na=[]
            for pos in poslist:
                na=na+list(self.ds.sa['dg'][pos[0]:pos[1]])
            self.dga=np.array(na)
            return 0
        if self.ds.sa[0]['rmsd']==0:
            self.withnative=True
        else:
            self.withnative=False
        sll=slevel.split(':')
        if len(sll)==2:
            self.setcombine=sll[0]
            slevel=sll[1]
        sll=slevel.split('+')
        criteriadict={}
        ctsn=0
        self.numofc=len(sll)
        self.criterialist=[] # dictionary and list stores the same information...
        for sl in sll:
            rer=re.search("([0-9\.]{1,10}x)",sl)
            if rer:
                ratio=float(rer.group(1)[:-1])
                sl=sl[len(rer.group(1)):]
            else:
                ratio=1
            if sl.startswith('NativeSelection'):
                criteriadict['NativeSelection']={'ctsn':ctsn,'ct':'NativeSelection','ratio':ratio}
                self.criterialist.append(criteriadict['NativeSelection'])
            elif sl.startswith('cc'):
                criteriadict['cc']={'ctsn':ctsn,'ct':'cc','ratio':ratio}
                self.criterialist.append(criteriadict['cc'])
            elif sl.startswith('dcg'):
                criteriadict['dcg']={'ctsn':ctsn,'ct':'dcg','ratio':ratio}
                self.criterialist.append(criteriadict['dcg'])
            elif sl.startswith('enrichment'):
                criteriadict['enrichment']={'ctsn':ctsn,'ct':'enrichment','ratio':ratio}
                self.criterialist.append(criteriadict['enrichment'])
                self.ds.sa['rmsd']=1-self.ds.sa['rmsd']#make 1 be postiive example, 0 negative example
            elif sl.startswith('lbslope'):
                criteriadict['lbslope']={'ctsn':ctsn,'ct':'lbslope','ratio':ratio}
                self.criterialist.append(criteriadict['lbslope'])
                self.initialize_lbslope(sl[7:])
            elif sl.startswith('top'):
                if not 'top' in criteriadict:
                    criteriadict['top']=[]
                criteriadict['top'].append(self.initialize_topmodels(sl))
                criteriadict['top'][-1]['ctsn']=ctsn
                criteriadict['top'][-1]['ratio']=ratio
                self.criterialist.append(criteriadict['top'][-1])
            elif sl=='bestrank':
                criteriadict['bestrank']={'ctsn':ctsn,'ct':'bestrank','ratio':ratio}
                self.criterialist.append(criteriadict['bestrank'])
            else:
                raise Exception('Type '+sl+' is not defined in dsscorestats.initialize_slevel()')
            ctsn=ctsn+1
        self.criteriadict=criteriadict
        if 'top' in self.criteriadict:
            self.initialize_topmodels_stats()
            #self.initialize_topmodels_sa()
        self.slevel=slevel
        if not self.withnative and slevel=='NativeSelection':
            raise Exception('The decoy set does not have native structure inside. You can not use NativeSelection as a criteria')


    def initialize_lbslope(self,lbpars):
        #initialize the rmsd bins, you should make sure the decoy set contain nativs
        if not self.withnative:
            raise Exception("No native structures, can not calculate lbslope")
        lbpl=lbpars.split('-')
        if len(lbpl)<2:
            raise Bugs("we need both the rmsd range and the number of bins to calculate the lowerbond")
        elif len(lbpl)==3:
            withnative=True
        else:
            withnative=False
        poslist=self.ds.indexlist
        self.rmsdbinlist=[]
        for i in range(0,len(poslist)):
            rmsds=self.ds.sa['rmsd'][poslist[i][0]:poslist[i][1]]
            rmsdmax=min(rmsds[-1],float(lbpl[0]))
            rmsdmin=0
            if withnative:
                raindex=[0,1]
                ramean=[0]
            else:
                raindex=[0]
                ramean=[]
            prmsdcutoff=rmsdmin
            ramean
            for rmsdcutoff in np.linspace(rmsdmin,rmsdmax,int(lbpl[1]))[1:]:
                cind=np.searchsorted(rmsds,rmsdcutoff)
                if cind>raindex[-1]+5:# each bin contain at leat
                    raindex.append(cind)
                    ramean.append(0.5*(rmsdcutoff+prmsdcutoff))# the center
                    prmsdcutoff=rmsdcutoff
            self.rmsdbinlist.append([raindex,np.array(ramean)])

    def initialize_topmodels(self,slevel):
        """
        slevel: define the stats we want to calculateate
            top+(cn)+'_'+(cp)+_+(cv)+_+(cf)
        cn: number of top models to look at
        cf: top model filter: defines which part of the top model we are looking at
            lessthan: bool values whether the rmsd is less than the specified value,"rmsd10"
            First:take the property of the first model
            None: or no filter
        cp: propertis of the filtered top model set
            (rmsd,rank,rmsddiff,irmsd,rrank,rrankr,rlrank,rlrankr,"")
        cs: combine the propertite, default
            (mean,min,sum,max)
            "": the value itself, only for first pass models.
            len:len()
            nonempty: len>0
            perc: percentage in total
        Examples::
            top1000_nonempty__rmsd10ORirmsd4 # whether there is a model sastifi
            top1000_len__rmsd10ORirmsd4 # the number of such models
            top1000_sum_revrank_rmsd10ORirmsd4FIRST
        cv: the values to look at for those models.
                The values of certain rmsd
                The value diff of certain rmsd
                The rank
                    Certain values combine the rank and the rmsd...
                The filter values
            rmsd
        """
        ctd={}# criteria dictionary 'tn':totoal number, 'cd':rmsd dictionary,'ct':criteria type
        sl=slevel.split('_')
        if len(sl)!=4:
            raise Exception('Do not know how to handle this benchmark criteria '+slevel)
        ctd['cn']=float(sl[0][3:])
        #if sl[1].endswith('rmsd'):
        #    self.sa_fields.append(sl[1])
        #elif sl[1].endswith('diff'):
        #    self.sa_fields.append(sl[1])
        ctd['cp']=sl[1]
        ctd['cs']=sl[2]
        ctd['cf']=sl[3]
        ctd['cfd']={}
        if sl[3]:
            filters=sl[3]
            if filters.endswith('FIRST'):
                ctd['cfd']['firstonly']=True
                filters=filters[:-5]
            else:
                ctd['cfd']['firstonly']=False
            fl=filters.split('OR')
            ctd['cfd']['cl']=[]
            for item in fl:
                if re.search('([a-z]+)([0-9\.]+)',item):#whether the structure complies to the criteria
                    rerl=re.search('([a-z]+)([0-9\.]+)',item)
                    #self.sa_fields.append(rerl.group(1))
                    ctd['cfd']['cl'].append([rerl.group(1),float(rerl.group(2))])
                else:# the rmsd ... value
                    raise Exception('Do not know how to handle this benchmark criteria '+slevel)
            #initialize total number for percentage calculation
            ctd['cfn']=[]
            poslist=self.ds.indexlist
            for i in range(0,len(poslist)):
                if self.withnative:
                    csa=self.ds.sa[poslist[i][0]+1:poslist[i][1]]
                else:
                    csa=self.ds.sa[poslist[i][0]:poslist[i][1]]
                bool=np.zeros(len(csa),dtype=np.bool)
                for rmsd, cutoff in ctd['cfd']['cl']:
                    bool.__ior__(csa[rmsd]<cutoff)
                #bni=np.nonzero(bool)[0] # index for the sortind
                ctd['cfn'].append(float(bool.sum())+0.00000000000001)
        return ctd

    def initialize_topmodels_stats(self):
        self.top_numarray=np.zeros([len(self.ds.indexlist),len(self.criteriadict['top'])+1],dtype=np.int)
        for i in range(0,len(self.ds.indexlist)):
            for sd,j in zip(self.criteriadict['top'],range(0,len(self.criteriadict['top']))):
                if sd['cn']<1:
                    self.top_numarray[i,j]=int(sd['cn']*(self.ds.indexlist[i][1]-self.ds.indexlist[i][0]))
                else:
                    self.top_numarray[i,j]=min(int(sd['cn']),self.ds.indexlist[i][1]-self.ds.indexlist[i][0])
            self.top_numarray[i,-1]=np.max(self.top_numarray[i,:-1])
        top_nummax=self.top_numarray.max()
        #self.top_valuearray=np.zeros(top_nummax,dtype=self.ds.sa.dtype) # the array for store changed sa value
        #self.top_propertyarray=np.zeros(top_nummax) #array for store intermidiate results.
        self.top_boolarray=np.zeros(top_nummax,dtype=np.bool) # the array for store bool values

    def initialize_topmodels_sa(self):
        for key in self.sa_fields:
            self.sa[key]=np.copy(self.ds.sa[key])
        poslist=self.ds.indexlist
        #loop through different decoy set
        for i in range(0,len(poslist)):
            self.sal.append({})
            for key in self.sa_fields:
                if self.withnative:
                    self.sal[-1][key]=self.sa[key][poslist[i][0]+1:poslist[i][1]]
                else:
                    self.sal[-1][key]=self.sa[key][poslist[i][0]:poslist[i][1]]

    def analyze_topmodels(self):
        #calc the success rate of top models to be in some rmsd range.slevel top10_rmsd2.0_irmsd2.0
        rmsdnpa=self.score
        poslist=self.ds.indexlist
        topcd=self.criteriadict['top']
        #loop through different decoy set
        for i in range(0,len(poslist)):
            #pdb.set_trace()
            if self.withnative:
                score=rmsdnpa[poslist[i][0]+1:poslist[i][1]]
                rmsds=self.ds.sa[poslist[i][0]+1:poslist[i][1]]
            else:
                score=rmsdnpa[poslist[i][0]:poslist[i][1]]
                rmsds=self.ds.sa[poslist[i][0]:poslist[i][1]]
            tn=self.top_numarray[i]
            statsarray=self.statsarray[i]
            #analyze single sequence
            ctn=tn[-1]
            if self.getidealvalue:
                rmsdcriteria='rmsd'
                for tpc,k in zip(topcd,range(len(topcd))):
                    if tpc['cp'].endswith('rmsd') or tpc['cp'].startswith('rmsd'):
                        rmsdcriteria=tpc['cp']
                lowestctn=bottleneck.argpartsort(rmsds[rmsdcriteria],ctn)
                lsind=np.argsort(rmsds[rmsdcriteria][lowestctn[:ctn]])
                sortind=lowestctn[:ctn][lsind]
            elif self.getworstvalue:
                rmsdcriteria='rmsd'
                for tpc,k in zip(topcd,range(len(topcd))):
                    if tpc['cp'].endswith('rmsd') or tpc['cp'].startswith('rmsd'):
                        rmsdcriteria=tpc['cp']
                lowestctn=bottleneck.argpartsort(-rmsds[rmsdcriteria],ctn)
                lsind=np.argsort(-rmsds[rmsdcriteria][lowestctn[:ctn]])
                sortind=lowestctn[:ctn][lsind]
            else:#testing code
                lowestctn=bottleneck.argpartsort(score,ctn)
                lsind=np.argsort(score[lowestctn[:ctn]])
                sortind=lowestctn[:ctn][lsind]
                #print rmsds[sortind]
            for k,tpc in enumerate(topcd):
                #set up filter:
                if tpc['cf']:
                    bool=self.top_boolarray[0:tn[k]]
                    bool[:]=False
                    for rmsd, cutoff in tpc['cfd']['cl']:
                        try:
                            bool.__ior__(rmsds[rmsd][sortind[0:tn[k]]]<cutoff)
                        except:
                            pdb.set_trace()
                    bni=np.nonzero(bool)[0] # index for the sortind
                    if tpc['cfd']['firstonly'] and len(bni)>0:
                        #if self.ds.codelist[i].endswith('1I2M'):
                        #    pdb.set_trace()
                    #    print self.ds.codelist[i]+' '+str(bni[0])+' '+str(rmsds['rmsd'][sortind[bni[0]]])+' '+str(rmsds['irmsd'][sortind[bni[0]]])
                        bni=bni[0:1]
                    #else:
                    #    print self.ds.codelist[i]+' '+str(-1)+' '+str(0)+' '+str(0)
                else:
                    bni=range(0,tn[k])# index list for the sortind
                #calculate property
                if tpc['cp']=='':
                    pa=bni
                elif tpc['cp'].endswith('rmsd') or tpc['cp'].startswith('rmsd'):
                    pa=[-rmsds[i][tpc['cp']] for i in sortind[bni]] #negative to make the large the better
                elif tpc['cp'].endswith('diff'):
                    pa=-(rmsds[tpc['cp'][:-4]][sortind[bni]].mean()-rmsds[tpc['cp'][:-4]][0:len(bni)]) #negative to make the large the better
                elif tpc['cp']=='rrank':
                    if len(bni)==0:
                        pa=np.array([0])
                    else:
                        pa=tn[k]+1-bni
                elif tpc['cp']=='rlrank':
                    if len(bni)==0:
                        pa=np.array([0])
                    else:
                        pa=np.log10(tn[k]+2)-np.log10(bni+1)
                elif tpc['cp']=='rrankr': #with reverse rmsd
                    if len(bni)==0:
                        pa=np.array([0])
                    else:
                        pa=tn[k]+1-bni+cutoff-rmsds[rmsd][sortind[bni]]
                elif tpc['cp']=='rlrankr':
                    if len(bni)==0:
                        pa=np.array([0])
                    else:
                        pa=np.log10(tn[k]+2+cutoff)-np.log10(bni+rmsds[rmsd][sortind[bni]]+1)
                else:
                    pdb.set_trace()
                    raise Exception('can not handle benchmark criteria cp '+str(tpc['cp']))
                #calculate return value by combining the property
                if tpc['cs']=='':
                    try:
                        statsarray[tpc['ctsn']]=pa
                    except:
                        pdb.set_trace()
                elif tpc['cs']=='mean':
                    statsarray[tpc['ctsn']]=sum(pa)/len(pa)
                elif tpc['cs']=='median':
                    statsarray[tpc['ctsn']]=np.median(pa)
                elif tpc['cs']=='sum':
                    statsarray[tpc['ctsn']]=np.sum(pa)
                elif tpc['cs']=='min':
                    statsarray[tpc['ctsn']]=np.min(pa)
                elif tpc['cs']=='len':
                    statsarray[tpc['ctsn']]=len(pa)
                elif tpc['cs']=='perc':
                    statsarray[tpc['ctsn']]=len(pa)/tpc['cfn'][i]
                elif tpc['cs']=='nonempty':
                    statsarray[tpc['ctsn']]=len(pa)>0
                else:
                    raise Exception('can not handle benchmark criteria cs '+str(tpc['cs']))
        #print self.statsarray

    def analyze_score_0(self):
    #calculate benchmark statistics based on scores
        rmsdnpa=self.score
        poslist=self.ds.indexlist
        sn=0
        ai=self.criteriadict['NativeSelection']['ctsn']
        for i in range(0,len(poslist)):
            score=rmsdnpa[poslist[i][0]:poslist[i][1]]
            if (np.argmin(score)==0 or self.getidealvalue) and (not self.getworstvalue):
                sn=sn+1
                self.statsarray[i,ai]=1

    def analyze_cc(self):
    #calculate benchmark statistics based on scores
        rmsdnpa=self.score
        poslist=self.ds.indexlist
        ai=self.criteriadict['cc']['ctsn']
        for i in range(0,len(poslist)):
            score=rmsdnpa[poslist[i][0]:poslist[i][1]]
            rmsds=self.ds.sa['rmsd'][poslist[i][0]:poslist[i][1]]
            if self.getidealvalue:
                self.statsarray[i,ai]=1
            elif self.getworstvalue:
                self.statsarray[i,ai]=-1
            else:
                self.statsarray[i,ai]=np.corrcoef(score,rmsds)[0,1]
                if np.isnan(self.statsarray[i,ai]):
                    self.statsarray[i,ai]=-1

    def analyze_dcg(self):#discounted cumulative gain
        rmsdnpa=self.score
        poslist=self.ds.indexlist
        ai=self.criteriadict['dcg']['ctsn']
        for i in range(0,len(poslist)):
            score=rmsdnpa[poslist[i][0]:poslist[i][1]]
            rmsds=self.ds.sa['rmsd'][poslist[i][0]:poslist[i][1]]
            if self.getidealvalue:
                sind=np.arange(poslist[i][1]-poslist[i][0])
            elif self.getworstvalue:
                sind=np.arange(poslist[i][1]-poslist[i][0]-1,-1,-1)
            else:
                sind=np.argsort(score)
            ta=(500-rmsds)/(sind+1)
            self.statsarray[i,ai]=(ta[ta>0]).sum()

    def analyze_enrichment(self):
        """
            Calculate the enrichment score
        """
        rmsdnpa=self.score
        poslist=self.ds.indexlist
        ai=self.criteriadict['enrichment']['ctsn']
        for i in range(0,len(poslist)):
            score=rmsdnpa[poslist[i][0]:poslist[i][1]]
            rmsds=self.ds.sa['rmsd'][poslist[i][0]:poslist[i][1]]
            if self.getidealvalue:
                sind=np.arange(poslist[i][1]-poslist[i][0])
            elif self.getworstvalue:
                sind=np.arange(poslist[i][1]-poslist[i][0]-1,-1,-1)
            else:
                sind=np.argsort(score)
            va=-np.log(np.nonzero(rmsds[sind])[0])
            self.statsarray[i,ai]=va.sum()/len(va)

    def analyze_lbslope(self):
    #calculate the slope of the lower bound line in the score vs rmsd file.
        # the lower bound rmsd bins are precalculated as saved in self.rmsdbins
        rmsdnpa=self.score
        poslist=self.ds.indexlist
        ai=self.criteriadict['lbslope']['ctsn']
        for i in range(0,len(poslist)):
            if len(self.rmsdbinlist[i][0])<4:
                continue
            score=rmsdnpa[poslist[i][0]:poslist[i][1]]
            if self.getidealvalue:
                self.statsarray[i,ai]=1000
            elif self.getworstvalue:
                self.statsarray[i,ai]=-1000
            else:
                y=[]
                for j in range(len(self.rmsdbinlist[i][0])-1):
                    y.append(score[self.rmsdbinlist[i][0][j]:self.rmsdbinlist[i][0][j+1]].min())
                self.statsarray[i,ai]=np.polyfit(self.rmsdbinlist[i][1],y,1)[0]

    def analyze_bestrank(self):
    #calculate benchmark statistics based on scores
        rmsdnpa=self.score
        poslist=self.ds.indexlist
        ai=self.criteriadict['bestrank']['ctsn']
        for i in range(0,len(poslist)):
            if self.withnative:
                score=rmsdnpa[poslist[i][0]+1:poslist[i][1]]
                rmsds=self.ds.sa['rmsd'][poslist[i][0]+1:poslist[i][1]]
                si=1
            else:
                score=rmsdnpa[poslist[i][0]:poslist[i][1]]
                rmsds=self.ds.sa['rmsd'][poslist[i][0]:poslist[i][1]]
                si=0
            if self.getidealvalue:
                self.statsarray[i,ai]=1
            elif self.getworstvalue:
                self.statsarray[i,ai]=len(score)
            else:
                self.statsarray[i,ai]=(score<score[si]).sum()+1

    def save(self):
        fh=open(self.dsscore.scorepath+'.stats.pickle','wb')
        cPickle.dump(self,fh)
        fh.close()

    def load(self):
        fh=open(self.score.scorepath+'.stats.pickle','rb')
        no=cPickle.load(fh)
        fh.close()
        return no

    def log(self):
        sdict={}
        result={}
        result['detailedstats']=self.sa
        result['setstats']=self.dssa
        result['overallstats']=self.asa
        sdict['catagory']='dsscore'
        sdict['result']=result
        sdict['method']=self.dsscore.method
        sdict['path']=self.dsscore.scorepath+'.stats.pickle'
        sdict['string']=self.dsscore.scorename+self.statsstring
        self.sdict=sdict
        runenv.log(sdict)
