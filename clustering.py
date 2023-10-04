"""
   SOAP clustering module, for clustering recovery functions.

"""
from __future__ import print_function
import numpy as np
import scipy.cluster

class singlecluster(object):
    """
    """
    def __init__(self, parindex=[],allclusters=None):
        self.pars=allclusters.pars[parindex,:]
        self.results=allclusters.results[parindex]
        self.bestresult=allclusters.results[parindex].max()
        self.worstresult=allclusters.results[parindex].min()
        self.parindex=parindex
        self.allclusters=allclusters
        self.min=self.pars.min(axis=0)
        self.max=self.pars.max(axis=0)
        self.middle=(self.min+self.max)/2
        self.mean=self.pars.mean(axis=0)
        self.range=self.max-self.min
        self.analyze()

    def find_overlap(self,region):
        max2=np.minimum(self.max,region[1])
        min2=np.maximum(self.min,region[0])
        if (max2>min2).all():
            return [min2,max2]
        else:
            return []

    def get_pars_in_range(self,min2,max2):
        print("not working")
        pdb.set_trace()
        ai=[]
        for i in range(0,len(self.pars)):
            if (np.logical_and(self.pars[i]>min2, self.pars[i]<max2)).all():
                ai.append(i)
        for i in range(0,len(self.otherpars)):
            if (np.logical_and(self.allclusters.pars[i]>min2, self.allclusters.pars[i]<max2)).all():
                ai.append(i)
        return ai

    def analyze(self):
        self.get_otherpars_within_cluster_boundary()
        self.calc_dist_to_middle()
        self.calc_avedist_vs_aveperformace()
        self.cluster_rankscore()

    def get_otherpars_within_cluster_boundary(self):
        ai=[]
        for i in range(0,len(self.allclusters.pars)):
            if i not in self.parindex and (np.logical_and(self.allclusters.pars[i]>self.min, self.allclusters.pars[i]<self.max)).all():
                ai.append(i)
        self.otherparindex=ai
        self.otherresults=self.allclusters.results[ai]
        if len(self.otherresults)>0:
            self.worstresult=self.otherresults.min()
        self.allpars=np.vstack((self.pars,self.allclusters.pars[self.otherparindex]))
        self.allresults=np.hstack((self.results,self.otherresults))

    def calc_dist_to_middle(self):
        self.distvsresult=np.zeros(len(self.allclusters.results),dtype=[('distance','f4'),('resultdiff','f4')])
        self.distvsresult['distance']=((np.abs(self.allclusters.pars-self.middle)).sum(axis=1))/(np.abs(self.middle).sum())
        self.distvsresult['resultdiff']=(self.allclusters.bestresult-self.allclusters.results)/(self.allclusters.bestresult-self.allclusters.worstresult)
        self.distvsresult.sort(order='distance')

    def get_closest_to(self,value):
        distance=((np.abs(self.pars-value)).sum(axis=1))
        parindex=np.argsort(distance)[0]
        return self.pars[parindex]

    def get_bestpar_close_to(self,value):
        print("accessing middle of the cluster")
        if self.allclusters.optscorer.assess_rrf(np.log(self.middle))==self.bestresult:
            return value
        else:
            return self.get_closest_to(value)

    def calc_avedist_vs_aveperformace(self):
        drds=self.distvsresult
        maxdist=drds['distance'].max()
        distl=[]
        resultl=[]
        startindex=0
        for i in range(len(drds)):
            if drds['resultdiff'][i]!=0:
                if i==0:
                    distl.append(drds['distance'][i])
                else:
                    distl.append(drds['distance'][i-1])
                resultl.append(0)
                startindex=i
                break
        sp=round(np.log10(drds['distance'][startindex]),2)
        sp=max(-5,sp)
        for i in np.arange(sp,1,0.1):
            rs=10**(i-0.05)
            re=10**(i+0.05)
            si,ei=np.searchsorted(drds['distance'],[rs,re])
            if si==ei:
                continue
            ard=drds['resultdiff'][si:ei]
            distl.append(10**i)
            resultl.append((ard.min()+ard.max())/2)
        self.slopearray=np.zeros(len(distl),dtype=[('distance','f4'),('resultdiff','f4')])
        self.slopearray['distance']=np.array(distl)
        self.slopearray['resultdiff']=np.array(resultl)

    def plot_dist_vs_result(self):
        phs=plt.semilogx(self.distvsresult['distance'],self.distvsresult['resultdiff'],'g.')
        phs2=plt.semilogx(self.slopearray['distance'],self.slopearray['resultdiff'])
        plt.xlabel('Percentage distance to the center of the cluster')
        plt.ylabel('Percentage performace difference to the best')
        plt.legend((phs[0],phs2[0]),('Different RRFs','Mean difference'),loc=0)

    def cluster_rankscore(self,diffratio=0.1):
        #determine how good the cluster is, so we can pick the best clusters
        sr=self.allclusters.scoreratio
        score=sr[0]*self.bestresult
        score+=sr[1]*self.worstresult/self.bestresult
        score+=sr[2]*1#self.allclusters.optscorer.assess_rrf(np.log(self.middle))/self.bestresult
        score+=sr[3]*((np.abs(self.max-self.min)).sum()/(np.abs(self.middle).sum()))
        ia=np.nonzero(self.slopearray['resultdiff']>diffratio)[0]
        if len(ia)==0:
            di=-1
        else:
            di=ia[0]
        score+=sr[4]*self.slopearray['distance'][di]
        self.rankscore=score

class optclustering(object):
    """
    """
    def __init__(self,optpars=[],optresults=[],scorer=[],clustermethod={}):
        self.pars=optpars
        self.results=optresults
        self.bestresult=self.results.max()
        self.worstresult=self.results.min()
        self.clustercutoff=2.5
        self.clusterdict={}
        self.bestpars=[]
        self.bestparsindex=[]
        self.clusterindexes=[]
        self.clusters=[]
        self.clusterdict={}
        self.clusteringperc=0.99999
        self.scoreratio=[1,1,0,3,1]
        self.pickwhichone='middle'
        self.singlebestpar=False
        self.hasscorer=False
        if len(scorer)>0:
            self.optscorer=scorer[0]
            self.testscorer=scorer[1]
            self.hasscorer=True
        if len(clustermethod)>0:
            for key in clustermethod:
                self.__dict__[key]=clustermethod[key]

    def get_pars_region(self):
        self.parregion=[0,0]
        self.parregion[0]=self.par.min(axis=0)
        self.parregion[1]=self.par.max(axis=1)
        return self.parregion

    def find_overlap(self,regions):
        alloverlap=[]
        for cl1 in self.clusters:
            for cl2 in regions:
                soverlap=cl1.find_overlap(cl2)
                if soverlap:
                    alloverlap.append(soverlap)
        return alloverlap

    def analyze(self):
        if len(set(self.results))==1:
            print("!!!!!!!!!!Search returns the same results for all the values tried!!!!!!!!!!!")
            self.pickedpar=np.log(self.pars.mean(axis=0))
            self.bestmodelresult=self.optscorer.assess_rrf(self.pickedpar,report='full')
            self.testresult=self.testscorer.assess_rrf(self.pickedpar)
            return self.pickedpar
        percresult=(self.results-self.worstresult)/(self.bestresult-self.worstresult)
        self.bestparsindex=np.nonzero(percresult>self.clusteringperc)[0]
        self.bestpars=self.pars[self.bestparsindex]
        self.clustering_bestpars()
        self.pick_best_cluster(self.clusters)
        self.get_pickedpar()
        print('Best result '+str(self.bestresult))
        self.optscorer.assess_rrf(self.pickedpar)
        if self.hasscorer:
            self.testresult=self.testscorer.assess_rrf(self.pickedpar)
            self.bestmodelresult=self.optscorer.assess_rrf(self.pickedpar,report='full')

    def reduce_size(self,lastone=False):
        if lastone:
            bestnum=1000
            othernum=4000
        else:
            bestnum=200
            othernum=10
            self.clusters=[]
        na=np.hstack([self.pars,self.results.reshape(len(self.results),1)])
        nna=get_representative_pars_forbest(na,maxcluster=bestnum,cutoffdistance=0.0001,clusteringperc=0.999, othermaxnumber=othernum)
        self.pars=nna[:,:-1]
        self.results=nna[:,-1]
        percresult=(self.results-self.worstresult)/(self.bestresult-self.worstresult)
        self.bestparsindex=np.nonzero(percresult>self.clusteringperc)[0]
        self.bestpars=self.pars[self.bestparsindex]

    def get_pickedpar(self):
        cl=self.clusters[self.pickedclusterindex]
        if self.pickwhichone=='middle':
            self.pickedpar=np.log(cl.get_closest_to(cl.middle))
        elif self.pickwhichone=='mean':
            self.pickedpar=np.log(cl.get_closest_to(cl.pars.mean(axis=0)))
        elif self.pickwhichone=='median':
            self.pickedpar=np.log(cl.get_closest_to(np.median(cl.pars,axis=0)))
        elif self.pickwhichone=='last':
            self.pickedpar=np.log(cl.pars[-1])
        elif self.pickwhichone=='rawmiddle':
            self.pickedpar=np.log(cl.get_bestpar_close_to(cl.middle))
        elif self.pickwhichone=='rawmean':
            self.pickedpar=np.log(cl.get_bestpar_close_to(cl.pars.mean(axis=0)))
        elif self.pickwhichone=='rawmedian':
            self.pickedpar=np.log(cl.get_bestpar_close_to(np.median(cl.pars,axis=0)))
        elif self.pickwhichone=='random':
            self.pickedpar=np.log(cl.pars[np.random.randint(0,cl.pars.shape[0],1)[0]])

    def pick_best_cluster(self,clusters):
        if len(self.clusters)==1:
            self.pickedclusterindex=0
            return 0
        rankscore=np.zeros(len(self.clusters))
        for i in range(len(clusters)):
            rankscore[i]=(self.clusters[i].rankscore)
        self.pickedclusterindex=np.argsort(rankscore)[-1]

    def plot(self):
        scale=1
        dc=0.2
        cct=0.1
        ph1=plt.figure(figsize=[12*scale,8*scale],dpi=90)
        self.plot_all_pars(diffcutoff=dc,cmapcutoff=cct,log=False)
        if len(self.clusters)==0:
            return [ph1,0]
        self.plot_clusters_boundary(self.clusters,log=False)
        ph2=plt.figure(figsize=[12*scale,8*scale],dpi=90)
        self.clusters[self.pickedclusterindex].plot_dist_vs_result()
        return [ph1,ph2]

    def plot_clusters_boundary(self,clusters,log=False):
        linestyle=['w--','w-.','w:','w.','w*','w+']
        lph=[]
        lt=[]
        k=0
        for cluster in clusters:
            minline=cluster.min
            maxline=cluster.max
            if log:
                minline=np.log(minline)
                maxline=np.log(maxline)
            plt.axes(self.currentaxis)
            ph=plt.plot(minline,linestyle[k],alpha=0.9,linewidth=3)
            lph.append(ph)
            if k==self.pickedclusterindex:
                lt.append('Cluster #'+str(k+1)+' (best)')
            else:
                lt.append('Cluster #'+str(k+1))
            plt.plot(maxline,linestyle[k],alpha=0.9,linewidth=3)
            k=k+1
        ph=plt.plot(np.exp(self.pickedpar),'-+g',alpha=0.8)
        lph.append(ph)
        lt.append('Picked')
        leg=plt.legend(lph,lt,loc=0,shadow=True)
        frame  = leg.get_frame()
        frame.set_facecolor('0.50')

    def plot_all_pars(self, diffcutoff=0.2,cmapcutoff=0.1,log=False,plotperc=0):
        totalplotnumber=2000
        result=np.round((self.results-self.worstresult)/(self.bestresult-self.worstresult),3)
        filtermatrix=(result>(1-diffcutoff))
        result=result[filtermatrix]
        pars=self.pars[filtermatrix]
        if log:
            pars=np.log(pars)
        self.cmap=get_scaledhsv_colormap(cscale=cmapcutoff/diffcutoff)
        na=np.array(list(set(result)))
        na.sort()
        pcn=totalplotnumber/(len(na)+3)
        for res in na:
            if res<plotperc:
                continue
            plotpars=pars[result==res]
            if res==1:
                plotpars=plotpars[get_representative_pars(plotpars,cutoffdistance=0.0001,maxcluster=pcn*5)]
            else:
                plotpars=np.random.permutation(plotpars)[:pcn]
            plotpars=plotpars.T
            plt.plot(plotpars,color=self.cmap((res-1+diffcutoff)/diffcutoff),alpha=1,linewidth=1)
        self.currentaxis=plt.gca()
        plt.xlabel('Parameter Index')
        plt.ylabel('Parameter Value')
        plt.title('The Restraints Recover Functions (RRF)')
        ax=matplotlib.colorbar.make_axes(self.currentaxis)
        cb=matplotlib.colorbar.ColorbarBase(ax[0], cmap=self.cmap)
        tp=np.array([0]+list(np.linspace(1-cmapcutoff/diffcutoff,1,10)))
        cb.set_ticks(tp)
        labels=np.round(((tp-1)*diffcutoff+1)*(self.results.max()-self.results.min())+self.results.min(),3)
        if (labels<0).all():
            labels=-labels
        cb.set_ticklabels(labels)
        plt.ylabel('The performace of the RRF')

    def plot_cluster_pars(self,linecolor='g',alphascale=5,log=False):
        totalnumber=2000
        pcn=totalnumber/(len(self.clusters)+4)
        rpcn=pcn
        k=0
        for cl in self.clusters:
            k=k+1
            pars=cl.pars
            if log:
                pars=np.log(pars)
            if k==len(self.clusters):
                rpcn=pcn*5
            pars=pars[get_representative_pars(pars,cutoffdistance=0.0002,maxcluster=rpcn)]
            plt.plot(pars.T,color=linecolor[0:3],linewidth=1)

    def clustering_bestpars(self):
        self.clusterindexes=self.clustering(self.bestpars,self.bestparsindex)
        self.get_cluster_objs(self.clusterindexes)

    def clustering(self,pars,parsindex):
        if len(pars)<3:
            return [parsindex]
        da=scipy.spatial.distance.pdist(pars)
        la=scipy.cluster.hierarchy.linkage(da)
        sp=-min(100,len(la)-1)
        la100=la[sp:,2]
        lar99=la100[1:]/la100[:-1]
        ct=1000
        for i in range(-1,sp,-1):
            if lar99[i]>self.clustercutoff:
                ct=la100[i-1]*self.clustercutoff
                break
        ca=scipy.cluster.hierarchy.fcluster(la,ct,criterion='distance')
        noc=ca.max()
        rac=[]
        for i in range(1,noc+1):
            rac.append(parsindex[np.nonzero(ca==i)[0]])
        return rac

    def get_cluster_objs(self,cil):
        self.clusters=[]
        for ci in cil:
            self.clusters.append(singlecluster(ci,self))

    def dump(self,path):
        self.optscorer=[]
        self.testscorer=[]
        with open(path,'wb') as fh:
            pickle.dump(self,fh)

    def plot_for_andrej1(self,ssp):
        #get the top and bottom curve
        scale=1
        minb=self.bestpars.min(axis=0)
        maxb=self.bestpars.max(axis=0)
        distance=((np.abs(self.bestpars-minb)).sum(axis=1))
        parindex=np.argsort(distance)[0]
        minpar=self.bestpars[parindex]
        distance=((np.abs(self.bestpars-maxb)).sum(axis=1))
        parindex=np.argsort(distance)[0]
        maxpar=self.bestpars[parindex]
        ph1=plt.figure(figsize=[12*scale,8*scale],dpi=90)
        #plot the two reference states
        ssp.plot_dist_for_andrej(np.log(minpar),np.log(maxpar))
        #save the plot
        ph1.savefig('fa1.eps')

    def plot_for_andrej2(self,refs,lt):
        scale=1
        ph1=plt.figure(figsize=[12*scale,8*scale],dpi=90)
        #plot the two reference states
        self.plot_all_pars(plotperc=0.999)
        linestyle=['--','-.']
        ph=[]
        plt.axes(self.currentaxis)
        for i in range(len(refs)):
            ref=refs[i]
            rawresult=self.optscorer.assess_rrf(ref)
            result=np.round((rawresult-self.worstresult)/(self.bestresult-self.worstresult),3)
            #ph.append(plt.plot(np.exp(ref),linestyle[i],alpha=0.9,linewidth=3))
            ph.append(plt.plot(np.exp(ref),linestyle=linestyle[i], color=self.cmap((result-1+0.2)/0.2),alpha=0.8,linewidth=3))#
        plt.legend(ph,lt)
        #save the plotcl
        ph1.savefig('fa2.eps')

class cvclustering(object):
    """
    """
    def __init__(self,optparlist=[],optresultlist=[],scorerlist=[],clustermethod=[], figpath=''):
        self.parlist=optparlist
        self.resultlist=optresultlist
        self.optclusterlist=[]
        self.scorerlist=scorerlist
        self.clustermethod=clustermethod
        self.figpath=figpath

    def analyze(self):
        for i in range(0,len(self.parlist)):
            self.optclusterlist.append(optclustering(optpars=self.parlist[i],optresults=self.resultlist[i], scorer=self.scorerlist[i],clustermethod=self.clustermethod))
        i=0
        for optcluster in self.optclusterlist:
            i=i+1
            print("Analyzing cluster #"+str(i))
            optcluster.analyze()

    def find_overlap(self):
        region=self.optclusterlist[-1].get_pars_region()
        for optcluster in self.optclusterlist:
            region=optcluster.clustering(region)
            if region==[]:
                print("No overlap regions")
                break
        return region

    def plot(self):
        print("ploting...... (can take a while)")
        ph=self.plot_all_bestpars()
        ph.savefig(self.figpath+'all_bestpars.eps')
        ph2,ph3=self.optclusterlist[-1].plot()
        if ph2!=0:
            ph2.savefig(self.figpath+'finalcluster.eps')
        if ph3!=0:
            ph3.savefig(self.figpath+'final_dist_vs_performance.eps')

    def plot_all_bestpars(self):
        scale=1
        ph=plt.figure(figsize=[12*scale,8*scale],dpi=90)
        cdict=matplotlib.pylab.cm.datad['hsv_r']
        cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap2',cdict,256)
        for i in range(len(self.optclusterlist)-1):
            self.optclusterlist[i].plot_cluster_pars(linecolor=cmap(float(i)/(len(self.optclusterlist)-1)))
        self.optclusterlist[i+1].plot_cluster_pars(linecolor=(0,0,0))
        plt.xlabel('Parameter Index')
        plt.ylabel('Parameter Value')
        plt.title('The Restraints Recover Functions (RRF)')
        ax=matplotlib.colorbar.make_axes(plt.gca())
        matplotlib.colorbar.ColorbarBase(ax[0], cmap=cmap)
        plt.ylabel('Different cross validation trials')
        return ph

    def find_concensus_pars(self):
        #find the pars that performs best on all different cv optimization runs
        #Meanwhile we should also notice whether there is large variations from different runs.
        #The fitting pars from the last optimization will play a more importatnt roles here.
        pass

    def dump(self,path):
        for cl in self.optclusterlist:
            cl.optscorer=[]
            cl.testscorer=[]
        self.scorerlist=[]
        with open(path,'wb') as fh:
            pickle.dump(self,fh)

def get_representative_pars(pars, cutoffdistance=0.001, maxcluster=0):
    #print("orinial number: "+str(len(pars)))
    if pars.shape[0]<10:
        return range(len(pars))
    pars=np.round(pars,4)
    da=scipy.spatial.distance.pdist(pars)
    #sda=scipy.spatial.distance.squareform(da)
    la=scipy.cluster.hierarchy.linkage(da,method='complete')
    cutoffdistance=cutoffdistance*pars.shape[1]
    if maxcluster>0:
        if len(la)>maxcluster+1 and la[-maxcluster,2]>cutoffdistance:
            cutoffdistance=la[-maxcluster,2]
    ca=scipy.cluster.hierarchy.fcluster(la,cutoffdistance,criterion='distance')
    noc=ca.max()
    rac=[]
    for i in range(1,noc+1):
        rl=np.nonzero(ca==i)[0]
        rac.append(rl[0])
    #print("representative par number: "+str(len(rac))+' at distance '+str(cutoffdistance))
    return rac

def get_representative_pars_maxlength(pars, cutoffdistance=0.001, maxlength=5000,maxcluster=5000):
    if pars.shape[0]<maxlength:
        return pars[get_representative_pars(pars[:,:-1],cutoffdistance,maxcluster)]
    parshape=pars.shape
    takepercentage=float(maxcluster)/parshape[0]
    numofrun=int(parshape[0]/float(maxlength)+0.8)
    prn=parshape[0]/numofrun+1
    pars=np.random.permutation(pars)
    nal=[]
    apl=range(0,parshape[0],prn)+[parshape[0]]
    print("total orinial number: "+str(parshape[0]))
    for i in range(0,numofrun):
        spars=pars[apl[i]:apl[i+1]]
        nal.append(spars[get_representative_pars(spars[:,:-1],cutoffdistance,maxcluster=int(takepercentage*(apl[i+1]-apl[i])))])
    na=np.vstack(nal)
    print("total representative par number: "+str(len(na))+' at distance '+str(cutoffdistance)+' max '+str(maxcluster))
    return na

def get_representative_pars_maxcluster(pars, cutoffdistance=0.001, maxlength=5000,maxcluster=5000):
    if pars.shape[0]<maxlength:
        return pars[get_representative_pars(pars[:,:-1],cutoffdistance,maxcluster)]
    parshape=pars.shape
    while pars.shape[0]>maxcluster*1.2:
        maxcluster0=max(int(pars.shape[0]*0.8),maxcluster*1.2)
        pars=get_representative_pars_maxlength(pars,cutoffdistance,maxlength,maxcluster0)
    pars=pars[get_representative_pars(pars[:,:-1],cutoffdistance,maxcluster=maxcluster)]
    return pars

def get_representative_pars_forall(na,distanceratio=0,maxcluster=5000,cutoffdistance=0.0001): #obselete
    #ar,apars=self.get_rrf(na)
    ar=na[:,-1]
    pars=na[:,:-1]
    parlen=pars.shape[1]
    parindex=np.arange(len(ar))
    result=np.round((ar-ar.min())/(ar.max()-ar.min()),3)
    rna=np.array(list(set(result)))
    rna.sort()
    pl=[]
    for res in rna:
        fm=(result==res)
        spars=pars[fm]
        pr=result[fm]
        sparindex=parindex[fm]
        print("reducing for performance pars "+str(res))
        na=get_representative_pars_maxcluster(spars,cutoffdistance=(cutoffdistance+(1-res)*distanceratio),maxcluster=maxcluster)
        pl.append(np.hstack((na,ar.min()+res*(ar.max()-ar.min())*np.ones([na.shape[0],1]))))
    return np.vstack(pl)

def get_representative_pars_forbest(na,maxcluster=5000,cutoffdistance=0.0001,clusteringperc=0.999, othermaxnumber=10000):
    #ar,apars=self.get_rrf(na)
    bestresult=na[:,-1].max()
    worstresult=na[:,-1].min()
    if bestresult==worstresult:
        return na[[0,-1],:]
    naperc=(na[:,-1]-worstresult)/(bestresult-worstresult)
    bestna=na[naperc>=clusteringperc]
    otherna=na[naperc<clusteringperc]
    bestna=get_representative_pars_maxcluster(bestna, cutoffdistance,5000,maxcluster)
    if othermaxnumber>0:
        stepsize=int(otherna.shape[0]/othermaxnumber)+1
        otherna=otherna[range(0,otherna.shape[0],stepsize)]
        return np.vstack([bestna,otherna])
    else:
        return bestna
