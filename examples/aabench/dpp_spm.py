from __future__ import print_function
from SOAP import *
import os
import profile
import sys
import pdb
import numpy as np

#runenv.queues=' -q lab.q '

ml=[]
rfeature=[['d',0,15,0.5]]

par={'object':rfeature, 'index1':0,'index2':2,'typelist':[['fixbin',0,1],['fixbin',1,2],['fixbin',1,3],['numafter',[1,3],2]],'parindex':0}

sfeature=[rfeature[0],'a158','as158']
pm=['pp5','','npendunsym']


slo={'key':'sacnf13',
    'valueset':{'sacnf13':{'sftype':'sacnf','par':[1.75,2.75,3.75,4.75,5.75,6.75,7.75,8.75,9.75,10.75,11.75,12.75,13.75],'parvalue':np.ones(15)},
                        'sacnf7':{'sftype':'sacnf','par':[2.25,3.75,5.25,6.75,8.75,10.75,12.75],'parvalue':np.ones(9)},
                        'sacnf5':{'sftype':'sacnf','par':[2.75,4.75,7.25,9.75,12.25],'parvalue':np.ones(7)},
                            #'psn':{'sftype':'psn','par':range(1,15,2),'parvalue':range(7)},
                        }}

slo={'key':'psn',
    'valueset':{'spline16':{'sftype':'pspline','par':range(16),'parvalue':np.ones(16)},
                        'splinen15':{'sftype':'psn','par':range(15),'parvalue':np.ones(15)},
                            'psn':{'sftype':'psn','par':range(1,15,4),'parvalue':range(4)},
                        }} 
par2={'uniform':2,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par3={'uniform':3,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par4={'uniform':4,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par5={'uniform':5,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par6={'uniform':6,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par7={'uniform':7,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par8={'uniform':8,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
slo2={'key':'sacnfflex3',
    'valueset':{
                    'sacnfflex3':{'sftype':'sacnf','par':par3,'parvalue':[1,1,1,1,1]},
                    'sacnfflex4':{'sftype':'sacnf','par':par4,'parvalue':[1,1,1,1,1,1]},
                    'sacnfflex5':{'sftype':'sacnf','par':par5,'parvalue':[1,1,1,1,1,1,1]},
                    'sacnfflex6':{'sftype':'sacnf','par':par6,'parvalue':[1,1,1,1,1,1,1,1]},
                    'sacnfflex7':{'sftype':'sacnf','par':par7,'parvalue':[1,1,1,1,1,1,1,1,1]},
                    'sacnfflex8':{'sftype':'sacnf','par':par8,'parvalue':[1,1,1,1,1,1,1,1,1,1]}}}
# just search parvalue for 7 anchor points with mcb2.0, mxmp1.0 seems really goooooooood - 30638.
#different step method

#different step method

ni=40
ssl={}

initvalues=list(np.arange(0,ni)/10.0+1)

initvalues=np.array(initvalues)
#,'zdock' 'top1000_rlrank__rmsd10ORirmsd4FIRST+top1000_rlrank__rmsd5ORirmsd2FIRST'
#top10__nonempty_rmsd10ORirmsd4FIRST
bmtype={'type':'dsscore','dslist':['aabench'],'criteria':'top1000_rlrank__rmsd10ORirmsd4FIRST+top1000_rlrank__rmsd5ORirmsd2FIRST',
    'combine':'scoresum','bm':'cs1','filters':[{'criteria':'has','criteriadict':{'rmsd':10,'irmsd':4}}],'finalcriteria':'top10__nonempty_rmsd10ORirmsd4'}#,{'criteria':'rep500'}
blocks=[]
kk=0
for apn in [36]:#57
    print(apn)
    noc=np.load(runenv.libdir+'ap'+str(apn)+'.npy').max()+1
    print(noc)
    scorers=[]
    searches=[]
    for i in range(noc):
        print(i)
        sbm=['cs1','ap'+str(apn)+'n',str(i)]
        ref1={'type':'sf','features':sfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0+0.00001*i],'bm':sbm}
    
        scaledsp1={'type':'scaledsp','pdbset':'X_2.2A','features':sfeature,'genmethod':'bs10dsp','pm':pm,'refs':[ref1],'ratio':ref1['ratio'],'bm':sbm}
        
        search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3,4,5,6],'InitialGenerator':{'type':'dfire','values':initvalues}}
        search2={'object':ref1,'key':'par','pos':[0,1,2,3,4,5,6],'InitialGenerator':{'type':'dfire','values':initvalues}}
        search3={'object':ref1,'key':'ratio','pos':[0],'InitialGenerator':[[1.0+0.00001*i] for j in range(40)]}
        scorers.append(ref1)
        scorers.append(scaledsp1)
        searches.append(search1)
        #searches.append(search2)
        searches.append(search3)
        blocks.append([kk,kk+1])
        kk=kk+2
            
    ssl['full']=[scorers,searches]
    #define the final mode


inner={'improved':2,'maxloop':10,'minloop':1}
outer={'improved':4,'maxloop':200}

inner={'improved':1,'maxloop':2,'minloop':1}
outer={'improved':2,'maxloop':2}

td=list(np.sqrt(np.linspace(1,float(10)**2,ni)))
td=np.array(td)

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

ssm0a={'sm':'mca','reset_betweenruns':2,'blockupdate':False, 'using_globalbest':False,
      'sample_schedule':sampleschedule,
      'stepmethod':'mxmp2.0','tune_interval':201,'add_bias':False,'temperature_distribution':np.zeros(ni)+2}

ssm1s={'sm':'mcs','reset_betweenruns':2,'blockupdate':False, 'using_globalbest':True,
      'sample_schedule':sampleschedule,
      'stepmethod':'mxmp2.0','tune_interval':201,'add_bias':False,'temperature_distribution':np.zeros(ni)+2}

sml=[ssm20,ssm2,ssm0,ssm,ssm0,ssm,ssm0,ssm, ssm1]
    


dsearch6={'object':ssm,'key':'temperature_distribution','valueset':[np.sqrt(np.linspace(1,float(10)**2,ni)),
                                                                    np.sqrt(np.linspace(1,float(20)**2,ni)),
                                                                    np.sqrt(np.linspace(1,float(30)**2,ni)),
                                                                    np.sqrt(np.linspace(2,float(20)**2,ni)),
                                                                    np.sqrt(np.linspace(3,float(30)**2,ni))]}

dsearch5={'object':ssm,'key':'add_bias','valueset':[False]}

dsearch4={'object':ssm,'key':'reset_betweenruns','valueset':[2.0]}#,1.56,1.57,1.58,1.59

dsearch3={'object':ssm,'key':'stepmethod','valueset':['mxmp2.0']}

dsearch2={'object':ssm,'key':'sm','valueset':['mcp']}

ssmodel={'index':'full','valueset':ssl}

dsearch1={'object':ssmodel,'key':'index','valueset':['full']}

dsearch6={'object':slo,'key':'key','valueset':['sacnfflex4','sacnfflex7']}#['sacnfflex4']}
dsearch8={'object':inner,'key':'improved','valueset':[1,2,3]}
dsearch9={'object':inner,'key':'maxloop','valueset':[10,30,100]}
dsearch10={'object':inner,'key':'minloop','valueset':[2,5,7,10]}

dsearch11={'object':outer,'key':'improved','valueset':[4,6,8]}

dsearch12={'object':ssm0,'key':'sm','valueset':['mcs','mca']}

dsearch13={'object':ssm0,'key':'using_globalbest','valueset':[True,False]}

dsearch14={'object':ssm0,'key':'temperature_distribution','valueset':[np.zeros(2*ni)+2,td]}

dsearch15={'object':ssm,'key':'exchange_method','valueset':[1]}

dsearch16={'object':[ssm2,ssm20],'key':['blockupdate','blockupdate'],
    'valueset':[[True,True],[False,False]]}

dsearch17={'object':[ssm2,ssm20,ssm,ssm0,ssm1],'key':['blockupdate','blockupdate','blockupdate','blockupdate','blockupdate'],
    'valueset':[[True,True,True,True,True],[False,False,False,False,False],[True,True,False,False,False]]}


scorerlist=ssl['full'][0]

#scorerlist=bm['scorers']
rfeature=[['b',0,1.0,0.1]]


sfeature=[rfeature[0],'a158']
pm=['pp5','','npend']
sratio=[3.8]
scaledsp2={'type':'scaledsp','pdbset':'X2_2.2A','features':sfeature,'genmethod':'','pm':'nbsum','refs':[],'ratio':sratio}

ref2={'type':'sf','features':rfeature,'sftype':'bins','par':range(10),'parvalue':[1 for i in range(10)],'bm':'','ratio':sratio}

search2={'object':ref2,'key':'parvalue','pos':range(0,10),'InitialGenerator':{'type':'dfire','values':list(np.zeros(40))}}

search3={'object':scaledsp2,'key':'ratio','pos':[0],'InitialGenerator':[[3.8] for i in range(40)]}

#dsearch15={'object':sratio,'key':0,'valueset':np.arange(3.5,4.5,0.1)}

scorerlist=scorerlist+[scaledsp2,ref2]

searchlist=ssl['full'][1]+[search2,search3]
blocks.append([kk,kk+1])
#scorerlist=scorerlist+[scaledsp2,ref2]

#searchlist=[search2]
#blocks=[kk,kk+1]

model1={'scorers':scorerlist,'bmtype':bmtype,'searches':searchlist,'thread':4,
        'dsearches':[dsearch6],'sml':sml,'cvk':0,'fold':2,'repeat':3,'initialpars':'all','runsperscorer':40,'testperc':0.0}#,
#dsearch1,searches,dsearch4,dsearch6,dsearch7,dsearch8, ,'testperc':0.333

if 1:
    spl=spss(model=model1)

    spl.currentcutoffpercinitial=0.0
    spl.currentcutoffpercratio=0.0
    spl.maximumround=200
    #spl.eval_allpars()
    spl.find_best_par()
    spl.log()
    pdb.set_trace()


elif 1:

    nm=convert2old(model1)


    for iii in range(3):
        for block in blocks:#[blocks[2],blocks[-1],blocks[3],blocks[5],blocks[10],blocks[8]]+blocks:#blocks:#
            #if k<16:
            #    continue
            #k+=1
            print(block)
            nm['searches']=[]
            if 'str' in nm:
                del nm['str']
            for i in block:
                nm['searches'].append(allsearches[i])
            spso=sps(modellist=[nm])
            spso.cv()
            if len(spso.logpathlist)>0:
                logpath=spso.logpathlist[0][-1]
            else:
                logpath=spso.cvlist[0].logpath
            bestmodel=pickle.load(open(os.path.join(logpath,'bestmodel.pickle')))
            k=-1
            for i in block:
                k=k+1
                try:
                    for j in range(len(allsearches[i]['object'][allsearches[i]['key']])):
                        allsearches[i]['object'][allsearches[i]['key']][j]=bestmodel['searches'][k]['object'][allsearches[i]['key']][j]
                    #use previous run result as initial conditions
                    #allsearches[i]['InitialGenerator']=[bestmodel['searches'][k]['object'][allsearches[i]['key']]]
                except Exception as e:
                    print(e)
                    pdb.set_trace()
