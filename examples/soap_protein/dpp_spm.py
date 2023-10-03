from __future__ import print_function
from sps import *
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
pm=['','npend']

par={'uniform':5,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
slo={'key':'sacnfflex',
    'valueset':{
                    'sacnfflex':{'sftype':'sacnf','par':par,'parvalue':[1,1,1,1]},
                        'psn3':{'sftype':'psn','par':range(1,15,5),'parvalue':range(3)},
                        'psn4':{'sftype':'psn','par':range(1,15,4),'parvalue':range(4)},
                        'psn5':{'sftype':'psn','par':range(1,15,3),'parvalue':range(5)},
                    'psn8':{'sftype':'psn','par':range(0,15,2),'parvalue':range(8)},
                        }} 


# just search parvalue for 7 anchor points with mcb2.0, mxmp1.0 seems really goooooooood - 30638.
#different step method

#different step method

ni=40
ssl={}

initvalues=list(np.arange(0,ni)/10.0+1)

initvalues=np.array(initvalues)
#,'zdock' 'top1000_rlrank__rmsd10ORirmsd4FIRST+top1000_rlrank__rmsd5ORirmsd2FIRST'
#top10__nonempty_rmsd10ORirmsd4FIRST
bmtype={'type':'dsscore','dslist':['DecoysRUS','casp58','baker1'],'criteria':'3xtop1_rmsd_mean_+10xNativeSelection+2xtop3_rmsd_mean_','combine':'scoresum','bm':'bs10dsp'}

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
        sbm=['bs10dsp','ap'+str(apn)+'n',str(i)]
        ref1={'type':'sf','features':sfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0+0.00001*i],'bm':sbm}
    
        scaledsp1={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree','features':sfeature,'genmethod':'bs10dsp','pm':pm,'refs':[ref1],'ratio':ref1['ratio'],'bm':sbm}
        
        search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3,4,5,6],'InitialGenerator':{'type':'dfire','values':initvalues}}
        #search2={'object':ref1,'key':'par','pos':[0,1,2,3,4,5,6],'InitialGenerator':{'type':'dfire','values':initvalues}}
        #search3={'object':ref1,'key':'ratio','pos':[0],'InitialGenerator':[[1.0+0.00001*i] for j in range(40)]}
        scorers.append(ref1)
        scorers.append(scaledsp1)
        searches.append(search1)
        #searches.append(search2)
        #searches.append(search3)
        blocks.append([kk,kk+1,kk+2])
        kk=kk+3
            
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


sml=[ssm20,ssm2,ssm0,ssm,ssm0,ssm,ssm0,ssm, ssm1]
sml=[ssm1]
#sml=[ssm]

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

dsearch6={'object':slo,'key':'key','valueset':slo['valueset'].keys()}
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

dsearch21={'object':par,'key':'uniform','valueset':[8,6,5,4,3]}


scorerlist=ssl['full'][0]

#scorerlist=bm['scorers']
rfeature=[['b',0,1.0,0.1]]


sfeature=[rfeature[0],'a158']
pm=['pp5','','npend']
sratio=[3.8]
scaledsp2={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree','features':sfeature,'genmethod':'','pm':'nbsum','refs':[],'ratio':sratio}

ref2={'type':'sf','features':rfeature,'sftype':'bins','par':range(10),'parvalue':[1 for i in range(10)],'bm':'','ratio':sratio}

search2={'object':ref2,'key':'parvalue','pos':range(0,10),'InitialGenerator':{'type':'dfire','values':list(np.zeros(40))}}

search3={'object':scaledsp2,'key':'ratio','pos':[0],'InitialGenerator':[[3.8] for i in range(40)]}

#dsearch15={'object':sratio,'key':0,'valueset':np.arange(3.5,4.5,0.1)}

scorerlist=scorerlist+[scaledsp2,ref2]

searchlist=ssl['full'][1]+[search2,search3]
blocks.append([kk,kk+1])



model1={'scorers':scorerlist,'bmtype':bmtype,'searches':searchlist,#'thread':4,
        'dsearches':[dsearch21],'sml':sml,'cvk':5,'fold':3,'repeat':1,'initialpars':'all','runsperscorer':40}
        #,'testperc':0.000}#,'testperc':0.0
#dsearch1,searches,dsearch4,dsearch6,dsearch7,dsearch8, ,'testperc':0.333
#'cvk':2,'fold':3,'repeat':10

if 0:
    so=scorer(model=convert2old(model1))
    print(so.assess_model())
    opt=optimizer(scorer=so)
    opt.get_initial_value()    
    print(so.assess(opt.initialvalue))
    opt.optimize()
    pdb.set_trace()
    spso=sps(modellist=[convert2old(model1)])
    #pdb.set_trace()
    spso.cv()#spl.eval_allpars()
    logpath=spso.cvlist[0].logpath    

elif 0:
    so=scorer(model=convert2old(model1))
    print(so.assess_model())
    opt=optimizer(scorer=so)
    opt.get_initial_value()    
    optres=opt.optimize()
    print(so.assess(opt.initialvalue))
    spso=sps(modellist=[convert2old(model1)])
    pdb.set_trace()
    spso.cv()#spl.eval_allpars()
    logpath=spso.cvlist[0].logpath
elif 1:
    spl=spss(model=model1)
    spl.eval_allpars()
    pdb.set_trace()


