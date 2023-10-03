from __future__ import print_function
from SOAP import *
import os
import profile
import sys
import pdb
import numpy as np

#recovery function

ml=[]
rfeature=[['d',0,6,0.2]]

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



#different step method

ni=40
ssl={}

initvalues=list(np.arange(0,ni)/10.0+1)

initvalues=np.array(initvalues)

bmtype={'type':'dsscore','dslist':['mhc2'],'criteria':'dcg','combine':'scoresum','bm':'cs1'}


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
    
        scaledsp1={'type':'scaledsp','pdbset':'X2_2.2A_0.25rfree_60','features':sfeature,'genmethod':'cs1','pm':pm,'refs':[ref1],'ratio':ref1['ratio'],'bm':sbm}
        
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

#sampling and searching parameters

inner={'improved':2,'maxloop':10,'minloop':1}
outer={'improved':4,'maxloop':200}

#inner={'improved':1,'maxloop':2,'minloop':1}
#outer={'improved':2,'maxloop':2}

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
#sml=[ssm1]
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

dsearch21={'object':par,'key':'uniform','valueset':[8,5]}



dsearch32={'object':rfeature[0],'key':2,'valueset':[15, 12,10,8,6]}#,10,8,7,6,5#15,12,10, 8,7,15,12,10, 8,7,6,

dsearch34={'object':par,'key':'uniform','valueset':[7,5,4]} #7,6,5,4,3

dsearch35={'object':rfeature[0],'key':3,'valueset':[0.05,0.1,0.5]}#0.025,0.05 #,0.25,0.1,0.05 #0.05,0.1,0.25,

dsearch39={'object':[scaledsp1,scaledsp1],'key':['genmethod','pdbset'],'valueset':[['cs1','X2_2.2A_0.25rfree'],['bs20dsp','X_2.2A_0.25rfree'],['bs10dsp','X_2.2A_0.25rfree'],['bs2dsp','X_2.2A_0.25rfree'],['bs5dsp','X_2.2A_0.25rfree']]}#,,'cs1'['bs20dsp','X_2.2A_0.25rfree'],



dsearch7={'object':pm,'key':0,'valueset':['','pp5','pp1','ks0.2']}#'gps0.6','gps0.3',,'ks0.2','ks0.4'

dsearch8={'object':pm,'key':1,'valueset':['nbsum']}


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


scorerlist=scorerlist+[scaledsp2,ref2]

searchlist=ssl['full'][1]+[search2,search3]
blocks.append([kk,kk+1])

model1={'scorers':scorerlist,'bmtype':bmtype,'searches':searchlist,#'thread':4,
        'dsearches':[dsearch32,dsearch34,dsearch35,dsearch39],'sml':sml,'cvk':2,'fold':3,'repeat':1,'initialpars':'all','runsperscorer':40}



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
elif 0:
    spl=spss(model=model1)
    spl.eval_allpars()
    pdb.set_trace()

spl=spss(model=model1) # set of statistical potentials/discrete searches

spl.currentcutoffpercinitial=0.0# no randomness on discrete searches
spl.currentcutoffpercratio=0.0
spl.maximumround=10
#spl.eval_allpars()
spl.find_best_par()#find the best parameters on discrete and continuous variables.
spl.log()
pdb.set_trace()                
