from __future__ import print_function
from SOAP import *
import os
import profile
import sys
import pdb
import numpy as np

rfeature=[['d',0,20,0.05]]
#'sftype':'psn3','par':range(1,15,3),'parvalue':range(4)

par={'uniform':4,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
slo={'key':'sacnfflex',
    'valueset':{
                    'sacnfflex':{'sftype':'sacnf','par':par,'parvalue':[1,1,1,1]},
                    'sacnf52':{'sftype':'sacnf','par':[2.75,3.75],'parvalue':[1,1,1,1]},
                    'sacnf53':{'sftype':'sacnf','par':[1.75,2.75,3.75],'parvalue':[1,1,1,1,1]},
                    'sacnf62':{'sftype':'sacnf','par':[2.25,3.75],'parvalue':[1,1,1,1]},
                    'sacnf63':{'sftype':'sacnf','par':[2.25,3.75,4.75],'parvalue':[1,1,1,1,1]},
                    'sacnf64':{'sftype':'sacnf','par':[1.75,2.75,3.75,4.75],'parvalue':[1,1,1,1,1,1]},
                    'sacnf73':{'sftype':'sacnf','par':[2.25,3.75,5.25],'parvalue':[1,1,1,1,1]},
                    'sacnf74':{'sftype':'sacnf','par':[2.75,3.75,4.75,5.75],'parvalue':[1,1,1,1,1,1]},
                    'sacnf75':{'sftype':'sacnf','par':[1.75,2.75,3.75,4.75,5.75],'parvalue':[1,1,1,1,1,1,1]},
                                        
                    'sacnf84':{'sftype':'sacnf','par':[1.75,3.25,4.75,6.25],'parvalue':[1,1,1,1,1,1]},
                    'sacnf104':{'sftype':'sacnf','par':[1.75,3.75,5.75,7.75],'parvalue':[1,1,1,1,1,1]},                                        
                    'sacnf105':{'sftype':'sacnf','par':[1.75,3.75,5.25,6.75,8.25],'parvalue':[1,1,1,1,1,1,1]},                                        
                    'sacnf124':{'sftype':'sacnf','par':[2.75,4.75,6.75,8.75],'parvalue':[1,1,1,1,1,1]},                                        
                    'sacnf125':{'sftype':'sacnf','par':[1.75,3.75,5.25,7.75,9.75],'parvalue':[1,1,1,1,1,1,1]},  
                    'sacnf154':{'sftype':'sacnf','par':[2.75,5.75,8.75,11.75],'parvalue':[1,1,1,1,1,1]},                                        
                    'sacnf156':{'sftype':'sacnf','par':[2.75,4.75,6.75,8.75,10.75,12.75],'parvalue':[1,1,1,1,1,1,1,1]},                                        
                    'sacnf158':{'sftype':'sacnf','par':[2.75,4.25,5.75,7.25,8.75,10.25,11.75,13.25],'parvalue':[1,1,1,1,1,1,1,1,1,1]}  
                        }} 
if 0:
                    for valuetype in ['aon','aoo','acn','aco','1oo','1cn','0oo']:
                            if valuetype.startswith('a'):
                                gl=['e','o','l','f','p','s']
                            elif valuetype.startswith('1'):
                                gl=['e','o','f','p']
                            elif valuetype.startswith('0'):
                                gl=['l','s']
                            for gentype in gl:
                                splinetype='s'+partype+valuetype+converttype+gentype
                                slo[splinetype+'2']={'sftype':splinetype,'par':[2.25,3.75]}
                                slo[splinetype+'3']={'sftype':splinetype,'par':[2.25,3.75,4.75]}
                                slo[splinetype+'4']={'sftype':splinetype,'par':[1.75,2.75,3.75,4.75]}

sfeature=[rfeature[0],'a158','as158']
pm=['','npend']

ni=3

initvalues=list(np.arange(0,ni)/10.0+1)

initvalues=np.array(initvalues)

#decoyset(dsname='ppd4s',sourcedir='/bell2/gqdong/rawdecoyfiles/ppd/ppd4s/').update_ppdt_pir()

ref1={'type':'sf','features':rfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0]}

scaledsp1={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree','features':sfeature,'genmethod':'bs10dsp','pm':pm,'refs':[ref1],'ratio':[1.0]}
#define the benchmark scheme
#dslist=['DecoysRUS','casp58']
bmtype={'type':'dsscore','dslist':['mhc2'],'criteria':'dcg','combine':'scoresum','bm':'cs1'}
#top10_rmsd10_irmsd4  10xtop10_rmsd10_irmsd4+10xtop10_rmsd5_irmsd2+  +100xtop10_rmsd5_irmsd2+top10_rmsd+top10_irmsd top10_rmsd10_irmsd4+top20_rmsd+
#define the search positions #3xtop1_rmsd_mean_+NativeSelection+2xtop3_rmsd_mean_
search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}
search2={'object':ref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}

#dsearch1={'object':[scaledsp1,bmtype],'key':['scaledsp1','bm'],'valueset':[['bs1dsp','bs1dsp'],['bs2dsp','bs2dsp'],['bs4dsp','bs4dsp'],['bs10dsp','bs10dsp']]}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'

dsearch2={'object':rfeature[0],'key':2,'valueset':[15,12,10,8,7,6,5,4]}#[6,5.4,5,4.8,4.6,4.4,4.2]#15, 7,6,5,4.5,4,3.5,3#,10,8,7,6,5#15,12,10, 8,7,15,12,10, 8,7,6,

#dsearch3={'object':[rfeature[0],sfeature[0]],'key':[3,3],'valueset':[[0.025,0.025],[0.05,0.05],[0.1,0.1],[0.25,0.25],[0.5,0.5]]}
dsearch4={'object':par,'key':'uniform','valueset':[6]}

dsearch5={'object':rfeature[0],'key':3,'valueset':[0.05,0.1,0.5]}#0.025,0.05 #,0.25,0.1,0.05 #0.05,0.1,0.25,

#dsearch6={'object':slo,'key':'key','valueset':slo['valueset'].keys()}#slo['valueset'].keys()'psn5','psn3','psn4','psn8'

#dsearch6={'object':[rfeature[0],slo],'key':[2,'key'],'valueset':[[7,'sacnf73'], [7,'sacnf74'],[7,'sacnf75'],[5,'sacnf52'], [5,'sacnf53'],[6,'sacnf64'], [6,'sacnf63'],[6,'sacnf62'],]}

#dsearch7={'object':pm,'key':0,'valueset':['','gps0.6','gps0.3','ks0.2','ks0.4','gps1.0','gps2.0']}

#dsearch10={'object':scaledsp1,'key':'pdbset','valueset':['X_2.1A_0.25rfree','X_2.1A_0.25rfree_30','X_2.1A_0.25rfree_60','X_2.1A_0.25rfree_95']}

#dsearch8={'object':pm,'key':1,'valueset':['npend','nbsum']}

dsearch9={'object':[scaledsp1,scaledsp1],'key':['genmethod','pdbset'],'valueset':[['bs3dsp','X_2.2A_0.25rfree'],['bs4dsp','X_2.2A_0.25rfree'],
                ['bs2dsp','X_2.2A_0.25rfree'],['bs9dsp','X_2.2A_0.25rfree'],
                ['bs10dsp','X_2.2A_0.25rfree'],['bs11dsp','X_2.2A_0.25rfree'],['cs1','X2_2.2A_0.25rfree']]}#,,'cs1'

inner={'improved':2,'maxloop':100,'minloop':2}
outer={'improved':4,'maxloop':5023}

td=list(np.sqrt(np.linspace(1,float(10)**2,ni)))
td=np.array(td)
tune_interval=200

sampleschedule={'inner':inner,'outer':outer}
ssm={'sm':'mcp','reset_betweenruns':2,'blockupdate':False, 'exchange_method':1,
      'sample_schedule':sampleschedule,
      'stepmethod':'mxmp2.0','tune_interval':tune_interval,'add_bias':False,'temperature_distribution':td}

ssm0={'sm':'mcs','reset_betweenruns':2,'blockupdate':False,'using_globalbest':True,
      'sample_schedule':sampleschedule,
      'stepmethod':'mxmp2.0','tune_interval':tune_interval,'add_bias':False,'temperature_distribution':td}

ssm2={'sm':'mcp','reset_betweenruns':2,'blockupdate':False, 'exchange_method':1,
      'sample_schedule':sampleschedule,
      'stepmethod':'mxmp2.0','tune_interval':tune_interval,'add_bias':False,'temperature_distribution':td}

ssm20={'sm':'mcs','reset_betweenruns':2,'blockupdate':False,'using_globalbest':True,
      'sample_schedule':sampleschedule,
      'stepmethod':'mxmp2.0','tune_interval':tune_interval,'add_bias':False,'temperature_distribution':td}

ssm1={'sm':'mca','reset_betweenruns':2,'blockupdate':False, 'using_globalbest':True,
      'sample_schedule':sampleschedule,
      'stepmethod':'mxmp2.0','tune_interval':tune_interval,'add_bias':False,'temperature_distribution':np.zeros(ni)+2}

ssm3={'sm':'powell','blockupdate':False}

sml=[ssm20,ssm2,ssm0,ssm,ssm0,ssm,ssm0,ssm, ssm1]
sml=[ssm3,ssm1]

#'sftype':'psn3','par':range(1,15,3),'parvalue':range(4)
#dsearch6,dsearch1,dsearch2,dsearch7,dsearch9
#define the final model
#define the final model
model1={'scorers':[scaledsp1,ref1],'bmtype':bmtype,'searches':[search1,search2], 'runsperscorer':ni,
    'dsearches':[dsearch4],'sml':sml,'cvk':0,'repeat':1,'testperc':0.0}#,dsearch2,dsearch5,dsearch6,dsearch7 #,'testperc':0.33
#dsearch1,dsearch2,dsearch6,dsearch7



if 0:
                #sml=[ssm3,ssm1] 
            
            
            #dsearch6,dsearch1,dsearch2,dsearch7,dsearch9
            #define the final model
            #define the final model
            
            
            so=scorer(model=convert2old(model1))
            opt=optimizer(scorer=so)
            opt.get_initial_value()
            optres=opt.optimize()
            #profile.run('opt.optimize()')
            so.assess_model()
            pdb.set_trace()
            
            spl=spss(model=model1)
            #spl.eval_allpars()
            pdb.set_trace()
            
            spl.currentcutoffpercinitial=0.0
            spl.currentcutoffpercratio=0.0
            spl.maximumround=30
            spl.find_best_par()
            spl.log()

#so=scorer(model=convert2old(model1))
#print(so.assess_ideal())
#print(so.assess_worst())

#pdb.set_trace()

#opt=optimizer(scorer=so)
#opt.get_initial_value()
#optres=opt.optimize()
#profile.run('opt.optimize()')
#so.assess_model()
#pdb.set_trace()
#spso=sps(modellist=[convert2old(model1)])
#pdb.set_trace()
#spso.cv()#spl.eval_allpars()

#pdb.set_trace()
spl=spss(model=model1)

spl.currentcutoffpercinitial=0.0
spl.currentcutoffpercratio=0.0
spl.maximumround=10
#spl.eval_allpars()
spl.find_best_par()
spl.log()
pdb.set_trace()

spl.currentcutoffpercinitial=0.0
spl.currentcutoffpercratio=0.0
spl.maximumround=30
spl.find_best_par()
spl.log()

