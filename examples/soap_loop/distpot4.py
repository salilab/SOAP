from sps import *
import os
import profile
import sys
import pdb
import numpy as np

rfeature=[['d',0,20,0.05]]
#'sftype':'psn3','par':range(1,15,3),'parvalue':range(4)
if 0:                    
                    slo={'key':'psn3',
                        'valueset':{
                                            'psn3':{'sftype':'psn3','par':[0.25,2.25,4.25],'parvalue':range(3)},
                                            'psn4':{'sftype':'psn','par':range(1,15,4),'parvalue':range(4)},
                                            'psn5':{'sftype':'psn','par':range(1,15,3),'parvalue':range(5)},
                                        'psn8':{'sftype':'psn','par':range(0,15,2),'parvalue':range(8)},
                                            }} 
                    
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

par={'uniform':5,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
slo={'key':'sacnfflex',
    'valueset':{
                    'sacnfflex':{'sftype':'sacnf','par':par,'parvalue':[1,1,1,1]},
                        'psn3':{'sftype':'psn','par':range(1,15,5),'parvalue':range(3)},
                        'psn4':{'sftype':'psn','par':range(1,15,4),'parvalue':range(4)},
                        'psn5':{'sftype':'psn','par':range(1,15,3),'parvalue':range(5)},
                    'psn8':{'sftype':'psn','par':range(0,15,2),'parvalue':range(8)},
                        }} 

sfeature=[rfeature[0],'a158','as158']
pm=['','npend']

ni=20

initvalues=list(np.arange(0,ni)/10.0+1)

initvalues=np.array(initvalues)

#decoyset(dsname='ppd4s',sourcedir='/bell2/gqdong/rawdecoyfiles/ppd/ppd4s/').update_ppdt_pir()

ref1={'type':'sf','features':rfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0]}

scaledsp1={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree_60','features':sfeature,'genmethod':'bs2dsp','pm':pm,'refs':[ref1],'ratio':[1.0]}

#define the benchmark scheme
looplist=['mgl'+str(i) for i in range(4,21)]
mjloop=['ml'+str(i) for i in range(4,13)]
#dslist=['DecoysRUS','casp58']
bmtype={'type':'dsscore','dslist':looplist[5:8]+mjloop[5:8],'criteria':'3xtop1_rmsd_mean_+10xNativeSelection+top3_rmsd_mean_+0.001xlbslope6-20-wn','combine':'scoresum','bm':'bs2dsp'}
#top10_rmsd10_irmsd4  10xtop10_rmsd10_irmsd4+10xtop10_rmsd5_irmsd2+  +100xtop10_rmsd5_irmsd2+top10_rmsd+top10_irmsd top10_rmsd10_irmsd4+top20_rmsd+
#define the search positions #3xtop1_rmsd_mean_+NativeSelection+2xtop3_rmsd_mean_
search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}
search2={'object':ref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}

dsearch1={'object':[scaledsp1,bmtype],'key':['genmethod','bm'],'valueset':[['bs1dsp','bs1dsp'],['bs2dsp','bs2dsp'],['bs3dsp','bs3dsp'],['bs4dsp','bs4dsp'],['bs10dsp','bs10dsp']]}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'

dsearch2={'object':rfeature[0],'key':2,'valueset':[20,15, 12,10,8,7,6,5]}#,10,8,7,6,5#15,12,10, 8,7,15,12,10, 8,7,6,

#dsearch3={'object':[rfeature[0],sfeature[0]],'key':[3,3],'valueset':[[0.025,0.025],[0.05,0.05],[0.1,0.1],[0.25,0.25],[0.5,0.5]]}
dsearch4={'object':par,'key':'uniform','valueset':[10,8,6,5,4,3]}

#dsearch4={'object':sfeature[0],'key':2,'valueset':[20,17,15,12,10,8,6]}

dsearch5={'object':rfeature[0],'key':3,'valueset':[0.05,0.1,0.2,0.5]}#0.025,0.05 #,0.25,0.1,0.05 #0.05,0.1,0.25,

dsearch6={'object':slo,'key':'key','valueset':slo['valueset'].keys()}#slo['valueset'].keys()'psn5','psn3','psn4','psn8'

dsearch7={'object':pm,'key':0,'valueset':['','pp1','pp5','ks0.2','ks0.1']}#'gps0.6','gps0.3','ks0.2','ks0.4','gps1.0','gps2.0']}

dsearch10={'object':scaledsp1,'key':'pdbset','valueset':['X_2.2A_0.25rfree','X_2.2A_0.25rfree_30','X_2.2A_0.25rfree_60','X_2.2A_0.25rfree_95']}


dsearch21={'object':scaledsp1,'key':'genmethod','valueset':['bs1dsp','bs2dsp','bs3dsp','bs4dsp','bs10dsp']}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'

dsearch22={'object':bmtype,'key':'bm','valueset':['bs1dsp','bs2dsp','bs3dsp','bs4dsp','bs10dsp']}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'


#dsearch8={'object':pm,'key':1,'valueset':['npend','nbsum']}

#dsearch9={'object':bmtype,'key':'bm','valueset':['ss1dsp','bs1dsp','bs2dsp','bs10dsp']}#,

inner={'improved':2,'maxloop':100,'minloop':2}
outer={'improved':4,'maxloop':5026}

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
#sml=[ssm3,ssm1]

#'sftype':'psn3','par':range(1,15,3),'parvalue':range(4)
#dsearch6,dsearch1,dsearch2,dsearch7,dsearch9
#define the final model
#define the final model

#'sftype':'psn3','par':range(1,15,3),'parvalue':range(4)
#dsearch6,dsearch1,dsearch2,dsearch7,dsearch9
#define the final model
#define the final model
model1={'scorers':[scaledsp1,ref1],'bmtype':bmtype,'searches':[search1,search2], 'runsperscorer':ni,
    'dsearches':[ dsearch2, dsearch4,dsearch5,dsearch7,dsearch10,dsearch21,dsearch22],'sml':sml,'cvk':2,'repeat':1,'fold':3}#,dsearch2,dsearch5,dsearch6,dsearch7 #,'testperc':0.33
#dsearch1,dsearch2,dsearch6,dsearch7
#so=scorer(model=convert2old(model1))
#print so.assess_ideal()
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
#spl.eval_allpars()
#pdb.set_trace()

spl.currentcutoffpercinitial=0.0
spl.currentcutoffpercratio=0.0
spl.maximumround=30
spl.find_best_par()
spl.log()

