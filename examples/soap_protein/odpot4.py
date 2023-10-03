from __future__ import print_function
from sps import *
import os
import profile
import sys
import pdb
import numpy as np

rfeature=[['dt',0,6,0.2]]
rrfeature=[['dt',0,6,6]]

runenv.otherrun=False
par={'uniform':8,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
slo={'key':'sacnfflex',
    'valueset':{
                    'sacnfflex':{'sftype':'sacnf','par':par,'parvalue':[1,1,1,1]},
                        'psn3':{'sftype':'psn','par':range(1,15,5),'parvalue':range(3)},
                        'psn4':{'sftype':'psn','par':range(1,15,4),'parvalue':range(4)},
                        'psn5':{'sftype':'psn','par':range(1,15,3),'parvalue':range(5)},
                    'psn8':{'sftype':'psn','par':range(0,15,2),'parvalue':range(8)},
                        }} 

sfeature=[rfeature[0],'g6#180','gs6#180','h12#180#-180','t306','ts306']
pm=['','nbsum']

ni=40

initvalues=list(np.arange(0,ni)/10.0+1)

initvalues=np.array(initvalues)

#decoyset(dsname='ppd4s',sourcedir='/bell2/gqdong/rawdecoyfiles/ppd/ppd4s/').update_ppdt_pir()

ref1={'type':'sf','features':rfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0]}
ref2={'type':'sf','features':[rrfeature[0],'g6#180'],'sftype':'spline','par':[0,60,120,180],'parvalue':[1,1,1,1],'ratio':[1.0]}
ref3={'type':'sf','features':[rrfeature[0],'gs6#180'],'sftype':'spline','par':[0,60,120,180],'parvalue':[1,1,1,1],'ratio':[1.0]}
ref4={'type':'sf','features':[rrfeature[0],'h12#180#-180'],'sftype':'spline','par':[-180,-60,60,180],'parvalue':[1,1,1,1],'ratio':[1.0]}

scaledsp1={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree','features':sfeature,'genmethod':'bs20dsp','pm':pm,'refs':[ref1,ref2,ref3,ref4],'ratio':[1.0]}

#define the benchmark scheme
looplist=['mgl'+str(i) for i in range(4,21)]
mjloop=['ml'+str(i) for i in range(4,13)]
#dslist=['DecoysRUS','casp58']
bmtype={'type':'dsscore','dslist':['DecoysRUS','baker1'],'criteria':'top1_rmsd_mean_+10xNativeSelection+cc','combine':'scoresum','bm':'bs20dsp'}
#top10_rmsd10_irmsd4  10xtop10_rmsd10_irmsd4+10xtop10_rmsd5_irmsd2+  +100xtop10_rmsd5_irmsd2+top10_rmsd+top10_irmsd top10_rmsd10_irmsd4+top20_rmsd+

#define the search positions
search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}
search2={'object':ref1,'key':'par','pos':[0,1,2],'InitialGenerator':{'type':'dfire','values':initvalues}}

#define the search positions
search3={'object':ref2,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}
search4={'object':ref3,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}
search5={'object':ref4,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}



dsearch1={'object':[scaledsp1,bmtype],'key':['genmethod','bm'],'valueset':[['bs1','bs1'],['bs2','bs2'],['bs3','bs3']
    ,['bs4','bs4'],['bs5','bs5'],['bs6','bs6'],['bs7','bs7'],['bs8','bs8'],['bs9','bs9'],['bs10','bs10']
    ,['bs10','bs1'],['bs10','bs2'],['bs10','bs3'],['bs10','bs4']]}

dsearch2={'object':ref1,'key':'features','valueset':['dt800#20','dt200#20','dt40#20','dt600#15','dt300#15','dt150#15','dt30#15','dt100#10']}
dsearch3={'object':ref2,'key':'features','valueset':['g24#180','g12#180','g6#180']}
dsearch4={'object':ref3,'key':'features','valueset':['gs24#180','gs12#180','gs6#180']}
dsearch5={'object':ref4,'key':'features','valueset':['h48#360','h24#360','gs12#360']}

dsearch6={'object':ref1,'key':'par','valueset':[[1,2,8,15],[1,3,8,15],[2,3,8,15],range(1,16,2),range(1,16),range(0,16)]}
dsearch7={'object':ref2,'key':'par','valueset':[[0,60,120,180],range(0,181,30),range(0,181,20),range(0,181,10)]}
dsearch8={'object':ref3,'key':'par','valueset':[[0,60,120,180],range(0,181,30),range(0,181,20),range(0,181,10)]}
#dsearch9={'object':ref4,'key':'par','valueset':[[0,120,240,360],range(0,361,60),range(0,361,30),range(0,361,10)]}

dsearch10={'object':scaledsp1,'key':'features','valueset':['dt30#6g6#180gs6#180h12#180#-180t306ts306','dt12#6g6#180gs6#180h12#180#-180t306ts306']}
dsearch11={'object':scaledsp1,'key':'pdbset','valueset':['X_2.2A','X_2.2A_nopsi']}



dsearch1={'object':[scaledsp1,bmtype],'key':['genmethod','bm'],'valueset':['bs20dsp','bs20dsp']}#[['bs2dsp','bs2dsp'],['bs4dsp','bs4dsp'],['bs9dsp','bs9dsp'],['bs10dsp','bs10dsp'],['bs12dsp','bs12dsp']]}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'

dsearch2={'object':[rfeature[0],rrfeature[0],rrfeature[0]],'key':[2,2,3],'valueset':[[7,7,7],[6,6,6]]}#,[5.4,5.4,5.4],[4.4,4.4,4.4],10,8,7,6,5#15,12,10, 8,7,15,12,10, 8,7,6,

#dsearch3={'object':[rfeature[0],sfeature[0]],'key':[3,3],'valueset':[[0.05,0.05],[0.1,0.1],[0.2,0.2]]}
dsearch4={'object':par,'key':'uniform','valueset':[8,6,5,4]}#[10,8,6,5,4,3]}



dsearch6={'object':slo,'key':'key','valueset':slo['valueset'].keys()}#slo['valueset'].keys()'psn5','psn3','psn4','psn8'

dsearch7={'object':pm,'key':0,'valueset':['','ks0.2','pp1','pp5']}#'gps0.6','gps0.3',,'ks0.2','ks0.4'

dsearch8={'object':pm,'key':1,'valueset':['nbsum']}

dsearch9={'object':[scaledsp1,scaledsp1],'key':['genmethod','pdbset'],'valueset':[['bs10dsp','X_2.2A_0.25rfree'],['cs1','X2_2.2A_0.25rfree']]}#,,'cs1'

#dsearch2={'object':[rfeature[0],rrfeature[0],rrfeature[0]],'key':[2,2,3],'valueset':[[6.4,6.4,6.4],[6.2,6.2,6.2],[6,6,6],[5.8,5.8,5.8],[5.6,5.6,5.6],[5.4,5.4,5.4],[5.2,5.2,5.2],[5,5,5]]}#,[5.4,5.4,5.4],[4.4,4.4,4.4],10,8,7,6,5#15,12,10, 8,7,15,12,10, 8,7,6,

#dsearch3={'object':[rfeature[0],sfeature[0]],'key':[3,3],'valueset':[[0.025,0.025],[0.05,0.05],[0.1,0.1],[0.25,0.25],[0.5,0.5]]}
#dsearch4={'object':par,'key':'uniform','valueset':[6,5,4,3]} 

#dsearch5={'object':rfeature[0],'key':3,'valueset':[0.2]}#0.025,0.05 #,0.25,0.1,0.05 #0.05,0.1,0.25,

#dsearch5={'object':rfeature[0],'key':3,'valueset':[0.05,0.1,0.5]}#0.025,0.05 #,0.25,0.1,0.05 #0.05,0.1,0.25,

#dsearch9={'object':bmtype,'key':'bm','valueset':['ss1dsp','bs1dsp','bs2dsp','bs10dsp']}#,

dsearch5={'object':scaledsp1,'key':'genmethod','valueset':['bs20dsp','bs15dsp','bs10dsp','bs30dsp']}#'bs4dsp','bs15dsp','bs5dsp',,'bs30dsp''bs9dsp',,'bs2dsp'

dsearch9={'object':scaledsp1,'key':'pdbset','valueset':['X_2.2A_0.25rfree','X_2.2A_0.25rfree_60','X_2.2A_0.25rfree_30']}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'


dsearch6={'object':bmtype,'key':'bm','valueset':['bs4dsp','bs8dsp','bs9dsp','bs10dsp','bs11dsp','bs12dsp','bs13dsp','bs14dsp','bs15dsp']}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'


inner={'improved':2,'maxloop':11,'minloop':1}
outer={'improved':4,'maxloop':2009}

td=list(np.sqrt(np.linspace(1,float(10)**2,ni)))
td=np.array(td)
tune_interval=200

sampleschedule={'inner':inner,'outer':outer}
ssm={'sm':'mcp','reset_betweenruns':2,'blockupdate':False, 'exchange_method':1,
      'sample_schedule':sampleschedule,
      'stepmethod':'mxmp2.0','tune_interval':tune_interval,'add_bias':False,'temperature_distribution':td}

ssm0={'sm':'mcs','reset_betweenruns':1,'blockupdate':False,'using_globalbest':True,
      'sample_schedule':sampleschedule,
      'stepmethod':'mxmp1.0','tune_interval':tune_interval,'add_bias':False,'temperature_distribution':td}

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
#sml=[ssm0] 

#dsearch6,dsearch1,dsearch2,dsearch7,dsearch9
#define the final model
#define the final model
if 0:
                    sml=[ssm20,ssm2,ssm0,ssm,ssm0,ssm,ssm0,ssm, ssm1]
                    #sml=[ssm3,ssm1] 
                    
                    
                    #dsearch6,dsearch1,dsearch2,dsearch7,dsearch9
                    #define the final model
                    #define the final model
                    
                    model1={'scorers':[scaledsp1,ref1,ref2,ref3,ref4],'bmtype':bmtype,'searches':[search1,search2,search3,search4,search5], 'runsperscorer':ni,
                        'dsearches':[dsearch7],'sml':sml,'cvk':0,'repeat':20,'testperc':0.00}#,dsearch2,dsearch5,dsearch6,dsearch7 #,'testperc':0.33
                    #dsearch1,dsearch2,dsearch6,dsearch7
                    
                    #so=scorer(model=convert2old(model1))
                    #opt=optimizer(scorer=so)
                    #opt.get_initial_value()
                    #optres=opt.optimize()
                    #profile.run('opt.optimize()')
                    #so.assess_model()
                    #pdb.set_trace()
                    
                    spl=spss(model=model1)
                    #spl.eval_allpars()
                    #pdb.set_trace()
                    
                    spl.currentcutoffpercinitial=0.0
                    spl.currentcutoffpercratio=0.0
                    spl.maximumround=30
                    spl.find_best_par()
                    spl.log()



model1={'scorers':[scaledsp1,ref1,ref2,ref3,ref4],'bmtype':bmtype,'searches':[search1,search2,search3,search4,search5], 'runsperscorer':ni,
    #'initialmodelpath':inipath,
    'dsearches':[dsearch5],'sml':sml,'cvk':3,'repeat':1}#,'testperc':0.00}#,dsearch2,dsearch5,dsearch6,dsearch7 #,'testperc':0.33
#dsearch1,dsearch2,dsearch6,dsearch7dsearch4,dsearch2, $#dsearch2,dsearch4,dsearch7,dsearch5

#so=scorer(model=convert2old(model1))
#opt=optimizer(scorer=so)
#opt.get_initial_value()

#optres=opt.optimize()
#profile.run('opt.optimize()')
#so.assess_model()
#pdb.set_trace()

#pdb.set_trace()
spl=spss(model=model1)

spl.currentcutoffpercinitial=0.0
spl.currentcutoffpercratio=0.0
spl.maximumround=10
#spl.eval_allpars()
spl.find_best_par()
spl.log()
pdb.set_trace()

logpath=spso.cvlist[0].logpath
print(logpath)
bm=pickle.load(open(logpath+'cvmodel.pickle', 'rb'))
bestpars=bm['allresult'][0]['repna']
so=scorer(model=bm['model'])
print(so.assess(bm['allresult'][0]['bestpar']))
print(so.assess_model())
bm['model']['bmtype']['dslist']=mjloop
so=scorer(model=bm['model'])
print(so.assess(bm['allresult'][0]['bestpar']))
print(so.assess_model())

pdb.set_trace()



