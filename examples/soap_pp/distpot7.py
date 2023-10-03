from __future__ import print_function
from sps import *
import os
import profile
import sys
import pdb
import numpy as np

rfeature=[['d',0,15,0.5]]
#'sftype':'psn3','par':range(1,15,3),'parvalue':range(4)
runenv.otherrun=False
par2={'uniform':2,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par3={'uniform':3,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par4={'uniform':4,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par5={'uniform':5,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par6={'uniform':6,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par7={'uniform':7,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
par8={'uniform':8,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}


slo={'key':'sacnfflex4',
    'valueset':{
                    'sacnfflex3':{'sftype':'sacnf','par':par3,'parvalue':[1,1,1,1,1]},
                    'sacnfflex4':{'sftype':'sacnf','par':par4,'parvalue':[1,1,1,1,1,1]},
                    'sacnfflex5':{'sftype':'sacnf','par':par5,'parvalue':[1,1,1,1,1,1,1]},
                    'sacnfflex6':{'sftype':'sacnf','par':par6,'parvalue':[1,1,1,1,1,1,1,1]},
                    'sacnfflex7':{'sftype':'sacnf','par':par7,'parvalue':[1,1,1,1,1,1,1,1,1]},
                    'sacnfflex8':{'sftype':'sacnf','par':par8,'parvalue':[1,1,1,1,1,1,1,1,1,1]},
                    'ig':{'sftype':'sacnf','par':[0],'parvalue':[2.0]},
                    'dfire':{'sftype':'sacnf','par':[0],'parvalue':[2.0]},
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

ni=40

initvalues=list(np.arange(0,ni)/10.0+1)

initvalues=np.array(initvalues)

#decoyset(dsname='ppd4s',sourcedir='/bell2/gqdong/rawdecoyfiles/ppd/ppd4s/').update_ppdt_pir()

ref1={'type':'sf','features':rfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0]}

scaledsp1={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree','features':sfeature,'genmethod':'bs9dsp','pm':pm,'refs':[ref1],'ratio':[1.0]}
#define the benchmark scheme
#dslist=['DecoysRUS','casp58']+2xtop3_rmsd_mean_ 
#dslist=['DecoysRUS','casp58']+2xtop3_rmsd_mean_
looplist=['mgl'+str(i) for i in range(4,21)]

bmtype={'type':'dsscore','dslist':['ppd4s'],'criteria':'top1000_rlrank__rmsd10ORirmsd4FIRST+top1000_rlrank__rmsd5ORirmsd2FIRST',
    'combine':'scoresum','bm':'cs1','filters':[{'criteria':'has','criteriadict':{'rmsd':10,'irmsd':4}}],'finalcriteria':'top10__nonempty_rmsd10ORirmsd4'}#,{'criteria':'rep500'}
#top10_rmsd10_irmsd4  10xtop10_rmsd10_irmsd4+10xtop10_rmsd5_irmsd2+  +100xtop10_rmsd5_irmsd2+top10_rmsd+top10_irmsd top10_rmsd10_irmsd4+top20_rmsd+
#define the search positions #3xtop1_rmsd_mean_+NativeSelection+2xtop3_rmsd_mean_
search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}
search2={'object':ref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}

dsearch1={'object':[scaledsp1,bmtype],'key':['genmethod','bm'],'valueset':[['bs1dsp','bs1dsp'],['bs2dsp','bs2dsp'],['bs3dsp','bs3dsp'],['bs4dsp','bs4dsp'],['bs6dsp','bs6dsp'],['bs8dsp','bs8dsp'],['bs9dsp','bs9dsp'],['bs10dsp','bs10dsp'],['ss1','ss1'],['ss2','ss2'],['ss3','ss3'],['ss4','ss4'],['ss9','ss9']]}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'

dsearch2={'object':rfeature[0],'key':2,'valueset':[20,17,15,12,10,8]}#[6,5.4,5,4.8,4.6,4.4,4.2]#15, 7,6,5,4.5,4,3.5,3#,10,8,7,6,5#15,12,10, 8,7,15,12,10, 8,7,6,

dsearch3={'object':[rfeature[0],sfeature[0]],'key':[3,3],'valueset':[[0.05,0.05],[0.1,0.1]]}
#dsearch4={'object':par,'key':'uniform','valueset':[10,9,8,7,6]}

dsearch5={'object':scaledsp1,'key':'genmethod','valueset':['bs12dsp','bs13dsp','bs14dsp','bs15dsp','bs16dsp','bs17dsp','bs20dsp']}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'

#dsearch6={'object':bmtype,'key':'bm','valueset':['bs12dsp','bs13dsp','bs14dsp','bs15dsp','bs16dsp','bs17dsp','bs20dsp']}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'

#dsearch5={'object':rfeature[0],'key':3,'valueset':[0.05,0.1,0.2]}#0.025,0.05 #,0.25,0.1,0.05 #0.05,0.1,0.25,

dsearch6={'object':slo,'key':'key','valueset':slo['valueset'].keys()}#slo['valueset'].keys()'psn5','psn3','psn4','psn8'

#dsearch6={'object':[rfeature[0],slo],'key':[2,'key'],'valueset':[[7,'sacnf73'], [7,'sacnf74'],[7,'sacnf75'],[5,'sacnf52'], [5,'sacnf53'],[6,'sacnf64'], [6,'sacnf63'],[6,'sacnf62'],]}

dsearch7={'object':pm,'key':0,'valueset':['','gps0.1','gps0.2','gps0.3','gps0.6','gps1.0','gps2.0','ks0.1','ks0.2','ks0.3','ks0.6','ks1.0','ks2.0']}

#dsearch10={'object':scaledsp1,'key':'pdbset','valueset':['X_2.1A_0.25rfree','X_2.1A_0.25rfree_30','X_2.1A_0.25rfree_60','X_2.1A_0.25rfree_95']}

dsearch8={'object':pm,'key':1,'valueset':['npend','nbsum']}

dsearch9={'object':[scaledsp1,scaledsp1],'key':['genmethod','pdbset'],'valueset':[['cs1','X2_2.2A_0.25rfree'],['bs9','X_2.2A_0.25rfree'],['bs8','X_2.2A_0.25rfree'],['bs4','X_2.2A_0.25rfree'],['bs3','X_2.2A_0.25rfree'],['bs2','X_2.2A_0.25rfree'],['bs1','X_2.2A_0.25rfree'],['bs9','X_2.2A_0.25rfree_95'],['bs9','X_2.2A_0.25rfree_60'],['bs9','X_2.2A_0.25rfree_30']]}#,,'cs1'

#dsearch9={'object':scaledsp1,'key':'pdbset','valueset':['X_2.2A_0.25rfree','X_2.2A_0.25rfree_95','X_2.2A_0.25rfree_60','X_2.2A_0.25rfree_30']}#,,'cs1'

inner={'improved':2,'maxloop':30,'minloop':2}
outer={'improved':4,'maxloop':5018}

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
#sml=[ssm2,ssm0]

#'sftype':'psn3','par':range(1,15,3),'parvalue':range(4)
#dsearch6,dsearch1,dsearch2,dsearch7,dsearch9
#define the final model
#define the final model
#

rfeature=[['b',0,1.0,0.1]]


sfeature=[rfeature[0],'a158']
pm=['pp5','','npend']
sratio=[3.8]
scaledsp2={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree','features':sfeature,'genmethod':'','pm':'nbsum','bm':'ss1','refs':[],'ratio':sratio}

ref2={'type':'sf','features':rfeature,'sftype':'bins','par':range(10),'parvalue':[1 for i in range(10)],'bm':'','ratio':sratio}

search3={'object':ref2,'key':'parvalue','pos':range(0,10),'InitialGenerator':{'type':'dfire','values':list(np.zeros(40))}}

search4={'object':scaledsp2,'key':'ratio','pos':[0],'InitialGenerator':[[3.8] for i in range(40)]}

#dsearch15={'object':sratio,'key':0,'valueset':np.arange(3.5,4.5,0.1)}




ssl={'distance':[[scaledsp1,ref1],[search1,search2]],'distancefixed':[[scaledsp1,ref1],[search1]],'distanceacc':[[scaledsp1,ref1,scaledsp2,ref2],[search1,search2,search3,search4]]}

ssmodel={'index':'distanceacc','valueset':ssl}

dsearch21={'object':ssmodel,'key':'index','valueset':ssl.keys()}

model1={'scorers':ssmodel,'bmtype':bmtype,'searches':ssmodel, 'runsperscorer':ni,
    'dsearches':[dsearch21,dsearch2,dsearch9],'sml':sml,'cvk':2,'repeat':1,'fold':3}#,dsearch2,dsearch5,dsearch6,dsearch7 #,'testperc':0.33
#dsearch1,dsearch2,dsearch6,dsearch7

#so=scorer(model=convert2old(model1))
#print(so.assess_ideal())
#print(so.assess_worst())
#print(so.assess_sascore(slevel='top1_rmsdbbif_mean_'))
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
spl.maximumround=200
#spl.eval_allpars()
spl.find_best_par()
spl.log()
pdb.set_trace()

spl.currentcutoffpercinitial=1.0
spl.currentcutoffpercratio=0.0
spl.maximumround=1
spl.find_best_par()
spl.log()

