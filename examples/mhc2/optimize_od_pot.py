from SOAP import *
import os
import profile
import sys
import pdb
import numpy as np

rfeature=[['dt',0,6,0.2]]
rrfeature=[['dt',0,6,6]]

runenv.otherrun=False
par={'uniform':10,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
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
ratio1=[1.0]
#decoyset(dsname='ppd4s',sourcedir='/bell2/gqdong/rawdecoyfiles/ppd/ppd4s/').update_ppdt_pir()
#protein-peptide interactions
ref1={'type':'sf','features':rfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':ratio1}
ref2={'type':'sf','features':[rrfeature[0],'g6#180'],'sftype':'spline','par':[0,60,120,180],'parvalue':[1,1,1,1],'ratio':ratio1}
ref3={'type':'sf','features':[rrfeature[0],'gs6#180'],'sftype':'spline','par':[0,60,120,180],'parvalue':[1,1,1,1],'ratio':ratio1}
ref4={'type':'sf','features':[rrfeature[0],'h12#180#-180'],'sftype':'spline','par':[-180,-60,60,180],'parvalue':[1,1,1,1],'ratio':ratio1}



scaledsp1={'type':'scaledsp','pdbset':'X2_2.2A_0.25rfree','features':sfeature,'genmethod':'cs1','pm':pm,'refs':[ref1,ref2,ref3,ref4],'ratio':ratio1}

ratio2=[1.04]

#peptide intra-atom interactions
ref21={'type':'sf','features':'dt30#6','sftype':slo,'par':slo,'parvalue':slo,'ratio':ratio2,'bm':'bs2dsp'}
ref22={'type':'sf','features':'dt1#6g6#180','sftype':'spline','par':[0,60,120,180],'parvalue':[1,1,1,1],'ratio':ratio2,'bm':'bs2dsp'}
ref23={'type':'sf','features':'dt1#6gs6#180','sftype':'spline','par':[0,60,120,180],'parvalue':[1,1,1,1],'ratio':ratio2,'bm':'bs2dsp'}
ref24={'type':'sf','features':'dt1#6h12#180#-180','sftype':'spline','par':[-180,-60,60,180],'parvalue':[1,1,1,1],'ratio':ratio2,'bm':'bs2dsp'}

scaledsp2={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree','features':'dt30#6g6#180gs6#180h12#180#-180t306ts306','genmethod':'bs2dsp','bm':'bs2dsp','pm':pm,'refs':[ref21,ref22,ref23,ref24],'ratio':ratio2}

ref31={'type':'sf','features':'r20','sftype':'bins','par':[],'parvalue':[0]*21,'ratio':[1.0],'bm':''}

#define the benchmark scheme


bmtype={'type':'dsscore','dslist':['mhc2'],'criteria':'dcg','combine':'scoresum','bm':'cs1'}


#top10_rmsd10_irmsd4  10xtop10_rmsd10_irmsd4+10xtop10_rmsd5_irmsd2+  +100xtop10_rmsd5_irmsd2+top10_rmsd+top10_irmsd top10_rmsd10_irmsd4+top20_rmsd+
#top1000_rlrank__rmsdallif1FIRST+top1000_rlrank__rmsdallif2FIRST+top1000_rlrank__rmsdallif3FIRST
#3xtop1_rmsdallif_mean_+2xtop3_rmsdallif_mean_
#define the search positions
search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}
search2={'object':ref1,'key':'par','pos':[0,1,2],'InitialGenerator':{'type':'dfire','values':initvalues}}

#define the search positions
search3={'object':ref2,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}
search4={'object':ref3,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}
search5={'object':ref4,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}

search21={'object':ref21,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}
search22={'object':ref21,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}

#define the search positions
search23={'object':ref22,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}
search24={'object':ref23,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}
search25={'object':ref24,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}
search9={'object':scaledsp2,'key':'ratio','pos':[0],'InitialGenerator':[[0.99] for i in range(40)]}

search31={'object':ref31,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
#search32={'object':ref31,'key':'ratio','pos':[0],'InitialGenerator':[[0.99] for i in range(ni)]}

#discrete search options
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



dsearch1={'object':[scaledsp1,bmtype],'key':['genmethod','bm'],'valueset':[['bs1dsp','bs1dsp'],['bs2dsp','bs2dsp'],['bs4dsp','bs4dsp'],['bs10dsp','bs10dsp']]}#'bs4dsp','bs5dsp','bs9dsp',,'bs2dsp'

dsearch2={'object':rfeature[0],'key':2,'valueset':[15, 12,10,8,6,5]}#,10,8,7,6,5#15,12,10, 8,7,15,12,10, 8,7,6,

#dsearch3={'object':[rfeature[0],sfeature[0]],'key':[3,3],'valueset':[[0.025,0.025],[0.05,0.05],[0.1,0.1],[0.25,0.25],[0.5,0.5]]}

#dsearch4={'object':sfeature[0],'key':2,'valueset':[20,17,15,12,10,8,6]}


dsearch6={'object':slo,'key':'key','valueset':slo['valueset'].keys()}#slo['valueset'].keys()'psn5','psn3','psn4','psn8'

dsearch7={'object':pm,'key':0,'valueset':['','pp5','pp1','ks0.2']}#'gps0.6','gps0.3',,'ks0.2','ks0.4'

dsearch8={'object':pm,'key':1,'valueset':['nbsum']}

dsearch9={'object':[scaledsp1,scaledsp1],'key':['genmethod','pdbset'],'valueset':[['cs1','X2_2.2A_0.25rfree'],['bs20dsp','X_2.2A_0.25rfree']]}#,,'cs1'['bs20dsp','X_2.2A_0.25rfree'],

dsearch2={'object':[rfeature[0],rrfeature[0],rrfeature[0]],'key':[2,2,3],'valueset':[[6,6,6],[8,8,8]]}#,,[10,10,10],[12,12,12][5.4,5.4,5.4],[4.4,4.4,4.4],10,8,7,6,5#15,12,10, 8,7,15,12,10, 8,7,6,
#[6.4,6.4,6.4],[6.2,6.2,6.2],[5.6,5.6,5.6],,[5,5,5]
#dsearch3={'object':[rfeature[0],sfeature[0]],'key':[3,3],'valueset':[[0.025,0.025],[0.05,0.05],[0.1,0.1],[0.25,0.25],[0.5,0.5]]}
dsearch4={'object':par,'key':'uniform','valueset':[10]} #7,6,5,4,3,[12,10,8,7,6,4]

dsearch5={'object':rfeature[0],'key':3,'valueset':[0.2]}#0.025,0.05 #,0.25,0.1,0.05 #0.05,0.1,0.25,

#dsearch5={'object':rfeature[0],'key':3,'valueset':[0.05,0.1,0.5]}#0.025,0.05 #,0.25,0.1,0.05 #0.05,0.1,0.25,

#dsearch9={'object':bmtype,'key':'bm','valueset':['ss1dsp','bs1dsp','bs2dsp','bs10dsp']}#,

#sampling and seaching parameters

inner={'improved':2,'maxloop':11,'minloop':1}
outer={'improved':4,'maxloop':2001}

td=list(np.sqrt(np.linspace(1,float(10)**2,ni)))
td=np.array(td)
tune_interval=202

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

dsearch99={'object':[ssm,ssm2],'key':['blockupdate','blockupdate'],'valueset':[[False]*2]}#,



model1={'scorers':[scaledsp1,ref1,ref2,ref3,ref4,scaledsp2,ref21,ref22,ref23,ref24,ref31],'bmtype':bmtype,'searches':[search1,search2,search3,search4,search5,search21,search22,search23,search24,search25,search9,search31], 'runsperscorer':ni,#number of replicas per scorers
    #'initialmodelpath':inipath,
    'dsearches':[dsearch4],#[dsearch4,dsearch2,dsearch9],,dsearch2,dsearch9
'sml':sml,#search schedule
'cvk':0,#number of repeation of cross validation
'repeat':4,#number of set of replicas
'fold':0,#fold of cross validation
'testperc':0.00}#,dsearch2,dsearch5,dsearch6,dsearch7 #,'testperc':0.33
#

spl=spss(model=model1) # set of statistical potentials/dicrete searches

spl.currentcutoffpercinitial=0.0# no randomness on dicrete searches
spl.currentcutoffpercratio=0.0
spl.maximumround=10
spl.eval_allpars()
#spl.find_best_par()#find the best parameters on discrete and continuous variables.
spl.log()
pdb.set_trace()
