from SOAP import *
import profile

#define recovery fucntions
rfeature=[['d',0,20,0.05]]
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

#use more recovery function models, not recommended unless you are sure what this means and think this will help
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

#
ni=40
initvalues=list(np.arange(0,ni)/10.0+1)
initvalues=np.array(initvalues)

#scoring table definitions using raw stats from PDB and recovery functions
sfeature=[rfeature[0],'a158','as158']
pm=['','npend']
ref1={'type':'sf','features':rfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0]}
scaledsp1={'type':'scaledsp','pdbset':'X2_2.2A_0.25rfree','features':sfeature,'genmethod':'cs1','pm':pm,'refs':[ref1],'ratio':[1.0]}

#define the benchmark scheme
bmtype={'type':'dsscore','dslist':['mhc2'],'criteria':'dcg','combine':'scoresum','bm':'cs1'}

#define the continuous parameter (model parameter) to be optimzied
search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}
search2={'object':ref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}


#define the models to be searched
dsearch2={'object':rfeature[0],'key':2,'valueset':[20,5.6,5.8,6,7]}#[6,5.4,5,4.8,4.6,4.4,4.2]#15, 7,6,5,4.5,4,3.5,3#,10,8,7,6,5#15,12,10, 8,7,15,12,10, 8,7,6,
dsearch4={'object':par,'key':'uniform','valueset':[7,6,5,4,3]}
dsearch5={'object':rfeature[0],'key':3,'valueset':[0.05]}#0.1,0.025,0.05 #,0.25,0.1,0.05 #0.05,0.1,0.25,
dsearch7={'object':pm,'key':0,'valueset':['','pp5','pp1','ks0.2','ks0.5']}#'gps0.5','gps0.2'
dsearch8={'object':pm,'key':1,'valueset':['npend','nbsum']}
dsearch9={'object':[scaledsp1,scaledsp1],'key':['genmethod','pdbset'],'valueset':[['cs1','X2_2.2A_0.25rfree']]}#,['cs1','X2_2.2A_0.25rfree_30'],['cs1','X2_2.2A_0.25rfree_60'],['cs1','X2_2.2A_0.25rfree_95'],['bs20dsp','X_2.2A_0.25rfree'],['bs20dsp','X_2.2A_0.25rfree_30'],['bs20dsp','X_2.2A_0.25rfree_60'],['bs15dsp','X_2.2A_0.25rfree'],['bs10dsp','X_2.2A_0.25rfree']]}#,,'cs1'


#optimization and sampling parameters
inner={'improved':2,'maxloop':100,'minloop':2}
outer={'improved':4,'maxloop':5023}

td=list(np.sqrt(np.linspace(1,float(10)**2,ni)))
td=np.array(td)
tune_interval=200

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

ssm3={'sm':'powell','blockupdate':False}

sml=[ssm20,ssm2,ssm0,ssm,ssm0,ssm,ssm0,ssm, ssm1]


#define the final model
model1={'scorers':[scaledsp1,ref1],'bmtype':bmtype,'searches':[search1,search2], 'runsperscorer':ni,
    'dsearches':[dsearch2,dsearch5,dsearch4,dsearch7,dsearch8,dsearch9],'sml':sml,'cvk':2,'repeat':1,'fold':3}#,dsearch2,dsearch5,dsearch6,dsearch7 #,'testperc':0.33

#test the model on local host
if 0:
                so=scorer(model=convert2old(model1))
                print so.assess_ideal()
                print so.assess_worst()
                print so.assess_sascore(slevel='top1_rmsdbbif_mean_')  
                pdb.set_trace()

#model selection
spl=spss(model=model1)

#setting up the parameters for model selection
spl.currentcutoffpercinitial=0.0
spl.currentcutoffpercratio=0.0
spl.maximumround=10

#search for the best model
spl.find_best_par()
spl.log()
pdb.set_trace()


