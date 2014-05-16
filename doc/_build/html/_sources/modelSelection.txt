.. SOAP documentation master file, created by
   sphinx-quickstart on Wed May 14 11:06:57 2014.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Select the best model and generate statistical potentials
================================

Define statistical potential model and model space 
---------------------------------
The model and model space are defined using python lists and dictionaries.

Define the recovery function features, distance only, from 0 to 20A with bin size 0.05::
    rfeature=[['d',0,20,0.05]]

Define the spliens used for generating recoery functions. The most important parameter is 'uniform'::
    par={'uniform':4,'featurebins':rfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
    slo={'key':'sacnfflex',
        'valueset':{
                        'sacnfflex':{'sftype':'sacnf','par':par,'parvalue':[1,1,1,1]},
                        'sacnf52':{'sftype':'sacnf','par':[2.75,3.75],'parvalue':[1,1,1,1]},
        }}
    
Define the recovery function::
    ref1={'type':'sf','features':rfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0]}

Define the features for the probablistic table calcualtion, a158 represent residue dependent atom type ::
    sfeature=[rfeature[0],'a158','as158']
    
Define the processing method of the probablistic table, 'npend' means normalized by the last bin value::
    pm=['','npend']

Define the probabilistic table using for scoring::
    scaledsp1={'type':'scaledsp','pdbset':'X2_2.2A_0.25rfree','features':sfeature,'genmethod':'cs1','pm':pm,'refs':[ref1],'ratio':[1.0]}


The benchmark criteria, defines how the statistical potentials are benchmarked::
    bmtype={'type':'dsscore','dslist':['fpdb','fpdu'],'criteria':'3xtop1_rmsdallif_mean_+2xtop3_rmsdallif_mean_','combine':'scoresum','bm':'cs1'}

Parameters to optimize::
    search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}
    search2={'object':ref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'dfire','values':initvalues}}

Discrete search dictionaries, defining the model space - the model variables/options to vary::
    dsearch4={'object':par,'key':'uniform','valueset':[7,6,5,4,3]}
    
    dsearch9={'object':[scaledsp1,scaledsp1],'key':['genmethod','pdbset'],'valueset':[['cs1','X2_2.2A_0.25rfree']]}#,['cs1','X2_2.2A_0.25rfree_30'],['cs1','X2_2.2A_0.25rfree_60'],['cs1','X2_2.2A_0.25rfree_95'],['bs20dsp','X_2.2A_0.25rfree'],['bs20dsp','X_2.2A_0.25rfree_30'],['bs20dsp','X_2.2A_0.25rfree_60'],['bs15dsp','X_2.2A_0.25rfree'],['bs10dsp','X_2.2A_0.25rfree']]}#,,'cs1'


Parameters controling sampling and optimization::
    ni=40
    
    initvalues=list(np.arange(0,ni)/10.0+1)
    
    initvalues=np.array(initvalues)
    
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


Define the final model::
    model1={'scorers':[scaledsp1,ref1],'bmtype':bmtype,'searches':[search1,search2], 'runsperscorer':ni,
        'dsearches':[dsearch2,dsearch5,dsearch4,dsearch7,dsearch8,dsearch9],'sml':sml,'cvk':2,'repeat':1,'fold':3}

.. automodule:: modelSelection
   :members:

