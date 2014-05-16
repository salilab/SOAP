from sps import *
import cProfile
runenv.otherrun=False
ni=40

##########bmtype section #################

#dslist=['fpdb','fpdu']
#dslist=['DecoysRUS','casp58']
dslist=['ml'+str(i) for i in range(4,13)]#
#looplist=['mgl'+str(i) for i in range(4,21)]
dslist=['DecoysRUS','baker1']
dslist=['wang','wangm','astex','astexm']
dslist=['dude']
criteria='3xtop1_rmsd_mean_+NativeSelection+2xtop3_rmsd_mean_'#'3xtop1_rmsd_mean_+NativeSelection+2xtop3_rmsd_mean_'
criteria='top1_rmsd_mean_+10xNativeSelection+cc'
criteria='enrichment'
bm=['cs1'] 
filters=[]
finalcriteria='enrichment'
if True:
    bmtype={'type':'dsscore','dslist':dslist,'criteria':criteria,'combine':'scoresum','bm':bm}
    if len(filters)>0:
        bmtype['filters']=filters
    if len(finalcriteria)>0:
        bmtype['finalcriteria']=finalcriteria
dsearchbm={'object':bmtype,'key':'bm','valueset':['bs2dsp','bs4dsp','bs9dsp','bs15dsp']}#,'bs20dsp','bs30dsp'

########### Features, searches and dsearches section ###########
ssl={} # scorers searches dictionary, specify which pair to use later

####Distance features for dsearch
dfeature=[['dl',0,7,0.05]]
dsfeature=[dfeature[0],'al184','asl184']
dpm=['pp1','npend']
dpdbset=['ligX']
dgenmethod=['cs1']
dgenmethod=bm
dratio=[1.0]

####svd distance scorers
noe=int(dfeature[0][2]/dfeature[0][3])
svdpm=dpm#['','npend']
if True:#57
    #noc=np.load(runenv.libdir+'svd'+str(apn)+'.npy').max()+1
    dssvd=[]
    svdspl=[]
    for i in range(noe):
        sbm=[bm,'svd'+str(i)]   
        scaledsp={'type':'scaledsp','pdbset':dpdbset,'features':dsfeature,'genmethod':dgenmethod,'pm':svdpm,'ratio':[1.0],'bm':sbm}
        search3={'object':scaledsp,'key':'ratio','pos':[0],'InitialGenerator':[[1.0+0.00001*i] for j in range(40)]}
        dssvd.append(search3)
        svdspl.append(scaledsp)   
ssl['dsvd20']=[svdspl[:50],dssvd[:50]]
ssl['dsvdall']=[svdspl,dssvd]
ssl['dsvd220']=[svdspl[1:50],dssvd[1:50]]
ssl['dsvd210']=[svdspl[1:],dssvd[1:10]]
ssl['dsvd110']=[svdspl[:50],dssvd[:10]]
ssl['dsvd10']=[svdspl,dssvd[0:10]]

####Distance scorers and searches
dSfList={'sacnf':range(4,9,2),'sacof':[4,8],'dfire':None,'rwr':[1,2,3],'normb':[1,2,3],'lognb':[1,2,3],'expb':[1,2,3],'rexpb':[1,2,3],'svdb':[range(0,10)],'pl':[3,4,6]}#'ig','dfire','bins','sinb':[],
dSfList={'sacnf':[4,8,10],'sacof':[4,8],'dfire':None,'expb':[1,3]}
if True:
    dscaledsp1={'type':'scaledsp','pdbset':dpdbset,'features':dsfeature,'genmethod':dgenmethod,'pm':dpm,'ratio':dratio}#'refs':[dref1]

    dsfl={'key':'sacnf10',
        'valueset':{
                            }}
    bnd={'rwr':3,'normb':3,'lognb':3,'expb':4,'rexpb':5}
    if 'bins' in dSfList:
        dsfl['valueset']['bins']={'sftype':'bins','par':[],'parvalue':None}
    for key in ['ig','dfire']:
        if key in dSfList:
            dsfl['valueset'][key]={'sftype':key,'par':[],'parvalue':[2.0]}
    if 'sacnf' in dSfList:
        for nb in dSfList['sacnf']:
            par={'uniform':nb,'featurebins':dfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
            dsfl['valueset']['sacnf'+str(nb)]={'sftype':'sacnf','par':par,'parvalue':[1.01 for i in range(nb+2)]}
    if 'sacf' in dSfList:
        for nb in dSfList['sacnf']:
            par={'uniform':nb,'featurebins':dfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
            dsfl['valueset']['sacf'+str(nb)]={'sftype':'sacf','par':par,'parvalue':[1.01 for i in range(nb+2)]}
    for key in bnd.keys():
        if key in dSfList:
            for nb in dSfList[key]:
                dsfl['valueset'][key+str(nb)]={'sftype':key,'par':[],'parvalue':[1 for i in range(nb*bnd[key])]}
    if 'svdb' in dSfList:
        for pvs in dSfList['svdb']:
            dsfl['valueset']['svdb']={'sftype':'svdb','par':[],'parvalue':pvs,'rawsp':scaledsp}
    if 'pl' in dSfList:
        for order in dSfList['pl']:
            dsfl['valueset']['pl'+str(order)]={'sftype':'pl','par':[],'parvalue':[1 for i in range(order)]}
    dref1={'type':'sf','features':dfeature,'sftype':dsfl,'par':dsfl,'parvalue':dsfl,'ratio':dratio}
    dsParvalue={'object':dref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    dsPar={'object':dref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    dsRatio={'object':dref1,'key':'ratio','pos':[0],'InitialGenerator':{'type':'random'}}

ssl['d']=[[dscaledsp1,dref1],[dsParvalue,dsPar]]
#ssl['d']=[[dscaledsp1],[]]
svdspl=copy.copy(svdspl)
svdspl.extend([dref1])
dssvd=copy.copy(dssvd)
dssvd.extend([dsParvalue,dsPar,dsRatio])
ssl['dsvdref']=[svdspl,dssvd]


####cluatered distance scorers
if False:
    blocks=[]
    kk=0
    for apn in [36]:#57
        print apn
        noc=np.load(runenv.libdir+'ap'+str(apn)+'.npy').max()+1
        print noc
        scorers=[]
        searches=[]
        for i in range(noc):
            print i
            sbm=['cs1','ap'+str(apn)+'n',str(i)]
            reft={'type':'sf','features':sfeature,'sftype':slo,'par':slo,'parvalue':slo,'ratio':[1.0+0.00001*i],'bm':sbm}
        
            scaledspt={'type':'scaledsp','pdbset':'X_2.2A_0.25rfree','features':sfeature,'genmethod':'bs10dsp','pm':pm,'ratio':ref1['ratio'],'bm':sbm}
            
            search1={'object':ref1,'key':'parvalue','pos':[0,1,2,3,4,5,6],'InitialGenerator':{'type':'dfire','values':initvalues}}
            #search2={'object':ref1,'key':'par','pos':[0,1,2,3,4,5,6],'InitialGenerator':{'type':'dfire','values':initvalues}}
            search3={'object':ref1,'key':'ratio','pos':[0],'InitialGenerator':[[1.0+0.00001*i] for j in range(40)]}
            scorers.append(ref1)
            scorers.append(scaledsp1)
            searches.append(search1)
            #searches.append(search2)
            searches.append(search3)
            blocks.append([kk,kk+1,kk+2])
            kk=kk+3
                
        ssl['full']=[scorers,searches]
        #define the final mode

####Distance discrete searches
ddsCutoff={'object':dfeature[0],'key':2,'valueset':[20,15,10,7,6]}
ddsBinsize={'object':dfeature[0],'key':3,'valueset':[0.05,0.1,0.2]}
ddsGenmethod={'object':dgenmethod,'key':0,'valueset':['bs2dsp','bs4dsp','bs9dsp','bs15dsp','bs20dsp','bs30dsp']}
ddsPdbset={'object':dpdbset,'key':0,'valueset':['X_2.2A_0.25rfree','X_2.2A_0.25rfree_95','X_2.2A_0.25rfree_60','X_2.2A_0.25rfree_30']}
ddsPm0={'object':dpm,'key':0,'valueset':['','pp5','pp1']}#'gps0.5','gps0.2'
ddsPm1={'object':dpm,'key':1,'valueset':['npend','nbsum']}#,'nbsum+svd120','nbsum+svd120+exp+npend'
#ddsPdbsetGenmetod={'object':[scaledsp1,scaledsp1],'key':['genmethod','pdbset'],'valueset':[['cs1','X2_2.2A_0.25rfree']]}#,['cs1','X2_2.2A_0.25rfree_30'],['cs1','X2_2.2A_0.25rfree_60'],['cs1','X2_2.2A_0.25rfree_95'],['bs20dsp','X_2.2A_0.25rfree'],['bs20dsp','X_2.2A_0.25rfree_30'],['bs20dsp','X_2.2A_0.25rfree_60'],['bs15dsp','X_2.2A_0.25rfree'],['bs10dsp','X_2.2A_0.25rfree']]}#,,'cs1'
ddsSf={'object':dsfl,'key':'key','valueset':dsfl['valueset'].keys()}


####OD pot scorers and searches
odrfeature=[['dt',0,8,0.2]]
odrrfeature=[['dt',0,8,8]]
odgf=[['g',0,180,30]]
odgsf=[['gs',0,180,30]]
odhf=[['h',-180,180,30]]
odpm=['','nbsum']

if True:
    sfeature=[odrfeature[0],'g6#180','gs6#180','h12#180#-180','t306','ts306']
    odscaledsp1={'type':'scaledsp','pdbset':dpdbset,'features':sfeature,'genmethod':dgenmethod,'pm':odpm,'ratio':dratio}
    odref1={'type':'sf','features':odrfeature,'sftype':dsfl,'par':dsfl,'parvalue':dsfl,'ratio':dratio}
    odref2={'type':'sf','features':[odrrfeature[0],odgf[0]],'sftype':'spline','par':[0,60,120,180],'parvalue':[1,1,1,1],'ratio':dratio}
    odref3={'type':'sf','features':[odrrfeature[0],odgsf[0]],'sftype':'spline','par':[0,60,120,180],'parvalue':[1,1,1,1],'ratio':dratio}
    odref4={'type':'sf','features':[odrrfeature[0],odhf[0]],'sftype':'spline','par':[-180,-60,60,180],'parvalue':[1,1,1,1],'ratio':dratio}
    odsParvalue={'object':odref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    odsPar={'object':odref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    odsRatio={'object':odref1,'key':'ratio','pos':[0],'InitialGenerator':{'type':'random'}}
    odsearch3={'object':odref2,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}
    odsearch4={'object':odref3,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}
    odsearch5={'object':odref4,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'sin','values':np.arange(0,ni)/10.0}}

ssl['od']=[[odscaledsp1,odref1,odref2,odref3,odref4],[odsParvalue,odsPar,odsearch3,odsearch4,odsearch5]]

####OD discrete dsearches

oddsCutoff={'object':[odrfeature[0],odrrfeature[0],odrrfeature[0]],'key':[2,2,3],'valueset':[[item,item,item] for item in [10,8,6]]}
oddsBinsize={'object':odrfeature[0],'key':3,'valueset':[0.2,0.4]}
oddsPm0={'object':odpm,'key':0,'valueset':['','pp5','pp1','ks0.2']}#'gps0.5','gps0.2'
oddsPm1={'object':odpm,'key':1,'valueset':['nbsum']}#,'nbsum+svd120','nbsum+svd120+exp+npend'

####OD-part scorers and searches
odppm=['','nbsumi1']
odpgratio=[1.0]
odphratio=[1.0]

if True:
    sfeature1=[odgf, odrfeature[0],'t306','ts306']
    sfeature2=[odhf, odrfeature[0],'t306','ts306']
    odpscaledsp1={'type':'scaledsp','pdbset':dpdbset,'features':sfeature1,'genmethod':dgenmethod,'pm':odppm,'ratio':odpgratio}
    odpscaledsp2={'type':'scaledsp','pdbset':dpdbset,'features':sfeature2,'genmethod':dgenmethod,'pm':odppm,'ratio':odphratio}
    odpref1={'type':'sf','features':[odrrfeature[0],odgf[0]],'sftype':'spline','par':[0,60,120,180],'parvalue':[1,1,1,1],'ratio':odpgratio}
    odpref2={'type':'sf','features':[odrrfeature[0],odhf[0]],'sftype':'spline','par':[-180,-60,60,180],'parvalue':[1,1,1,1],'ratio':odphratio}
    odpsParvalue1={'object':odpref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    odpsPar1={'object':odpref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    odpsParvalue2={'object':odpref2,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    odpsPar2={'object':odpref2,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    odpsRatio1={'object':odpref1,'key':'ratio','pos':[0],'InitialGenerator':{'type':'random'}}
    odpsRatio2={'object':odpref2,'key':'ratio','pos':[0],'InitialGenerator':{'type':'random'}}

ssl['odp']=[[dscaledsp1,dref1,odpscaledsp1,odpscaledsp2,odpref1,odpref2],[dsParvalue,dsPar,dsRatio,odpsParvalue1,odpsParvalue2,odpsPar1,odpsPar2,odpsRatio1,odpsRatio2]]

####OD-part discrete searches
odpdsPm0={'object':odppm,'key':0,'valueset':['','pp5','pp1']}#'gps0.5','gps0.2'
odpdsCutoff={'object':[odgf[0],odgsf[0],odhf[0]],'key':[3,3,3],'valueset':[[item,item,item] for item in [30,20,10,5]]}

####Fractional Surface assessibility feature
bfeature=[['b',0,1.0,0.1]]
bsfeature=[bfeature[0],'a184']
bpm=['pp1','nbsum']
bratio=[3.8]
dSfList={'sacnf':[4],'sacf':[4],'bins':True,'expb':[1,2]}#'ig','dfire','bins','sinb':[],

if True:
    bscaledsp2={'type':'scaledsp','pdbset':dpdbset,'features':bsfeature,'genmethod':'','pm':bpm,'refs':[],'ratio':bratio}

    bsfl={'key':'expb1',
        'valueset':{
                            }}
    if 'bins' in dSfList:
        bsfl['valueset']['bins']={'sftype':'bins','par':[],'parvalue':None}

    if 'sacnf' in dSfList:
        for nb in dSfList['sacnf']:
            par={'uniform':nb,'featurebins':dfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
            bsfl['valueset']['sacnf'+str(nb)]={'sftype':'sacnf','par':par,'parvalue':[1 for i in range(nb+2)]}
    if 'sacf' in dSfList:
        for nb in dSfList['sacnf']:
            par={'uniform':nb,'featurebins':dfeature[0],'distribute':'lownumberpriority','firstpoint':2.25}
            bsfl['valueset']['sacf'+str(nb)]={'sftype':'sacf','par':par,'parvalue':[1 for i in range(nb+2)]}
    for key in bnd.keys():
        if key in dSfList:
            for nb in dSfList[key]:
                bsfl['valueset'][key+str(nb)]={'sftype':key,'par':[],'parvalue':[1 for i in range(nb*bnd[key])]}
    bref1={'type':'sf','features':dfeature,'sftype':bsfl,'par':bsfl,'parvalue':bsfl,'ratio':bratio}
    bsParvalue={'object':bref1,'key':'parvalue','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    #bsPar={'object':bref1,'key':'par','pos':[0,1,2,3],'InitialGenerator':{'type':'random'}}
    bsRatio={'object':bref1,'key':'ratio','pos':[0],'InitialGenerator':{'type':'random'}}

ssl['d+b']=[[dscaledsp1,dref1,bscaledsp2,bref1],[dsParvalue,dsPar,dsRatio,bsParvalue,bsRatio]]

####surface discrete searches


####distance-ss scorers and searches
dsspm=['pp5','npendi1']

if True:
    dssfeature=[dfeature[0], 'm5','ms5','a158','as158']
    dssScaledsp1={'type':'scaledsp','pdbset':dpdbset,'features':dssfeature,'genmethod':dgenmethod,'pm':dsspm}

ssl['dss']=[[dssScaledsp1,dref1],[dsParvalue,dsPar]]

####distance-ss scorers discrete searches
dssdsPm0={'object':dsspm,'key':0,'valueset':['','pp5','pp1']}#'gps0.5','gps0.2'


####distance-aa scorers and searches
daapm=['pp5','nbsum']

if True:
    daafeature=[dfeature[0], bfeature[0],bfeature[0],'a158','as158']
    daaScaledsp1={'type':'scaledsp','pdbset':dpdbset,'features':daafeature,'genmethod':dgenmethod,'pm':daapm}

ssl['daa']=[[daaScaledsp1,dref1],[dsParvalue,dsPar]]

####distance-aa scorers discrete searches
daadsPm0={'object':daapm,'key':0,'valueset':['','pp5','pp1']}#'gps0.5','gps0.2'

####distance-aass scorers and searches

if True:
    daassfeature=[dfeature[0], 'b2','bs2','m5','ms5','a158','as158']
    daassScaledsp1={'type':'scaledsp','pdbset':dpdbset,'features':daassfeature,'genmethod':dgenmethod,'pm':dsspm}

ssl['daass']=[[daassScaledsp1,dref1,bref1],[dsParvalue,dsPar,dsRatio,bsParvalue,bsRatio]]

####distance-aa scorers discrete searches
daadsPm0={'object':daapm,'key':0,'valueset':['','pp5','pp1']}#'gps0.5','gps0.2'


####features to use and search search

#ssl contains all possible scorers and their parameters to optimize
ssmodel={'index':'d','valueset':ssl}

#dsSc defines the scorers to try
dsSs={'object':ssmodel,'key':'index','valueset':['d','dsvd20','dsvdall','dsvd220','dsvdref','dsvd210','dsvd110','dsvd10']}#'d','b'

#dsList defins the model space to search for
dsList=[ddsCutoff,ddsBinsize,ddsSf]#[dsSs,ddsPm0,ddsPm1]#ddsSf

##########optimization method section ##############
step=1001
if True:
    inner={'improved':2,'maxloop':100,'minloop':2}
    outer={'improved':4,'maxloop':5023}
    
    td=list(np.sqrt(np.linspace(1,float(10)**2,ni)))
    td=np.array(td)
    
    sampleschedule={'inner':inner,'outer':outer}
    ssm={'sm':'mcp','reset_betweenruns':2,'blockupdate':False, 'exchange_method':1,
          'sample_schedule':sampleschedule,
          'stepmethod':'mxmp2.0','tune_interval':step,'add_bias':False,'temperature_distribution':td}
    
    ssm0={'sm':'mcs','reset_betweenruns':2,'blockupdate':False,'using_globalbest':True,
          'sample_schedule':sampleschedule,
          'stepmethod':'mxmp2.0','tune_interval':step,'add_bias':False,'temperature_distribution':td}
    
    
    ssm2={'sm':'mcp','reset_betweenruns':2,'blockupdate':False, 'exchange_method':1,
          'sample_schedule':sampleschedule,
          'stepmethod':'mxmp2.0','tune_interval':step,'add_bias':False,'temperature_distribution':td}
    
    ssm20={'sm':'mcs','reset_betweenruns':2,'blockupdate':False,'using_globalbest':True,
          'sample_schedule':sampleschedule,
          'stepmethod':'mxmp2.0','tune_interval':step,'add_bias':False,'temperature_distribution':td}
    
    
    ssm1={'sm':'mca','reset_betweenruns':2,'blockupdate':False, 'using_globalbest':True,
          'sample_schedule':sampleschedule,
          'stepmethod':'mxmp2.0','tune_interval':step,'add_bias':False,'temperature_distribution':np.zeros(ni)+2}

    ssm3={'sm':'powell','blockupdate':False}
    
sml=[ssm20,ssm2,ssm0,ssm,ssm0,ssm,ssm0,ssm, ssm1]
#sml=[ssm3,ssm1]
sml2=[ssm3, ssm1]
sml3=[ssm1]
###############Define the final model#########
cvk=2
repeat=1
fold=3
testperc=None
initialmodelpath=''
optm=sml
if True:
    model1={'scorers':ssmodel,'bmtype':bmtype,'searches':ssmodel, 'runsperscorer':ni,
        'dsearches':dsList,'sml':optm,'cvk':cvk,'repeat':repeat,'fold':fold}
    
    if testperc!=None:
        model1['testperc']=testperc
        
    if len(initialmodelpath)>0:
        model['initialmodelpath=']=initialmodelpath

##############run the model###############
rrn=0
#0: test optimizer
#1: find best par--default
#2: find best par -- no elimination
#3: eval all pars
#4: blockbyblock search
inimp=''#copy values from initial model
if len(inimp)>0:
    #blp='/bell3/gqdong/statpot/results/ppd4s_zdock_ppd4/top1000_rlrank__rmsd10ORirmsd4FIRST+top1000_rlrank__rmsd5ORirmsd2FIRST/runs/38772-nan,nan_nan+nan_nan+nan_2.836_0.000/'
    blp=inimp
    bm=pickle.load(open(blp+'cvmodel.pickle'))
    bestpars=bm['allresult'][0]['repna']
    bm['model']['bmtype']['dslist']=['ppd4s']    
    so=scorer(model=bm['model'])
    print so.assess_model()
    
    print so.assess(bm['allresult'][0]['bestpar'])
    #set the default values for nm
    nm=convert2old(model1)
    allsearches=nm['searches']
    
    for nscorer,scorerdict in zip(nm['scorers'],so.model['scorers']):
        if scorerdict['type']=='sf':
            for key in ['par','parvalue','ratio']:
                for i in range(len(scorerdict[key])):
                    nscorer[key][i]=scorerdict[key][i]
    
    for i in range(len(allsearches)):
        for j in range(len(allsearches[i]['object'][allsearches[i]['key']])):
            bms=allsearches[i]
            allsearches[i]['InitialGenerator']=copy.deepcopy(bms['object'][bms['key']])
    #pdb.set_trace()

    so1=scorer(model=nm)
    print so1.assess_model()
    opt=optimizer(scorer=so1)
    print so1.assess_model()

    opt.get_initial_value()    
    print so1.assess(opt.initialvalue)
    
    spso=sps(modellist=[nm])
    spso.cv()#spl.eval_allpars()
    #pdb.set_trace()
    logpath=spso.cvlist[0].logpath
    pdb.set_trace()
     
if rrn==0:
    tm=convert2old(copy.deepcopy(model1))
    so=scorer(model=tm)
    print so.assess_ideal()
    print so.assess_model()
    opt=optimizer(scorer=so)
    opt.get_initial_value()
    pdb.set_trace()
    cProfile.run('opt.optimize()')
    optres=opt.optimize()
    print so.assess(opt.initialvalue)
    pdb.set_trace()
elif rrn==1:
    spl=spss(model=model1)
    spl.find_best_par()
    pdb.set_trace()
elif rrn==2:
    spl=spss(model=model1)
    spl.currentcutoffpercinitial=0.0
    spl.currentcutoffpercratio=0.0
    spl.maximumround=200
    #spl.eval_allpars()
    spl.find_best_par()
    spl.log()
    pdb.set_trace()    
elif rrn==3:
    spl=spss(model=model1)
    spl.eval_allpars()
    pdb.set_trace()
elif rrn==4:
    #blp='/bell3/gqdong/statpot/results/ppd4s_zdock_ppd4/top1000_rlrank__rmsd10ORirmsd4FIRST+top1000_rlrank__rmsd5ORirmsd2FIRST/runs/38772-nan,nan_nan+nan_nan+nan_2.836_0.000/'
    blp=inimp
    bm=pickle.load(open(blp+'cvmodel.pickle'))
    bestpars=bm['allresult'][0]['repna']
    so=scorer(model=bm['model'])
    bm['model']['bmtype']['dslist']=['ppd4s']
    print so.assess(bm['allresult'][0]['bestpar'])
    
    print so.assess_model()
        
    if 1:
        for iii in range(3):
            for block in blocks:
                print block
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
                    except Exception,e:
                        print e
                        pdb.set_trace()
