"""
   SOAP utility module

"""
from env import *



def squeeze_dict(d):
    nd={}
    for key in d:
        if d[key]:
            nd[key]=copy.deepcopy(d[key])
    return nd

def squeeze_list(dl):
    ndl=[]
    for item in dl:
        collect=True
        for item2 in ndl:
            if item==item2:
                collect=False
                break
        if collect:
            ndl.append(item)
    return ndl

def get_array_mem(na):
    nas=na.shape
    ams=8
    for item in nas:
        ams=ams*item
    return ams

def get_residues_atoms():
    fh=open(runenv.basedir+'lib/atmcls-mf.lib')
    fhc=fh.read()
    fh.close()
    resname1=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V']
    resname3=[]
    resatomdict={}
    atoms=[]
    allresidues=re.findall('\\n\s*ATOM\s*\'([A-Z]{3})\'\s*\'([A-Z0-9]+)\'',fhc)
    resatoms=[]
    allatoms=[]
    cr=''
    k=-1
    for item in allresidues:
        if item[0]=='HSD' or item[0]=='CSS':
            continue
        if item[0]==cr:
            resatoms.append(item[1])
            allatoms.append(item[1])
        else:
            k=k+1
            if resatoms:
                atoms.append(resatoms)
                print k
                resatomdict[cr]=resatoms
            resatoms=[item[1]]
            allatoms.append(item[1])
            resname3.append(item[0])
            cr=item[0]
            print item[0]
    atoms.append(resatoms)
    resatomdict[cr]=resatoms
    resname3.append('HSD')
    resname3.append('CSS')
    resatomdict['HSD']=resatomdict['HIS']
    resatomdict['CSS']=resatomdict['CYS']
    return {'resname1':resname1,'resname3':resname3,'resatoms':atoms,'atomset':set(allatoms),'resatomdict':resatomdict}

def get_sidechainatom_index():
    fh=open(runenv.basedir+'lib/atmcls-mf.lib')
    fhc=fh.read()
    fh.close()
    na=np.zeros(158)
    gl=fhc.split('ATMGRP')
    for i in range(158):
        if '\'N\'' in gl[i+1] or '\'C\'' in gl[i+1] or '\'CA\'' in gl[i+1] or '\'O\'' in gl[i+1]:
            na[i]=0 #mainchain
        else:
            na[i]=1 #sidechain
    return na
                    
def get_ph(code):
    bdir='/salilab/park2/database/pdb/divided/'
    fh=gzip.open(bdir+code[1:3]+'/pdb'+code[0:4]+'.ent.gz')
    fc=fh.read()
    fh.close
    try:
        phr=re.search('REMARK 200\s+PH\s+\:\s*([0-9\.]+)',fc)
        ph=float(phr.group(1))
    except:
        return 0
    return ph

def float2int(tbca):
    amax=tbca.max()
    print 'array max '+str(amax)
    if amax<=255:
        na=tbca.astype(np.uint8)
    elif amax<=65535:
        na=tbca.astype(np.uint16)
    elif amax<=4294967295:
        na=tbca.astype(np.uint32)
    else:
        print 'the data is too large to be converted to int to save memory'
        na=tbca
    return na

def metalions():
    return ['K' , 'Y' , 'ZN', 'CD', 'PA', 'GD', 'PR', 'NA', 'MO', 'RB' ,'YB'   , 'MN'   , 'RE'  , 'HO'   , 'AU'   , 'PU'   , 'AL'   , 'BE'   , 'IR'   , 'PT' ,'FE'   , 'CA'   , 'AG'   , 'CM'   , 'HG'   , 'SI'   , 'CO'   , 'SR'   ,'MG'   , 'PB' ,'GA'   , 'SM'   , 'BA'   , 'CU'   , 'NI'   , 'TE'   , 'CR'   , 'TB'   , 'PD'   , 'OS' ,'EU'   , 'W'   , 'B'   , 'TL'   , 'LI'   , 'U'   , 'V'   , 'RU'   , 'AS'   , 'LU' ,'CS' ]    

def isfirstlowest(array1):
    al=len(array1)
    print al
    re=1
    c2=array1[0]
    print c2
    code="""
         int i;
         for(i=0;i<al;i++){
            if (array1(i)<=c2){
                re=0;
                break;
            }
         }
         """
    weave.inline(code,['array1','al','re','c2'],type_converters=converters.blitz, compiler='gcc')
    return re

def filter_lib(atomtype='CA',ofp='atmcls-mf.lib',ffp='caatm.lib'):
    fh=open(ofp)
    fc=fh.read()
    fh.close()
    fcl=fc.split('ATMGRP')
    nfcl=[]
    for item in fcl:
        if item[3:5]=='CA':
            nfcl.append(item)
    fh=open(ffp,'w')
    nls='ATMGRP'+'ATMGRP'.join(nfcl)
    fh.write(nls)
    fh.close()
    
def extendsum(array1,array2):
    ref=array1
    reft=array2
    ref=ref.reshape(list(ref.shape)+list(np.ones(len(reft.shape))))+(reft.reshape(list(np.ones(len(ref.shape)))+list(reft.shape)))
    return ref

def read_bndgrp(libpath=''):
    if not libpath:
        libpath='bndgrp.lib'
    fh=open(libpath)
    fhc=fh.read()
    fh.close()
    bndgrps=re.findall('BOND \'([A-Z0-9]+)\' \'([A-Z0-9]+)\' \'([A-Z0-9]+)\'',fhc)
    return bndgrps

class bndsep(object):
    def __init__(self):
        self.load_bndsepdict()
        
    def load_bndsepdict(self):
        fh=open(runenv.libdir+'bndsep.pickle')
        self.bndsepdict=pickle.load(fh)
    
    def get_bndsep(self,atm1n,atm2n,res1, res2,ss):
        if ss<0:
            ss=-ss
            temp=atm1n
            atm1n=atm2n
            atm2n=temp
            temp=res1
            res1=res2
            res2=temp
        if ss==0:
            return self.bndsepdict[res1][atm1n+atm2n]
        else:
            bs=self.bndsepdict[res1][atm1n+'C']
            bs=bs+(ss-1)*3+self.bndsepdict[res2][atm2n+'N']+1
            return bs
    
class mypickle(object):
    """ this class will save the numpy array seperately when pickling the object to speed up performance and avoid pickle problems with large array, /////old version: only works when the numpy array is a varialbe of an object. WILL NOT WORK if you have a list of numpy arrays.
    """
    def __init__(self):
        self.npyobjlist=[]
        self.npyobjidlist=[]
        self.npyobjnamelist=[]
        self.classidlist=[]
        self.arraysize=0
        self.indexlen=0
        self.indexlist=[]
    
    def dump(self,obj,nameprefix):
        self.arraysize=0
        self.dumpnparray(obj,nameprefix)
        fh=open(nameprefix+'.pickle','wb')
        pickle.dump(obj,fh)
        fh.close()
        print "dump finished"
        self.get_array_size(nameprefix)
        return self.arraysize
        
    def dumpnparray(self,obj,nameprefix):
        #print "checking object "+obj.__repr__
        if isinstance(obj,list) or isinstance(obj,tuple):
            if len(obj)>3000:
                return 0
            self.dumpkeys(obj,range(0,len(obj)),nameprefix)
        elif isinstance(obj,dict):
            if len(obj)>3000:
                return 0
            self.dumpkeys(obj,obj.keys(),nameprefix)          
        elif '__dict__' in dir(obj) and str(type(obj)).startswith('<class'):
            self.dumpkeys(obj.__dict__,obj.__dict__.keys(),nameprefix)
        
    def dumpkeys(self,od,keys,nameprefix):
        if id(od) in self.classidlist:
            return 0
        self.classidlist.append(id(od))
        for key in keys:
            item=od[key]
            if isinstance(item,np.ndarray):
                ns=item.shape
                nl=1
                for nsi in ns:
                    nl=nl*nsi
                if nl<1000:
                    continue
                self.savesinglearray(od,key,nameprefix)
            else:
                self.dumpnparray(item,nameprefix)
                
    def savesinglearray(self,obj,key,nameprefix):
        dtype=obj[key].dtype
        if not obj[key].flags['OWNDATA']:
            obj[key]=None#you should relink the array afterwords
            return 0
        elif id(obj[key]) in self.npyobjidlist:
            cid=self.npyobjidlist.index(id(obj[key]))
        else:
            cid=len(self.npyobjidlist)
            print 'saving numpy array to hdf5'+' '+str(key)+' shape '+str(obj[key].shape)
            tf=h5py.File(nameprefix+str(cid),'w')
            self.arraysize=self.arraysize+obj[key].nbytes
            na=tf.create_dataset('na',data=obj[key],compression='lzf')
            #np.save(nameprefix+str(cid),obj[key])
            tf.close()
            self.npyobjidlist.append(id(obj[key]))
        obj[key]=savednpy(nameprefix+str(cid))
        obj[key].dtype=dtype
                
    def loadnparray(self,obj):
        if isinstance(obj,list) or isinstance(obj,tuple):
            if len(obj)>1000:
                return 0
            self.loadkeys(obj,range(0,len(obj)))
        elif isinstance(obj,dict):
            if len(obj)>1000:
                return 0
            self.loadkeys(obj,obj.keys())
        elif '__dict__' in dir(obj) and str(type(obj)).startswith('<class'):
            self.loadkeys(obj.__dict__,obj.__dict__.keys())

    def loadkeys(self,od,keys):
        if id(od) in self.classidlist:
            return 0
        self.classidlist.append(id(od))
        for key in keys:
            item=od[key]
            if 'savednpy' in str(type(item)):
                self.loadsinglearray(od,key)
            else:
                self.loadnparray(item)
                
    def loadsinglearray(self,obj,key):
        sao=obj[key]
        if sao.name in self.npyobjnamelist:
            ind=self.npyobjnamelist.index(sao.name)
        else:
            self.npyobjnamelist.append(sao.name)
            print "loading numpy array"+sao.name
            if sao.name.endswith('npy'):
                na=np.load(sao.name)
            else:
                tf=h5py.File(sao.name,'r')
                if len(self.indexlist)==0 or not self.indexlen in tf['na'].shape:
                    na=np.empty(tf['na'].shape,dtype=sao.dtype)
                    tf['na'].read_direct(na)
                else:
                    na=self.load_part_indexlist(tf,sao.dtype)
                tf.close()
            self.npyobjlist.append(na)
            ind=len(self.npyobjlist)-1
        obj[key]=self.npyobjlist[ind]
                
    def load_part_indexlist(self,tf,dtype):
        nindexlist=squeeze_indexlist(self.indexlist)
        if tf['na'].shape[0]==self.indexlen:
            nashape=[nindexlist[-1][-1]]+list(tf['na'].shape)[1:]
            na=np.zeros(nashape,dtype=dtype)
            for ind, nind in zip(self.indexlist,nindexlist):
                na[nind[0]:nind[1],...]=tf['na'][ind[0]:ind[1]]
        elif tf['na'].shape[1]==self.indexlen and len(tf['na'].shape)==2:
            nashape=[tf['na'].shape[0],nindexlist[-1][-1]]
            na=np.zeros(nashape,dtype=dtype)
            for ind, nind in zip(self.indexlist,nindexlist):
                na[:,nind[0]:nind[1]]=tf['na'][:,ind[0]:ind[1]]
        else:
            raise Exception('can not load_part_index in mypickle')
        return na
                
    def load(self,nameprefix):
        print "loading pickle file"
        name=nameprefix if nameprefix.endswith('.pickle') else nameprefix+'.pickle'
        fh=open(name,'rb')
        obj=pickle.load(fh)
        fh.close()
        print "finished loading pickle file"
        #fl=os.listdir('./')
        #fl=[item for item in fl if item.startswith(nameprefix) and item.endswith('npy')]
        #if len(fl)>0:
        self.loadnparray(obj)
        print "load finished"
        return obj
    
    def loadpickleonly(self,nameprefix):
        fh=open(nameprefix+'.pickle','rb')
        obj=pickle.load(fh)
        fh.close()
        return obj
    
    def get_array_size(self,nameprefix):
        fl=os.listdir('./')
        fl=[item for item in fl if item.startswith(nameprefix) and item.endswith('pickle')]
        for f in fl:
            self.arraysize=self.arraysize+os.stat(f).st_size
        return self.arraysize

class savednpy(object):
    def __init__(self,name):
        self.name=name
        self.dtype=[]
        #if not self.name.endswith('.npy'):
        #    self.name=self.name+'.npy'
def calc_zscore(sl):
    std=np.std(sl)
    mean=np.mean(sl)
    return (sl-mean)/std
    
def decode_genmethod(genmethod):
    """
    Decode genmethod str
    """
    csmin=0
    csmax=99
    if re.search('cs([0-1]+)',genmethod):
         errsre=re.search('cs([0-1]+)',genmethod)
         csv=int(errsre.group(1))
         if csv==0:
             csmax=0
         elif csv==1:
             csmin=1
    bs=-1
    bsm=-1
    if re.search('bs([0-9]+)',genmethod):
        bs=int(re.search('bs([0-9]+)',genmethod).group(1))
        bsm=99999999
    if re.search('bsm([0-9]+)',genmethod):
        bsm=int(re.search('bsm([0-9]+)',genmethod).group(1))
        if bs==-1:
            bs=0
    kcut=0
    if re.search('ss([0-9]*)',genmethod):
        kcs=re.search('ss([0-9]*)',genmethod)
        kcs=kcs.group(1)
        kcut=int(kcs)
    if re.search('ssm([0-9]*)',genmethod):
        kcsm=re.search('ssm([0-9]*)',genmethod)
        kcsm=kcsm.group(1)
        kcut1=int(kcsm)
    else:
        kcut1=9999
    errsv=0;
    if re.search('err([0-9\.]*)',genmethod):
        errsre=re.search('err([0-9\.]*)',genmethod)
        errs=errsre.group(1)
        errsv=float(errs)
    if re.search('dsp',genmethod):
        ssp=True
    else:
        ssp=False
    return csmin, csmax,kcut, kcut1, bs,bsm, errsv, ssp

def check_remote_file_exist(rs,rd,rfl):
    proc1=subprocess.Popen('ssh '+rs+' '+'\'cd '+rd+'; ls \'',shell=True,stdout=subprocess.PIPE)
    runstats=proc1.communicate()[0]
    fl=runstats.split('\n')
    fl=[item.strip() for item in fl]
    efl=[]
    nefl=[]
    for f in rfl:
        if f.endswith('*'):
            if re.search(f[:-1],runstats):
                efl.append(f)
            else:
                nefl.append(f)
        else:
            if f.strip() in fl:
                efl.append(f)
            else:
                nefl.append(f)
    return efl,nefl
    
def list2array(tl):
    #This only works if all the element of list have the same length except the most inner list.
    ldim=[]
    le=tl
    lastflat=True
    while isinstance(le,list):
        if isinstance(le[0],list):
            ldim.append(len(le))
            ole=le
            le=le[0]
        else:
            ll=0
            for stl in ole:
                ll=ll+len(stl)
                if len(stl)!=len(ole[0]):
                    lastflat=False
            if not lastflat:
                ldim[-1]=ll
            else:
                ldim.append(len(stl))
            break
    npa=np.zeros(ldim)
    list2array_recursive(lastflat,[],tl,npa)
    return npa.squeeze()
         
def list2array_recursive(lastflat,cil,cl0,npa):
    if not isinstance(cl0[0],list):
        for i in range(0,len(cl0)):
            try:
                npa[tuple(cil+[i])]=cl0[i]
            except:
                pdb.set_trace()
        return True
    elif (not isinstance(cl0[0][0], list)) and (not lastflat):
        li=0
        for cl1 in cl0:
            for cl2 in cl1:
                npa[tuple(cil+[li])]=cl2
                li=li+1
        return True
    else:
        for i in range(0,len(cl0)):
            ncil=copy.deepcopy(cil)
            ncil.append(i)
            list2array_recursive(lastflat,ncil,cl0[i],npa)
             

def make_tmpdir():
    rp='/scratch/gqd'+str(time.time())
    if os.path.isdir(rp) or os.path.isfile(rp):
        time.sleep(0.1)
        rp='/scratch/gqd'+str(time.time())
    os.mkdir(rp)
    return rp+'/'
        
def get_scaledhsv_colormap_inverted(cscale=0.1):
    cscale=0.1/0.84127
    cdict= {'blue': ((0.0, 0.0, 0.0),
                    (0.333333*cscale, 0.0, 0.0),
                    (0.349206*cscale, 0.0625, 0.0625),
                    (0.507937*cscale, 1.0, 1.0),
                    (0.84127*cscale, 1.0, 0.5),
                    (1.0, 0.5, 0.5)),
            'green': ((0.0, 0.0, 0.0),
                    (0.15873*cscale, 0.9375, 0.9375),
                    (0.174603*cscale, 1.0, 1.0),
                    (0.507937*cscale, 1.0, 1.0),
                    (0.666667*cscale, 0.0625, 0.0625),
                    (0.68254*cscale, 0.0, 0.0),
                    (0.84127*cscale, 0.0, 0.5),                        
                    (1.0, 0.5, 0.5)),
            'red': ((0.0, 1.0, 1.0),
                    (0.15873*cscale, 1.0, 1.0),
                    (0.174603*cscale, 0.96875, 0.96875),
                    (0.333333*cscale, 0.03125, 0.03125),
                    (0.349206*cscale, 0.0, 0.0),
                    (0.666667*cscale, 0.0, 0.0),
                    (0.68254*cscale, 0.03125, 0.03125),
                    (0.84127*cscale, 0.96875, 0.5),
                    (1.0, 0.5, 0.5))}
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap
        
def get_scaledhsv_colormap(cscale=0.1):
    s=cscale/(1-0.31745999999999996)
    b=(1-cscale-0.31745999999999996)/(1-0.31745999999999996)
    cdict={'blue': [(0.0, 0.5, 0.5),
                    (s*0.31745999999999996+b, 0.5, 1.0),
                    (s*0.49206300000000003+b, 1.0, 1.0),
                    (s*0.650794+b, 0.0625, 0.0625),
                    (s*0.666667+b, 0.0, 0.0),
                    (1.0, 0.0, 0.0)],
            'green': [(0.0, 0.5, 0.5),
                    (s*0.31745999999999996+b, 0.5, 0.0),
                    (s*0.333333+b, 0.0625, 0.0625),
                    (s*0.49206300000000003+b, 1.0, 1.0),
                    (s*0.8253969999999999+b, 1.0, 1.0),
                    (s*0.84127+b, 0.9375, 0.9375),
                    (1.0, 0.0, 0.0)],
            'red': [(0.0, 0.5, 0.5),
                    (s*0.31745999999999996+b, 0.5, 0.03125),
                    (s*0.333333+b, 0.0, 0.0),
                    (s*0.650794+b, 0.0, 0.0),
                    (s*0.666667+b, 0.03125, 0.03125),
                    (s*0.8253969999999999+b, 0.96875, 0.96875),
                    (s*0.84127+b, 1.0, 1.0),
                    (1.0, 1.0, 1.0)]}
    my_cmap = matplotlib.colors.LinearSegmentedColormap('my_colormap',cdict,256)
    return my_cmap        

def convert2ndarray(values):
    if isinstance(values, np.ndarray):
        return values
    elif isinstance(values, list) or isinstance(values, tuple):
        return np.array(values)
    elif isinstance(values, (int, long, float)):
        return np.array([values])
    else:
        try:
            return np.array(values)
        except:
            raise Exception('Don\'t know how to convert '+str(values)+' to numpy array')
            
def dummy(*args):
    pass

def any2str(anything):
    if isinstance(anything, (np.ndarray,float,int,str)):
        return str(anything)
    elif anything==None:
        return 'None'
    elif isinstance(anything, (list,tuple)):
        nl=[]
        for item in anything:
            nl.append(any2str(item))
        return '['+', '.join(nl)+']'
    elif isinstance(anything,dict):
        keys = anything.keys()
        keys.sort()
        dl=[any2str(anything[key]) for key in keys]
        dls=any2str(dl)
        nds=[]
        for key,value in zip(keys,dl):
            #if key=='searches':
            #    pdb.set_trace()
            nds.append(str(key)+':'+str(value))
        return '{'+', '.join(nds)+'}'
    else:
        pdb.set_trace()
        raise Exception('Donot how to convert anything to string: '+str(anything))
        
def build_scorer(scorerpath='./'):
    os.chdir(scorerpath)
    fh=open(scorerpath+'bestmodel.pickle')
    model=pickle.load(fh)
    fh.close()
    for scorer in model['scorers']:
        if scorer['type']=='scaledsp':
            ssp=scaledsp(model=scorer)
            print os.system('ln -s '+ssp.ssppath+'.hdf5 '+ssp.sspname+'.hdf5')
    print os.system('cp ~/Dropbox/Code/calc_score.py ./')
    


    
def sameobjectIn2list(list1,list2):
    for item1 in list1:
        for item2 in list2:
            if item1 is item2:
                return True
    return False
    
def findsameobjectkeys(d,list2):
    # find the keys in d, whose values are the same objects in list2
    kl=[]
    il=[]
    for key in d:
        for i in range(len(list2)):
            item2=list2[i]
            if d[key] is item2:
                #print "same key " +str(key)+' value '+str(d[key])
                kl.append(key)
                il.append(i)
    return [kl,il]
    
def convert2float(l):
    for i in range(len(l)):
        l[i]=float(l[i])
        
def checknodememory(taskid):
    fl=os.listdir('./')
    fl=[f for f in fl if f.startswith('s.'+str(taskid))]
    runsonsamenode
    
def get_system_freemem():
    process = subprocess.Popen("free -m | grep buffers/cache",
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        )
    result = process.communicate()[0].split()[3]
    return float(result)
    
def get_starting_jobs():
    fl=os.listdir(runenv.serverUserPath+'runstatus/')
    fl=[f for f in fl if f.startswith(hcn+'#')]
    return fl

def squeeze_indexlist(indexlist):
    nindexlist=[]
    cp=0
    for item in indexlist:
        nindexlist.append([cp,cp+item[1]-item[0]])
        cp=nindexlist[-1][-1]
    return nindexlist

def load_scorer(originalscorer,indexlen):
    """loading scorer, used for tranfering scoring across computers"""
    print "loading scorer"
    gc.collect()
    mp=mypickle()
    scorer=copy.deepcopy(originalscorer)
    mp.indexlen=indexlen
    mp.indexlist=scorer.dsss.ds.indexlist
    mp.loadnparray(scorer) #load only revelent part of the array
    scorer.relink_scorearray()
    scorer.dsss.ds.indexlist=squeeze_indexlist(originalscorer.dsss.ds.indexlist)
    #intialize the scorer for parallel job if specified. Get scorer list.
    scorer.singlescorerlist=scorer.get_ds_list()
    #pdb.set_trace()
    try:
        if scorer.loadfilter!=None:
            scorer.filter(scorer.loadfilter)
        scorer.loadfilter=None
    except:
        pass
    return scorer

class Mymodel(model):
    """
    Modeller model class with customized functions. 
    """
    def select_loop_atoms(self,loops):
        s=selection()
        if not isinstance(loops[0],list):
            self.add_loop2selection(loops,s)
        else:
            for loop in loops:
                self.add_loop2selection(loop,s)
        return s
   
    def add_loop2selection(self,loop,s):
        if loop[0]=='':
            loop[0]=self.chains[0].name
        try:
            s.add(selection(self.residue_range(str(loop[2])+':'+loop[0],str(loop[3])+':'+loop[0])))
        except:
            lind=self.residues[str(loop[2])+':'+loop[0]].index
            s.add(selection(self.residues[lind:(lind+loop[1])]))    
    
    def load_model_file(self,nativef):
        if nativef[-2:]=='gz':
            fh=gzip.open(nativef)
        else:
            fh=open(nativef)
        fhc=fh.read()
        self.fname=nativef
        self.atomlines=re.findall('\\n(ATOM.*)','\n'+fhc+'\n')
        self.hetatomlines=re.findall('\\n(HETATM.*)','\n'+fhc+'\n')
        self.gen_atomlist()
        
    def gen_atomlist(self):#[chainname,residue number,residue name,atom name, alternative indicator, x,y,x,occupancy, B-factor, elem symbol, charge,index]
        self.atomlist=[]
        k=-1
        for atomline in self.atomlines:
            k=k+1
            atompars=[atomline[21],int(atomline[22:26]),atomline[17:20].strip(),atomline[12:16].strip(),atomline[16], float(atomline[30:38]),float(atomline[38:46]),float(atomline[46:54]),float(atomline[54:60]),float(atomline[60:66]),k]
            #if len(atomline)>=80:
            #        atompars.append(atomline[76:78].strip())
            #        atompars.append((atomline[78:80]))
            self.atomlist.append(atompars)
        self.atomlist.sort(key=itemgetter(0,1)) #sort the atom file by the chainids and the residue numbers, because some of the pdb files is messed up.
        newatomlines=[]
        for atompars in self.atomlist:
            newatomlines.append(self.atomlines[atompars[-1]])
        self.atomlines=newatomlines
    
    def gen_chainlist(self): #get the chain name and pos [chain name, chainstartpos,chainendpos, [[residuename,residue num, residue startpos, residue endpos]...]]
        currentchain=''
        chainlist=[]
        currentresidue=''
        residuelist=[]
        try:
            k=-1
            for atompars in self.atomlist:
                k=k+1
                if atompars[0]!=currentchain:
                    currentchain=atompars[0]
                    chainlist.append([atompars[0],k])
                    if len(chainlist)>1:
                        chainlist[-2].append(k-1)
                        chainlist[-2].append(residuelist)
                        chainlist[-2][3][-1].append(k-1)
                        residuelist=[]
                if atompars[1]!=currentresidue:
                    currentresidue=atompars[1]
                    residuelist.append([atompars[2],atompars[1],k])
                    if len(residuelist)>1:
                        residuelist[-2].append(k-1)
            chainlist[-1].append(k)
            chainlist[-1].append(residuelist)
            chainlist[-1][3][-1].append(k)
        except Exception,e:
            traceback.print_exc()
            pdb.set_trace()
        self.chainlist=chainlist
                
    def gen_atomindexdict(self):
        atomindexdict={}
        for k in range(0,len(self.atomlist)):
            atompars=self.atomlist[k]
            atomid=atompars[0]+str(atompars[1])+atompars[2]+atompars[3]
            atomindexdict[atompars[0]+str(atompars[1])+atompars[2]+atompars[3]]=k
        self.atomindexdict=atomindexdict
            
    def only_std_residue_heavy_atom(self):
        resatoms=get_residues_atoms()
        atomlines=[]
        for k in range(0,len(self.atomlist)):
            atompars=self.atomlist[k]
            if atompars[2] in set(resatoms['resname3']) and atompars[3] in set(resatoms['resatomdict'][atompars[2]]):
                atomlines.append(self.atomlines[k])
        self.atomlines=atomlines
        self.gen_atomlist()
        self.gen_chainlist()
        self.gen_atomindexdict()
    
    def transfer_model_atoms(self,mdl2):
        self.gen_atomindexdict()
        mdl2.gen_atomindexdict()
        #copy mdl2 atoms into model 1
        for key in mdl2.atomindexdict:
            try:
                self.atomlines[self.atomindexdict[key]]=self.atomlines[self.atomindexdict[key]][0:11]+mdl2.atomlines[mdl2.atomindexdict[key]][11:]
            except:
                pdb.set_trace()
        
    def transfer_residues(self, mdl2,chain, residuerange):
        #copy coordinates in the residuerange from mdl2 to current model
        selfstartpos=self.get_residue_index(chain,residuerange[0])[2]
        selfendpos=self.get_residue_index(chain,residuerange[1])[3]
        mdl2startpos=mdl2.get_residue_index(chain,residuerange[0])[2]
        mdl2endpos=mdl2.get_residue_index(chain,residuerange[1])[3]
        self.atomlines[selfstartpos:selfendpos+1]=mdl2.atomlines[mdl2startpos:mdl2endpos+1]
        self.gen_atomlist()
        
    def save_model(self,filename):
        fh=open(filename,'w')
        fh.write('\n'.join(self.atomlines))
        fh.close()
    
    def get_residue_index(self,chainid,residuenum):
        for item in self.chainlist:
            if item[0]!=chainid:
                continue
            for res in item[3]:
                if residuenum==res[1]:
                    return res

    def copy_native(self,code,tdir):
        code=code.lower()
        print os.system('cp '+runenv.opdbdir+code[1:3]+'/pdb'+code+'.ent.gz '+tdir)
        print os.system('gunzip '+tdir+'pdb'+code+'.ent.gz ')

    def clash_contact(self,code):
        if os.path.isdir('clashtemp'):
            print os.system('rm -r clashtemp')
        os.mkdir('clashtemp')
        os.chdir('clashtemp')
        self.copy_native(code,'./')
        print os.system('cp ~/Dropbox/Code/crystal_contacts.py ./')
        print os.system('chimera --nogui crystal_contacts.py')
        fh=open('selfclash','r')
        scs=fh.read()
        fh.close()
        fh=open('otherclash','r')
        ocs=fh.read()
        fh.close()
        fh=open('othercontact','r')
        octs=fh.read()
        fh.close()
        clashlist=[]
        contactlist=[]
        rep=re.compile('#([0-9]+) ([A-Z]{3}) ([0-9]+)\.([A-Z]{1}) ([A-Z0-9\.\']+)\s* #([0-9]+) ([A-Z]{3}) ([0-9]+)\.([A-Z]{1}) ([A-Z0-9\.\']+)\s* ([0-9\-\.]+)\s* ([0-9\-\.]+)')
        allcontacts=rep.findall(octs)
        allclashes=rep.findall(scs)+rep.findall(ocs)
        clashlist=[]
        for item in allclashes:
            clashlist.append([[item[0],item[3],item[2],item[1],item[4]],[item[5],item[8],item[7],item[6],item[9]],item[10],item[11]])
        contactlist=[]
        for item in allcontacts:
            contactlist.append([[item[0],item[3],item[2],item[1],item[4]],[item[5],item[8],item[7],item[6],item[9]],item[10],item[11]]) 
        os.chdir('..')
        print os.system('rm -r clashtemp')
        return clashlist,contactlist

    def select(self,loops):
        s=selection()
        for loop in loops:
            try:
                s.add(selection(self.residue_range(str(loop[2])+':'+loop[0],str(loop[3])+':'+loop[0])))
            except:
                lind=self.residues[str(loop[2])+':'+loop[0]].index
                s.add(selection(self.residues[lind:(lind+loop[1])]))
        return s
    
    def get_residues(self,loops):
        return [loop[0]+str(lnum) for loop in loops for lnum in range(loop[2],loop[3]+1)]
    
    def has_nonstd(self,loops):
        s=self.select(loops)
        l=self.get_residues(loops)
        return len(l), len(set(a.residue for a in s))

    def has_clash(self,loops,code):
        clash,_=self.clash_contact(code)
        loopres=self.get_residues(loops)
        clashlist=[]
        for item in clash:
            if item[0][0]=='0':
                clashlist.append(item[0][1]+item[0][2])
            if item[1][0]=='0':
                clashlist.append(item[1][1]+item[1][2])
        rl=[]
        for item in clashlist:
            if item in loopres:
                print 'clash with '+','.join(item)+'  '
                rl.append('clash with '+','.join(item)+'  ')
        return rl
    
     

    def gen_base_file(self,chain,residuerange,basefilename,othername):
        selfstartpos=self.get_residue_index(chain,residuerange[0])[2]
        selfendpos=self.get_residue_index(chain,residuerange[1])[3]
        otheratomlines=self.atomlines[selfstartpos:selfendpos+1]
        self.atomlines[selfstartpos]='replacethis'
        for i in range(selfendpos,selfstartpos,-1):
            self.atomlines.pop(i)
        self.save_model(basefilename)
        self.atomlines=otheratomlines
        self.save_model(othername)
        self.gen_atomlist()


def load_pickle(fn):
    pass




class MemoryMonitor(object):
    """Monitor the memory usage of the current program"""
    def __init__(self):
        """Create new MemoryMonitor instance."""
        self.mul=[]

    def usage(self, name ):
        """Return int containing memory used by user's processes."""
        self.process = subprocess.Popen("ps -aux | grep "+name,
                                        shell=True,
                                        stdout=subprocess.PIPE,
                                        )
        
        self.mm = self.process.communicate()[0].split('\n')[0].split()[5]
        print self.mm
        return int(self.mm)
        
    def monitor(self,name,fn='mlog'):
        fh=open(fn,'a+')
        mu=11
        while mu>10:
            mu=self.usage(name)
            self.mul.append(mu)
            fh.write(str(mu)+'\n')
            time.sleep(5)
            fh.flush()
        fh.close()
        
def load_log(path):
    from collections import defaultdict
    from operator import itemgetter
    dl=os.listdir(path)
    logs=[]
    ld=defaultdict(list)
    for d in dl:
        logpath=os.path.join(path,d,'log')
        if os.path.isfile(logpath):
            ll=open(logpath).read().split('\n')
            for l in ll:
                if re.match('-[1-9]+.*',l) and (l.split()[0],l.split('0.000  ')[-1]) not in logs:
                    logs.append((l.split()[0],l.split('0.000  ')[-1]))
                    
                    #ld[l.split()[0]].append(l)
    logs.sort(key=itemgetter(0),reverse=True)
    #for key in ld:
    #    print key
    #    for v in ld[key]:
    #        print v
    for l in logs:
        print '    '.join(l)
    return None