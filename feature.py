"""
   SOAP feature module for defining MDT features.

"""

import mdt
import mdt.features
from env import *
         
             
class feature(object):
    """
    Define the features for calculating statistics from structures
    
    :param str features: the string representation of features.
    
    Feature definition mini language::
      
        features => feature [,features]
        feature  => type, number of bins, #, range, [,#, start value] #start value will be zero if not defined
    
    Example::
        
        features='d30#15'# 0-30A, 0.5 bin size    
            
    
    Types are from MDT:
    
    .. literalinclude:: ../feature.py 
       :pyobject: feature._gen_single_feature
        

        
    """
    def __init__(self,features,mlib=None):#mlib does not need to be supplied anymore
        self.features=features
        #self.preprocess_features()
        self.featurelist=[]
        self.featurenames=[]
        self.mlib=mlib
        self.atomclasslib='atmcls-mf.lib'
        self.fclist=re.findall('[a-z]{1,4}[0-9\#\-\.]{1,20}',self.features) #feature  code list
        self.ffclist=copy.deepcopy(self.fclist)
        self.frlist=[] # feature range list
        self.fdlist=[] # feature bin number list
        self.fslist=[] # feature start post list
        for i in range(len(self.fclist)):
            sfc=self.fclist[i]
            if sfc.count('#')==2:
                rer=re.search('([a-z]{1,5})([0-9]{1,8})\#([0-9\-\.]{1,5})\#([0-9\-\.]{1,5})',sfc)
                self.fclist[i]=rer.group(1)
                self.fdlist.append(int(rer.group(2)))
                self.frlist.append(float(rer.group(3)))
                self.fslist.append(float(rer.group(4)))
            if sfc.count('#')==1:
                rer=re.search('([a-z]{1,5})([0-9]{1,8})\#([0-9\.]{1,5})',sfc)
                self.fclist[i]=rer.group(1)
                self.fdlist.append(int(rer.group(2)))
                self.frlist.append(float(rer.group(3)))
                self.fslist.append(0)
            elif sfc.count('#')==0:
                self.frlist.append(0)
                rer=re.search('([a-z]{1,5})([0-9]{1,8})',sfc)
                self.fdlist.append(int(rer.group(2)))
                self.fslist.append(0)
        if not mlib:
            env = runenv.env
            #log.minimal()
            self.mlib=mdt.Library(env)
      
    def _genfeaturelist(self):
        self.read_lib()
        for f,r,d,s in zip(self.fclist,self.frlist,self.fdlist,self.fslist):
            self.featurelist.append(self._gen_single_feature(f,r,d,s))
        return self.featurelist
        
    def _gen_single_feature(self,sfc,end,numofbin,start):
        mlib=self.mlib
        if end > 0:
            ub=mdt.uniform_bins(numofbin, start, (float(end)-float(start))/numofbin)
        if sfc[0:2]=='as':
            return mdt.features.AtomType(mlib, pos2=True)
        elif sfc=='b2':
            return mdt.features.FractionalAtomAccessibility(mlib,bins=mdt.uniform_bins(1, 0,0.1)+mdt.uniform_bins(1, 0.1,0.9))
        elif sfc=='bs2':
            return mdt.features.FractionalAtomAccessibility(mlib,pos2=True, bins=mdt.uniform_bins(1, 0,0.1)+mdt.uniform_bins(1, 0.1,0.9))
        elif sfc[0:1]=='a':
            return mdt.features.AtomType(mlib)        
        elif sfc[0:2]=='rs':
            return mdt.features.ResidueType(mlib,delta=numofbin, pos2=True)
        elif sfc[0:1]=='r':
            return mdt.features.ResidueType(mlib,delta=numofbin)
        elif sfc in ['dr']:
            return mdt.features.ResidueDistance(mlib, bins=ub)
        elif sfc[0:1]=='m':
            return mdt.features.MainchainConformation(mlib)  
        elif sfc[0:2]=='ms':
            return mdt.features.MainchainConformation(mlib,pos2=True)  
        elif sfc[0:2]=='ca':
            return mdt.features.AngleType(mlib)
        elif sfc[0:2]=='cb':
            return mdt.features.BondType(mlib)
        elif sfc[0:2]=='cd':
            return mdt.features.DihedralType(lib)
        elif sfc[0:2]=='co':
            return mdt.features.BondLength(mlib,bins=ub)
        elif sfc[0:2]=='cn':
            return mdt.features.Angle(mlib,bins=ub)
        elif sfc[0:2]=='ci':
            return mdt.features.Dihedral(mlib,bins=ub)
        elif sfc[0:2]=='ps':
            return mdt.features.PsiDihedral(mlib,bins=ub)
        elif sfc[0:2]=='ph':
            return mdt.features.PhiDihedral(mlib,bins=ub)        
        elif sfc in ['d','dca','dmc','dsc','dl']:
            return mdt.features.AtomDistance(mlib,bins=ub)
        elif sfc=='l':
            return mdt.features.ResidueAccessibility(mlib,bins=ub)
        elif sfc=='b':
            return mdt.features.FractionalAtomAccessibility(mlib,bins=ub)
        elif sfc=='bs':
            return mdt.features.FractionalAtomAccessibility(mlib,bins=ub,pos2=True)               
        elif sfc=='ba30':
            return mdt.features.AtomAccessibility(mlib,bins=mdt.uniform_bins(30, 0, 0.5))
        elif sfc=='bas30':
            return mdt.features.AtomAccessibility(mlib,bins=mdt.uniform_bins(30, 0, 0.5),pos2=True)    
        elif sfc=='s2':
            return mdt.features.ResidueIndexDifference(mlib, bins=mdt.uniform_bins(1, -1,1)+mdt.uniform_bins(1, 1,1), absolute=False)
        elif sfc=='s33':
            return mdt.features.ResidueIndexDifference(mlib, bins=mdt.uniform_bins(1, -10015,10000)+mdt.uniform_bins(31, -15,1)+mdt.uniform_bins(1, 16,100000), absolute=False)
        elif sfc in ['dt','dtl']:
            return mdt.features.TupleDistance(mlib,bins=ub)
        elif sfc=='g':
            return mdt.features.TupleAngle1(mlib,bins=ub)
        elif sfc=='gs':
            return mdt.features.TupleAngle2(mlib,bins=ub)
        elif sfc=='h':
            return mdt.features.TupleDihedral1(mlib,bins=ub)
        elif sfc[0]=='t':            
            return mdt.features.TupleType(mlib)
        elif sfc[:2]=='ts':
            return mdt.features.TupleType(mlib, pos2=True)
        else:
            raise Exception('can not find feature '+sfc+' in module feature._gen_single_feature()')
        
        #remember to update feature type and the getlib feature when adding features

    def define_feature_type(self,sfc):
        if sfc[0] in ['a','t','r'] or sfc[:2] in ['ca','cb','cd']:#bug need fix
            return 0 #independent feature
        elif sfc[0] in ['d','b','s','g','h','l','m','c','p']:
            return 1
        else:
            raise Exception('The feature type of '+sfc+' is undefined in feature.define_feature_type()')

    def get_lib(self):
        if self.features in ['posescore']:
            return 'atmcls-plas.lib'
        if re.search('dl|al',self.features):
            return 'atmcls-plas.lib'
        elif re.search('tl|dtl',self.features):
            return 'dicls-lig.lib'
        elif re.search('dca',self.features):
            return 'caatm.lib'
        elif re.search('dmc',self.features):
            return 'atmcls-mc.lib'
        elif re.search('dsc',self.features):
            return 'atmcls-sc.lib'        
        elif re.search('cn',self.features):
            return 'anggrp.lib'
        elif re.search('co',self.features):
            return 'bndgrp.lib'
        elif re.search('ci',self.features):
            return 'impgrp.lib'  
        elif self.features[0] in ['h','g','t'] or self.features[0:2] in ['gs','dt']:
            return 'dicls-all.lib'
        else:
            return 'atmcls-mf.lib'

    def issubset(self, f2):
        shrinkindex=[[] for item in f2.fclist]
        permuteindex=[]
        nclist=[c2 for c2 in f2.fclist if c2 in self.fclist]
        for f,r,d,s,i in zip(self.fclist,self.frlist,self.fdlist,self.fslist,range(len(self.fclist))):
            if not f in f2.fclist:
                return []
            fi=f2.fclist.index(f)
            si=self.issubset_singlefeature(r,d,s,f2.frlist[fi],f2.fdlist[fi],f2.fslist[fi])
            if len(si)==0:
                return []
            else:
                permuteindex.append(nclist.index(f))
                shrinkindex[fi]=si
        return [shrinkindex,permuteindex]
        
    def issubset_singlefeature(self,r,d,s,r2,d2,s2):
        if r==0:
            return [0]
        bs1=(float(r)-s)/d
        bs2=(float(r2)-s2)/d2
        br=np.rint(bs1/bs2)
        if np.allclose(float(br),bs1/bs2) and s2<=s and r2>=r:
            sb=(float(s)-s2)/bs2
            eb=(float(r)-s2)/bs2
            sbn=np.rint(sb)
            ebn=np.rint(eb)
            if np.allclose(sb,sbn) and np.allclose(eb,ebn):
                return range(int(sbn),int(ebn)+1,int(br))
            else:
                return []
        else:
            return []
   
    def get_featurelist(self):
        if self.featurelist:
            return self.featurelist
        featurelist=self._genfeaturelist()
        self.featurelist=featurelist
        return featurelist

    def get_dimension(self):
        self.get_featurelist()
        return len(self.featurelist)

    def get_all_bins_centers(self):
        featurelist=self.get_featurelist()
        allbins=[]
        #pdb.set_trace()
        for pos in range(0,self.get_dimension()):
            allbins.append(self.get_bin_centers(pos))
        return allbins
    
    def get_featurenames(self):
        if self.featurenames:
            return self.featurenames
        for f in self.get_featurelist():
            self.featurenames.append(str(f).split(' ')[0].split('.')[-1])
        return self.featurenames
    
    def get_bins(self,pos=0):
        featurelist=self.get_featurelist()
        #if isinstance(featurelist,(list,tuple)):
        #    sf=featurelist[pos]
        #else:
        #    sf=featurelist
        if self.frlist[pos]>0:
            bins=[]
            binwidth=(self.frlist[pos]-float(self.fslist[pos]))/self.fdlist[pos]
            for i in range(self.fdlist[pos]):
                bins.append([self.fslist[pos]+binwidth*i,self.fslist[pos]+binwidth*(i+1)])
            return bins
        else:
            return []

    def get_runmem(self):
        s=8
        for item in self.fdlist:
            s=s*(item)
        runmem=float(s)/1000000000+0.6
        if runmem<0.5:
            runmem=0.5
        runmem=int(float(runmem)*1000)/1000.0
        if 'l' in self.fclist:
            runmem=1.5
        return runmem
      
    def get_bin_centers(self,pos=0):
        bins=self.get_bins(pos)
        rp=[]
        if len(bins)>0:
            for i in range(len(bins)):
                rp.append((bins[i][0]+bins[i][1])/2)
            return rp
        else:
            return []

    def read_lib(self,mlib=[]):
        if mlib:
            mlibr=mlib
        else:
            mlibr=self.mlib
        libfile=self.get_lib()
        if re.search('dicls',libfile):
            mlibr.tuple_classes.read(runenv.libdir+libfile)
        else:
            mlibr.atom_classes.read(runenv.libdir+libfile) 


    def issymetry(self):
        fn=self.get_featurenames()
        tl=[]
        dl=[]
        for i in range(0,len(self.fdlist)):
            if fn[i] in ['AtomType','ResidueType']:
                tl.append(fn[i])
                dl.append(self.fdlist[i])
        if len(tl)>2:
            raise Exception('feature.issymetry() don\'t know how to handle features with three typesss')
        elif len(tl)<2:
            return False
        else:
            if tl[0]==tl[1] and dl[0]==dl[1]:
                return True
            else:
                print "attention: feature "+self.features+" is not symmetric"
                return False
        
    def get_featuretypepos(self):
        self.ifp=[]
        self.dfp=[]
        rl2=[]
        for i in range(0,len(self.fclist)):
            if self.define_feature_type(self.fclist[i])==1:
                self.dfp.append(i)
            elif self.define_feature_type(self.fclist[i])==0:
                self.ifp.append(i)
        return self.ifp, self.dfp
     
    def get_independent_feature(self):
        fs=''
        for i in self.ifp:
            fs=fs+self.ffclist[i]
        return fs
     
    def get_feature_names(self):
        if len(self.fclist)>1:
            raise Exception('we can only get feature names for single feature')
        fl=self.get_featurelist()
        m=mdt.Table(self.mlib,features=fl)
        sl=[]
        for bin in m.features[0].bins:
            sl.append(bin.symbol)
        return sl
       
class atmsprop(object):
    def __init__(self,libfile):
        self.vdrd={'C':1.8,'N':1.64,'O':1.46,'S':1.77,'P':1.88,'F':1.56,'Cl':1.74,'Br':1.98,'I':2.09}
        self.libfile=libfile
        self.atomNameList=self.read_libfile()
        self.atomTypeList=[self.get_atomtype(an) for an in self.atomNameList]
        self.atomrl=[self.vdrd[at] for at in self.atomTypeList]
        
    def get_atomtype(self,atomname):
        if 'Cl' in atomname:
            return 'Cl'
        for key in self.vdrd:
            if key in atomname:
                return key
        raise exception('atom type not known')
    
    def read_libfile(self):
        fh=open(runenv.basedir+'lib/'+self.libfile)
        fhc=fh.read()
        fh.close()
        if 'DBLGRP' in fhc:
            agl=[l[8:-1] for l in fhc.split('\n') if l.startswith('DBLGRP')]
        elif 'ATMGRP' in fhc:
            agl=[l[8:-1] for l in fhc.split('\n') if l.startswith('ATMGRP')]
        return [a for a in agl if len(a)>0]
    
    def get_repulsion(self,ind,bins):
        atm1ind=ind/len(self.atomTypeList)
        atm2ind=ind%len(self.atomTypeList)
        radius=self.atomrl[atm1ind]+self.atomrl[atm2ind]
        pv=np.zeros(len(bins))
        for i in range(len(bins)):
            if bins[i]>=radius:
                break
            pv[i]=np.exp(-0.5*((bins[i]-radius)/0.05)**2)/(0.05*2.51)
        return pv
