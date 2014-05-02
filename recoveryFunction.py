"""
   SOAP recovery function module

"""

from statsTable import *
import numpy as np
import scipy.interpolate

class sf(object):
    def __init__(self,data=None,sftype='',par=[],parvalue=[],features='',type='',model=[],rawsp=None,*args,**args2):
        if model:
            sftype=model['sftype']
            features=model['features']
            par=model['par']
            parvalue=model['parvalue']
        self.data=data
        self.type=sftype
        self.par=par
        self.parvalue=parvalue
        self.features=features
        if features:
            self.initialize()
        if rawsp!=None:
            self.svdu,v=scaledsp(model=rawsp).get_svd()
            del v

    def write_potential(self,affix,genmethod, ratio=1):
        ro=rawsp(pdbset='X_2.2A',features=self.features,genmethod=genmethod)
        ro.write_hdf5(self.features+'.'+affix+'.hdf5', ratio*self.get_sf(), permute=False)

    def initialize(self):
        if self.features.endswith('a158as158'):
            features=self.features[:-9]
        else:
            features=self.features
        self.bins=feature(features).get_all_bins_centers()
        self.rs=[]
        self.rsi=[]
        self.types=self.type.split('_')
        self.sfdim=np.zeros(len(self.types))
        self.ntypes=[]
        self.stypes=[]
        self.ntpos=[]
        self.stpos=[]
        self.nlen=1
        self.slen=1
        self.spos=[0]
        self.stypepardim=[]
        if self.type.startswith('sippl'):
            return 0
        for i in range(0,len(self.types)):
            if re.match('\d+',self.types[i]):
                self.types[i]=int(self.types[i])
                self.ntypes.append(self.types[i])
                self.nlen=self.nlen*self.types[i]
                self.ntpos.append(i)
            else:
                self.stypes.append(self.types[i])
                type=self.types[i]
                if type in ['ig','dfire']:
                    self.stypepardim.append(1)#length of sf
                elif type.startswith('spline'):
                    if len(type)>6:
                        self.stypepardim.append(int(type[6:]))
                    else:
                        self.stypepardim.append(len(self.par))
                elif type.startswith('pspline'):
                    if len(type)>7:
                        self.stypepardim.append(int(type[7:]))
                    else:
                        self.stypepardim.append(len(self.par))
                elif type in ['bins','logbins']:
                    self.stypepardim.append(len(self.parvalue))
                else:
                    self.stypepardim.append(len(self.parvalue))
                self.slen=self.slen*self.stypepardim[-1]
                self.spos.append(self.spos[-1]+self.stypepardim[-1])
                self.stpos.append(i)
        for i in range(0,len(self.bins)):
            bin=self.bins[i]
            self.rs.append(len(bin))
            if i in self.stpos:
                self.rsi.append(range(0,len(bin)))
            else:
                self.rsi.append(0)
        self.sfv=np.zeros(self.rs)
        if self.nlen*self.slen!=len(self.parvalue):
            print 'The length of parvalue is not matched'
            pdb.set_trace()
            raise Exception('The length of parvalue is not the same as needed for building the reference distribution')
                 
    def get_sf(self, returnreft=False):
        if self.type.startswith('sippl'):
            return np.log(self.data.sum(0)+0.000000000000000001)
        for i in range(0,self.nlen):
            if len(self.ntypes)==0:
                ai=i
            else:
                ai=np.unravel_index(i,self.ntypes)
            if len(self.par)>0:
                par=self.par[i*self.slen:(i+1)*self.slen]
            else:
                par=[]
            parvalue=self.parvalue[i*self.slen:(i+1)*self.slen]
            reft=self.getmd(par,parvalue)
            if returnreft or self.rsi==[[]]:
                self.sfv=reft
                return reft
            rsi=copy.copy(self.rsi)
            for k in range(0,len(self.ntpos)):
                rsi[self.ntpos[k]]=ai[k]
            self.sfv[rsi]=reft
        if np.isnan(self.sfv).sum()>0:
            #pdb.set_trace()
            raise NanInScore('nan in sf')
        self.sfv[self.sfv<-1000000]=-100000
        self.sfv[self.sfv>1000000]=100000
        return self.sfv
            
    def getmd(self,par,parvalue):#should not be useful in any situation.
        ref=np.array(0,dtype=np.float32)
        for k in range(0,len(self.stypes)):
            stype=self.stypes[k]
            if len(par)>0:
                spar=par[self.spos[k]:self.spos[k+1]]
            else:
                spar=[]
            sparvalue=parvalue[self.spos[k]:self.spos[k+1]]
            reft=self.get1d(stype,spar,sparvalue)
            reft=reft.astype(np.float32)
            ref=ref.reshape(list(ref.shape)+list(np.ones(len(reft.shape))))+(reft.reshape(list(np.ones(len(ref.shape)))+list(reft.shape)))
        return ref

    def get1d(self,ptype,par=[],parvalue=[]):
        #remember to update initial value generation routine
        pi=self.types.index(ptype)
        if ptype=='dfire':
            return parvalue[0]*(np.log(self.bins[pi])-np.log(self.bins[pi][-1]))
        elif ptype.startswith('pl'):
            res=numpy.polynomial.polynomial.polyval(self.bins[pi],parvalue)
            if ptype.endswith('s'):
                res=res/res.sum()
            else:
                res=res/res[-1]            
            return np.log(np.abs(res))
        elif ptype.startswith('rwr'):
            n=max(int(parvalue[0]*200),1)
            la=parvalue[1]
            po=parvalue[2]
            r2=np.reshape(np.array(self.bins[pi])**2,[len(self.bins[pi]),1])
            ni=np.reshape(1.0/np.arange(1,n+1),[1,n])
            res=(r2**po)*np.sum(np.exp(-la*r2*ni)*(ni**1.5),axis=1)
            if ptype.endswith('s'):
                res=res/res.sum()
            else:
                res=res/1.0
            return np.log(np.abs(res))  
        elif ptype.startswith('rw'): #random walk reference state from Yang zhang.
            n=int(parvalue[0])
            la=parvalue[1]
            r2=np.reshape(np.array(self.bins[pi])**2,[len(self.bins[pi]),1])
            ni=np.reshape(1/np.arange(1,n+1),[1,n])
            res=r2*np.sum(np.exp(-la*r2*ni)*(ni**1.5),axis=1)
            if ptype.endswith('s'):
                res=res/res.sum()
            else:
                res=res/res[-1]
            return np.log(np.abs(res))          
        elif ptype.startswith('normb'):
            n=len(parvalue)/3 # mean ,std, weight
            res=np.zeros(len(self.bins[pi]))
            for i in range(n):
                res+=parvalue[i*3+2]*scipy.stats.norm.pdf(self.bins[pi],parvalue[i*3],np.abs(parvalue[i*3+1]))
            if ptype.endswith('s'):
                res=res/res.sum()
            else:
                res=res/1.0            
            return np.log(np.abs(res))  
        elif ptype.startswith('lognb'):
            n=len(parvalue)/3 # mean ,std, weight
            res=np.zeros(len(self.bins[pi]))
            for i in range(n):
                res+=parvalue[i*3+2]*scipy.stats.lognorm.pdf(self.bins[pi],np.abs(parvalue[i*3]),np.abs(parvalue[i*3+1]))
            if ptype.endswith('s'):
                res=res/res.sum()
            else:
                res=res/1.0
            return np.log(np.abs(res+0.0000001))  
        elif ptype.startswith('sinb'):
            n=len(parvalue)/3 # mean ,std, weight
            res=np.zeros(len(self.bins[pi]))
            for i in range(n):
                res+=parvalue[i*3+2]*np.sin(parvalue[i*3+1]*(np.array(self.bins[pi])-parvalue[i*3]))
            if ptype.endswith('s'):
                res=res/res.sum()
            else:
                res=res/1.0           
            return np.log(np.abs(res)) 
        elif ptype.startswith('expb'):
            n=len(parvalue)/4 # mean ,std, weight,exponential
            res=np.zeros(len(self.bins[pi]))
            for i in range(n):
                res+=parvalue[i*4+3]*np.exp(-parvalue[i*4+2]*np.abs(np.array(self.bins[pi])-parvalue[i*4])**parvalue[i*4+1])
            if ptype.endswith('s'):
                res=res/res.sum()
            else:
                res=res/1.0            
            return np.log(np.abs(res))         
        elif ptype.startswith('rexpb'):
            n=len(parvalue)/5 # mean ,std, weight
            res=np.zeros(len(self.bins[pi]))
            for i in range(n):
                res+=(np.array(self.bins[pi])**parvalue[i*5+4])*parvalue[i*5+3]*np.exp(-parvalue[i*5+2]*np.abs(np.array(self.bins[pi])-parvalue[i*5])**parvalue[i*5+1])
            if ptype.endswith('s'):
                res=res/res.sum()
            else:
                res=res/1.0            
            return np.log(np.abs(res))   
        elif ptype.startswith('svdb'):
            k=0
            for pvs in parvalue:
                if k==0:
                    res=self.svdu[:,k]*pvs
                else:
                    res+=self.svdu[:,k]*pvs
                k+=1                
            if ptype.endswith('s'):
                res=res/res.sum()
            else:
                res=res/1.0          
            return res
        elif ptype.startswith('spline'):
            us=scipy.interpolate.InterpolatedUnivariateSpline(par,parvalue)
            return np.log(np.abs(us(self.bins[pi]))+0.00000000000000000001)
        elif ptype.startswith('pspline'):
            us=scipy.interpolate.InterpolatedUnivariateSpline(par,np.abs(parvalue))
            return np.log(np.abs(us(self.bins[pi]))+0.00000000000000000001)
        elif ptype.startswith('s'):
            pl=list(par)
            newparvalue=list(parvalue)
            #whether par and parvalue are all specified #a - all specified
            pl=[self.bins[pi][0]]+pl+[self.bins[pi][-1]]
            if ptype[1]!='a':
                newparvalue=newparvalue+[int(ptype[1])]
            par=np.array(pl)
            #how to convert parvalue to control point values
            if ptype[2]=='o':
                parvalue=np.array(newparvalue)
            elif ptype[2]=='c':
                parvalue=np.cumsum(newparvalue)
            if ptype[3]=='n':
                if ptype[4] in ['l','s']:
                    parvalue=parvalue-parvalue[-1]
                else:
                    parvalue=parvalue/parvalue[-1]
            #how to generate the splines
            print par
            print parvalue
            if ptype[4]=='e':
                us=scipy.interpolate.pchip(par,np.log(parvalue))
                svalue=us(self.bins[pi])
            elif ptype[4]=='o':
                us=scipy.interpolate.pchip(par,parvalue)           
                svalue=np.log(np.abs(us(self.bins[pi]))+0.00000000000000000001)
            elif ptype[4]=='l':
                us=scipy.interpolate.pchip(par,parvalue)           
                svalue=us(self.bins[pi])
            elif ptype[4]=='f':
                us=scipy.interpolate.InterpolatedUnivariateSpline(par,np.log(np.abs(parvalue)))
                svalue=us(self.bins[pi])               
            elif ptype[4]=='p':
                us=scipy.interpolate.InterpolatedUnivariateSpline(par,parvalue)           
                svalue=np.log(np.abs(us(self.bins[pi]))+0.00000000000000000001)
            elif ptype[4]=='s':
                us=scipy.interpolate.InterpolatedUnivariateSpline(par,parvalue)           
                svalue=us(self.bins[pi])
            return svalue
        elif ptype.startswith('psn4'): #
            #positive monotonic spline, not normailized, search interval
            pl=list(par)
            newparvalue=list(parvalue)
            pl.append(self.bins[pi][-1])
            par=np.array(pl)    
            parvalue=np.cumsum(newparvalue)
            us=scipy.interpolate.pchip(par,np.abs(parvalue))
            return np.log(np.abs(us(self.bins[pi]))+0.00000000000000000001) 
        elif ptype.startswith('psn3'):
            #positive monotonic spline, normailized at last bin=1, search interval
            #print par
            #print parvalue
            pl=list(par)
            newparvalue=list(parvalue)+[1]
            pl.append(self.bins[pi][-1])
            par=np.array(pl)    
            parvalue=np.cumsum(newparvalue)
            parvalue=parvalue/parvalue[-1]
            #print par 
            #print parvalue
            #valuelist=pypchip.pchip(par,np.abs(parvalue),self.bins[pi])
            #print self.bins[pi]
            #print valuelist
            #return np.log(np.abs(valuelist)+0.00000000000000000001)
            us=scipy.interpolate.pchip(par,np.abs(parvalue))
            return np.log(np.abs(us(self.bins[pi]))+0.00000000000000000001)             
            #us=scipy.interpolate.InterpolatedUnivariateSpline(par,np.log(np.abs(parvalue)))
            #return us(self.bins[pi])
        elif ptype.startswith('psn2'):
            #positive  spline, not normalized, 
            #print par
            #print parvalue
            pl=list(par)
            newparvalue=list(parvalue)
            pl.append(self.bins[pi][-1])
            par=np.array(pl)    
            parvalue=np.array(newparvalue)
            #print par
            #print parvalue
            us=scipy.interpolate.InterpolatedUnivariateSpline(par,np.abs(parvalue))
            return np.log(np.abs(us(self.bins[pi]))+0.00000000000000000001) 
        elif ptype.startswith('psn'):
            #positive spline, normailized at last bin=1, 
            #print par
            #print parvalue
            pl=list(par)
            newparvalue=list(parvalue)
            pl.append(self.bins[pi][-1])
            newparvalue.append(1)
            par=np.array(pl)    
            parvalue=np.array(newparvalue)
            #print par
            #print parvalue
            us=scipy.interpolate.InterpolatedUnivariateSpline(par,np.abs(parvalue))
            return np.log(np.abs(us(self.bins[pi]))+0.00000000000000000001)       
        elif ptype.startswith('ps'):
            #print par
            #print parvalue
            pl=list(par)
            newparvalue=list(parvalue)
            par=np.array(pl)    
            parvalue=np.array(newparvalue)
            #print par
            #print parvalue
            us=scipy.interpolate.InterpolatedUnivariateSpline(par,np.abs(parvalue))
            return np.log(np.abs(us(self.bins[pi]))+0.00000000000000000001) 
        elif ptype=='ig':
            xi=np.arange(0.5,30,1)
            aprpr=[   144.434703,   1412.820011,   3834.829326,   7202.871092,
             11341.975029,  16096.371468,  21329.303897,  26919.548441,
             32742.984165,  38692.779078,  44664.067503,  50566.74941 ,
             56312.288375,  61810.089186,  66974.927745,  71737.414983,
             76029.089773,  79803.212568,  83024.246634,  85671.647628,
             87740.847002,  89233.801466,  90173.44354 ,  90561.097026,
             90429.358327,  89795.74556 ,  88703.604215,  87182.984942,
             85264.625698,  83004.011034]
            us=scipy.interpolate.InterpolatedUnivariateSpline(xi*parvalue[0],aprpr)
            ref=us(self.bins[pi])
            ref=np.abs(ref/ref[-1])
            return np.log(ref)
        elif ptype.startswith('logbins'):
            return parvalue
        elif ptype.startswith('bins'):
            sfv=np.array(parvalue)#np.log(np.abs(parvalue)+0.000000000000000000001)
            return sfv
        else:
            raise Bugs('Unknown reference function type: '+ptype)

