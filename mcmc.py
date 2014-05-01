"""
   SOAP Markov Chain Monte Carlo

"""
from env import *

class MyMCMC(MCMC):
    def restore_sm_state(self):
        pass

class MMetropolis(Metropolis): #discrete anchor points and lower limit on the scale factor
    def __init__(self, stochastic, blockind=-999, scale=1., proposal_sd=None, proposal_distribution=None, verbose=1, tally=True, check_before_accepting=True):
        # Metropolis class initialization
        # Initialize superclass
        StepMethod.__init__(self, [stochastic], tally=tally)
        self.blockind=blockind
        # Initialize hidden attributes
        self.adaptive_scale_factor = 1.
        self.accepted = 0.
        self.rejected = 0.
        self._state = ['rejected', 'accepted', 'adaptive_scale_factor', 'proposal_sd', 'proposal_distribution', 'check_before_accepting']
        self._tuning_info = ['adaptive_scale_factor']
        self.check_before_accepting = check_before_accepting

        # Set public attributes
        self.stochastic = stochastic
        if verbose is not None:
            self.verbose = verbose
        else:
            self.verbose = stochastic.verbose

        if proposal_distribution != "Prior":
            # Avoid zeros when setting proposal variance
            if proposal_sd is not None:
                self.proposal_sd = proposal_sd
            else:
                if all(self.stochastic.value != 0.):
                    self.proposal_sd = np.ones((self.stochastic.value.shape)) * np.abs(self.stochastic.value) * scale
                else:
                    self.proposal_sd = ones(shape(self.stochastic.value)) * scale

            # Initialize proposal deviate with array of zeros
            self.proposal_deviate = np.zeros((self.stochastic.value.shape), dtype=float)

            # Determine size of stochastic
            if isinstance(self.stochastic.value, np.ndarray):
                self._len = len(self.stochastic.value.ravel())
            else:
                self._len = 1
        self.proposal_distribution='Normal'
        self.blockind=blockind
        
    def propose(self):
        """
        This method is called by step() to generate proposed values
        if self.proposal_distribution is "Normal" (i.e. no proposal specified).
        """
        runenv.stepblockind=self.blockind
        if self.proposal_distribution == "Normal":
            self.stochastic.value = rnormal(self.stochastic.value, self.adaptive_scale_factor * self.proposal_sd, size=self.stochastic.value.shape)
        elif self.proposal_distribution == "Prior":
            self.stochastic.random()
        
class MyMetropolis(Metropolis):
    def __init__(self, stochastic, parindex, proposemethod, proposepars, blockind=-999, searchind=[],scale=1, proposal_sd=None, proposal_distribution=None, verbose=1, tally=True, check_before_accepting=True):
        # Metropolis class initialization
        # Initialize superclass
        StepMethod.__init__(self, [stochastic], tally=tally)
        self.blockind=blockind
        if self.blockind>=0:
            self.searchind=searchind
            parindexstart=parindex[searchind[0]]
            parindex=[par-parindexstart for par in parindex]
            proposemethod=proposemethod[searchind[0]:searchind[-1]+1]
            proposepars=proposepars[searchind[0]:searchind[-1]+1]
        else:
            self.searchind=range(len(proposemethod))
        
        # Initialize hidden attributes
        self.adaptive_scale_factor = 0.5
        self.accepted = 0.
        self.rejected = 0.
        self._state = ['rejected', 'accepted', 'adaptive_scale_factor', 'proposal_sd', 'proposal_distribution', 'check_before_accepting']
        self._tuning_info = ['adaptive_scale_factor']
        self.check_before_accepting = check_before_accepting

        # Set public attributes
        self.stochastic = stochastic
        if verbose is not None:
            self.verbose = verbose
        else:
            self.verbose = stochastic.verbose
        self.parindex=parindex
        self.proposemethod=proposemethod
        self.proposepars=proposepars
        self.proposevalue=np.zeros(self.stochastic.value.shape)
        for i in range(len(self.parindex)-1):
             if self.proposemethod[i] == "splineparn":
                numofpars=self.proposepars[i][2]
                parmean=self.proposepars[i][3]
                parstd=parmean*scale
                cov=np.zeros([numofpars-1,numofpars-1])+parstd**2/(numofpars)
                for j in range(numofpars-1):
                    cov[j,j]-=parstd**2
                self.proposepars[i].append(cov)

        if proposal_distribution != "Prior":
            # Avoid zeros when setting proposal variance
            if proposal_sd is not None:
                self.proposal_sd = proposal_sd
            else:
                if all(self.stochastic.value != 0.):
                    self.proposal_sd = np.ones((self.stochastic.value.shape)) * np.abs(self.stochastic.value) * scale
                else:
                    self.proposal_sd = ones(shape(self.stochastic.value)) * scale

            # Initialize proposal deviate with array of zeros
            self.proposal_deviate = np.zeros((self.stochastic.value.shape), dtype=float)

            # Determine size of stochastic
            if isinstance(self.stochastic.value, np.ndarray):
                self._len = len(self.stochastic.value.ravel())
            else:
                self._len = 1
        self.proposal_distribution='Normal'
    
    def tune(self, divergence_threshold=1e10, verbose=0):
        """
        Tunes the scaling parameter for the proposal distribution
        according to the acceptance rate of the last k proposals:

        Rate    Variance adaptation
        ----    -------------------
        <0.001        x 0.1
        <0.05         x 0.5
        <0.2          x 0.9
        >0.5          x 1.1
        >0.75         x 2
        >0.95         x 10

        This method is called exclusively during the burn-in period of the
        sampling algorithm.

        May be overridden in subclasses.
        """

        if self.verbose is not None:
            verbose = self.verbose

        if self.verbose is not None:
            verbose = self.verbose

        # Verbose feedback
        if verbose > 0:
            print '\t%s tuning:' % self._id
            print '\t\tadaptive scale factor(old):', self.adaptive_scale_factor
            
        # Flag for tuning state
        tuning = True

        # Calculate recent acceptance rate
        if not (self.accepted + self.rejected): return tuning
        acc_rate = self.accepted / (self.accepted + self.rejected)
        print self.accepted
        print self.rejected

        print "Tunning with acc_rate_ratio "+str(runenv.acc_rate_ratio)
        acc_rate=acc_rate*runenv.acc_rate_ratio
        # Switch statement
        if acc_rate<0.001:
            # reduce by 90 percent
            self.adaptive_scale_factor *= 0.1
        elif acc_rate<0.05:
            # reduce by 50 percent
            self.adaptive_scale_factor *= 0.5
        elif acc_rate<0.2:
            # reduce by ten percent
            self.adaptive_scale_factor *= 0.9
        elif acc_rate>0.95:
            # increase by factor of ten
            self.adaptive_scale_factor *= 10.0
        elif acc_rate>0.75:
            # increase by double
            self.adaptive_scale_factor *= 2.0
        elif acc_rate>0.5:
            # increase by ten percent
            self.adaptive_scale_factor *= 1.1
        else:
            tuning = False

        if self.adaptive_scale_factor>0.5:
            self.adaptive_scale_factor=0.5
        elif self.adaptive_scale_factor<0.00001:
            self.adaptive_scale_factor=0.00001

        # Re-initialize rejection count
        self.rejected = 0.
        self.accepted = 0.

        # More verbose feedback, if requested
        if verbose > 0:
            if hasattr(self, 'stochastic'):
                print '\t\tvalue:', self.stochastic.value
            print '\t\tacceptance rate:', acc_rate
            print '\t\tadaptive scale factor:', self.adaptive_scale_factor
            print
        return tuning

    def propose(self):
        """
        This method is called by step() to generate proposed values
        if self.proposal_distribution is "Normal" (i.e. no proposal specified).
        """
        runenv.stepblockind=self.blockind
        pvalue=np.copy(self.proposevalue)
        for i in range(len(self.proposemethod)):
            if self.proposemethod[i] == "pnormal":
                while True:
                    value=np.random.normal(self.stochastic.value[self.parindex[i]:self.parindex[i+1]], self.adaptive_scale_factor * self.proposal_sd[self.parindex[i]:self.parindex[i+1]])
                    if np.all(value>0):
                        break
                pvalue[self.parindex[i]:self.parindex[i+1]]=value
                #print value
            elif self.proposemethod[i] == "splineparn":
                while True:
                    sv=self.stochastic.value[self.parindex[i]:self.parindex[i+1]]
                    svr=np.copy(sv)
                    svr[:-1]=sv[1:]-sv[:-1]
                    value=np.random.multivariate_normal(svr[:-1],self.proposepars[i][-1]*(self.adaptive_scale_factor**2))
                    if np.any(value<=0):
                        continue
                    parvalue=value.cumsum()+self.proposepars[i][0]
                    if parvalue[-1]<self.proposepars[i][1]:
                        break
                #print value
                #print parvalue
                pvalue[self.parindex[i]+1:self.parindex[i+1]]=parvalue
                pvalue[self.parindex[i]]=self.proposepars[i][0]
        self.stochastic.value=pvalue
        
class MxMetropolis(MyMetropolis): #discrete anchor points
    def __init__(self, stochastic, parindex, proposemethod, proposepars, blockind=-999, searchind=[],scale=1., proposal_sd=None, proposal_distribution=None, verbose=mcverbose, tally=True, check_before_accepting=True):
        # Metropolis class initialization
        # Initialize superclass
        StepMethod.__init__(self, [stochastic], tally=tally)
        self.blockind=blockind
        if self.blockind>=0:
            self.searchind=searchind
            parindexstart=parindex[searchind[0]]
            parindex=[par-parindexstart for par in parindex]
            parindex=parindex[searchind[0]:searchind[-1]+2]
            proposemethod=proposemethod[searchind[0]:searchind[-1]+1]
            proposepars=proposepars[searchind[0]:searchind[-1]+1]
        else:
            self.searchind=range(len(proposemethod))
            
        # Initialize hidden attributes
        self.adaptive_scale_factor = 1.
        self.accepted = 0.
        self.rejected = 0.
        self._state = ['rejected', 'accepted', 'adaptive_scale_factor', 'proposal_sd', 'proposal_distribution', 'check_before_accepting']
        self._tuning_info = ['adaptive_scale_factor']
        self.check_before_accepting = check_before_accepting

        # Set public attributes
        self.stochastic = stochastic
        if verbose is not None:
            self.verbose = verbose
        else:
            self.verbose = stochastic.verbose
        self.parindex=parindex
        self.proposemethod=proposemethod
        self.proposepars=proposepars
        print self.proposepars
        self.proposevalue=np.zeros(self.stochastic.value.shape)
        for i in range(len(self.parindex)-1):
             if self.proposemethod[i] == "splineparn":
                numofpars=self.proposepars[i][2]
                parmean=self.proposepars[i][3]-self.proposepars[i][4]
                parstd=parmean*scale
                cov=np.zeros([numofpars,numofpars])+parstd**2/(numofpars+1)
                for j in range(numofpars):
                    cov[j,j]-=parstd**2
                self.proposepars[i].append(cov)
        if proposal_distribution != "Prior":
            # Avoid zeros when setting proposal variance
            if proposal_sd is not None:
                self.proposal_sd = proposal_sd
            else:
                if all(self.stochastic.value != 0.):
                    self.proposal_sd = np.ones((self.stochastic.value.shape)) *np.abs(self.stochastic.value) * scale
                else:
                    self.proposal_sd = ones(shape(self.stochastic.value)) * scale

            # Initialize proposal deviate with array of zeros
            self.proposal_deviate = np.zeros((self.stochastic.value.shape), dtype=float)

            # Determine size of stochastic
            if isinstance(self.stochastic.value, np.ndarray):
                self._len = len(self.stochastic.value.ravel())
            else:
                self._len = 1
        self.proposal_distribution='Normal'
    
    def propose(self):
        """
        This method is called by step() to generate proposed values
        if self.proposal_distribution is "Normal" (i.e. no proposal specified).
        """
        runenv.stepblockind=self.blockind
        pvalue=np.copy(self.proposevalue)
        #print self.proposal_sd
        for i in range(len(self.proposemethod)):
            if self.proposemethod[i] == "pnormal":
                #print "pnormal"
                #print self.stochastic.value[self.parindex[i]:self.parindex[i+1]]
                #print self.proposal_sd[self.parindex[i]:self.parindex[i+1]]
                value=np.random.normal(list(self.stochastic.value[self.parindex[i]:self.parindex[i+1]]), list(self.adaptive_scale_factor * self.proposal_sd[self.parindex[i]:self.parindex[i+1]]))
                #print value
                while np.any(value<0):
                    boolfilter=(value<0)
                    partvalue=np.random.normal(self.stochastic.value[self.parindex[i]:self.parindex[i+1]][boolfilter], self.adaptive_scale_factor * self.proposal_sd[self.parindex[i]:self.parindex[i+1]][boolfilter])
                    value[boolfilter]=partvalue
                    #print value
                pvalue[self.parindex[i]:self.parindex[i+1]]=value
                #print value
            elif self.proposemethod[i] == "normal":
                #print "pnormal"
                #print self.stochastic.value[self.parindex[i]:self.parindex[i+1]]
                #print self.proposal_sd[self.parindex[i]:self.parindex[i+1]]
                value=np.random.normal(list(self.stochastic.value[self.parindex[i]:self.parindex[i+1]]), list(self.adaptive_scale_factor * self.proposal_sd[self.parindex[i]:self.parindex[i+1]]))
                pvalue[self.parindex[i]:self.parindex[i+1]]=value
                #print value
            elif self.proposemethod[i] == "splineparn":
                #print "splineparn"
                #print self.stochastic.value[self.parindex[i]:self.parindex[i+1]]
                #print self.proposepars[i][-1]
                binsize=self.proposepars[i][4]
                k=0
                success=False
                while k<200:
                    k+=1
                    sv=self.stochastic.value[self.parindex[i]:self.parindex[i+1]]
                    svr=np.array([self.proposepars[i][0]]+list(sv))
                    svr[:-1]=svr[1:]-svr[:-1]
                    svr=svr-binsize
                    value=np.random.multivariate_normal(svr[:-1],self.proposepars[i][-1]*(self.adaptive_scale_factor**2))
                    #print value
                    if np.any(value<=0):
                        continue
                    value=value+binsize
                    parvalue=value.cumsum()+self.proposepars[i][0]
                    if parvalue[-1]<(self.proposepars[i][1]-binsize):
                        success=True
                        break
                #print value
                #print parvalue
                if success:
                    parvalue=np.round((parvalue-self.proposepars[i][0])/binsize)*binsize+self.proposepars[i][0]
                    pvalue[self.parindex[i]:self.parindex[i+1]]=parvalue
                else:
                    pvalue[self.parindex[i]:self.parindex[i+1]]=self.stochastic.value[self.parindex[i]:self.parindex[i+1]]
                #print parvalue
                #pvalue[self.parindex[i]]=self.proposepars[i][0]
        self.stochastic.value=pvalue
 
class MzMetropolis(MxMetropolis): #discrete anchor points and lower limit on the scale factor
    def __init__(self, stochastic, parindex, proposemethod, proposepars, blockind=-999,scale=1., proposal_sd=None, proposal_distribution=None, verbose=1, tally=True, check_before_accepting=True):
        # Metropolis class initialization
        # Initialize superclass
        MxMetropolis.__init__(self, stochastic, parindex, proposemethod, proposepars,blockind, scale, proposal_sd, proposal_distribution, verbose, tally, check_before_accepting)
        self.scale_factor_min=[]
        for i in range(len(self.parindex)-1):
            if self.proposemethod[i] == "splineparn":
                parmean=self.proposepars[i][3]-self.proposepars[i][4]
                self.scale_factor_min.append(0.5*(self.proposepars[i][4]/(scale*parmean)))
            else:
                self.scale_factor_min.append(0.001)
   
    def propose(self):
        """
        This method is called by step() to generate proposed values
        if self.proposal_distribution is "Normal" (i.e. no proposal specified).
        """
        runenv.stepblockind=self.blockind
        pvalue=np.copy(self.proposevalue)
        for i in range(len(self.proposemethod)):
            if self.adaptive_scale_factor<self.scale_factor_min[i]:
                print "using minscalefactor"
                scalefactor=self.scale_factor_min[i]
            else:
                scalefactor=self.adaptive_scale_factor
            if self.proposemethod[i] == "pnormal":
                while True:
                    value=np.random.normal(self.stochastic.value[self.parindex[i]:self.parindex[i+1]], scalefactor * self.proposal_sd[self.parindex[i]:self.parindex[i+1]])
                    if np.all(value>0):
                        break
                pvalue[self.parindex[i]:self.parindex[i+1]]=value
                #print value
            elif self.proposemethod[i] == "splineparn":
                binsize=self.proposepars[i][4]
                while True:
                    sv=self.stochastic.value[self.parindex[i]:self.parindex[i+1]]
                    svr=np.copy(sv)
                    svr[:-1]=sv[1:]-sv[:-1]
                    svr=svr-binsize
                    value=np.random.multivariate_normal(svr[:-1],self.proposepars[i][-1]*(scalefactor**2))
                    if np.any(value<=0):
                        continue
                    value=value+binsize
                    parvalue=value.cumsum()+self.proposepars[i][0]
                    if parvalue[-1]<(self.proposepars[i][1]-binsize):
                        break
                #print value
                #print parvalue
                parvalue=np.round((parvalue-self.proposepars[i][0])/binsize)*binsize+self.proposepars[i][0]
                #print parvalue
                pvalue[self.parindex[i]+1:self.parindex[i+1]]=parvalue
                pvalue[self.parindex[i]]=self.proposepars[i][0]
        self.stochastic.value=pvalue   
   