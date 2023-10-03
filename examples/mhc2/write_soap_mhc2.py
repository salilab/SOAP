"""
    Write out the SOAP tables for use in Modeller and IMP.
    
    Also write our raw scores for comparison with scores from IMP or Modeller.
"""

from __future__ import print_function
from SOAP.scorer import *

import pdb

#the path to the optimal potential
logpath=runenv.basedir+'results/mhc2/dcg/runs/57466-nan,nan_nan+nan_nan+nan_2297.581_0.000__nan+nan_nan+nan_0.000_0.000 /'
print(logpath)

bm=pickle.load(open(logpath+'cvmodel.pickle'))

so=scorer(model=bm['model'])

#assess the model using the best parameter (necessary for writing the optimal potential)
print(so.assess(bm['allresult'][0]['bestpar'],slevel='dcg'))
print(so.assess_model())

#write out the potential in HDF5 format; the path to the HDF5 files are printed out in log 
so.write_potential(filetype='hdf5')

#write out detailed score for comparison with other implementation of SOAP_mhc2
np.save('allscore.npy',so.score)
jdso=pickle.load(open(runenv.basedir+'Decoys/mhc2/mhc2.pickle'))
rl=[', '.join([dso.dnlist[i],str(so.score[i]),str(so.scorearray[9,i]*so.ratioarray[9]),str(sum(so.scorearray[4:9,i]*so.ratioarray[4:9])),str(sum(so.scorearray[:4,i])+so.scorearray[10,i])]+map(str,list(so.scorerlist[-1].idist[i]))) for  i in range(len(dso.dnlist))]
open('scores','w').write('\n'.join(rl))
print(os.getcwd())
#enter debug mode, so that you can examine the detail scores if necessary
pdb.set_trace()
