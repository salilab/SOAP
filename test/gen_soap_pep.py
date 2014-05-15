"""
    Write out the SOAP tables for use in Modeller and IMP.
    
    Also write our raw scores for comparasion with scores from IMP or Modeller.
"""

from SOAP.scorer import *

import pdb

logpath=runenv.basedir+'results/mhc2/dcg/runs/57466-nan,nan_nan+nan_nan+nan_2297.581_0.000__nan+nan_nan+nan_0.000_0.000 /'
print logpath

bm=pickle.load(open(logpath+'cvmodel.pickle'))

so=scorer(model=bm['model'])


print so.assess(bm['allresult'][0]['bestpar'],slevel='dcg')
print so.assess_model()

so.write_potential(filetype='hdf5')
#np.save('allscore.npy',so.score)
jdso=pickle.load(open(runenv.basedir+'Decoys/mhc2/mhc2.pickle'))


rl=[', '.join([dso.dnlist[i],str(so.score[i]),str(so.scorearray[9,i]*so.ratioarray[9]),str(sum(so.scorearray[4:9,i]*so.ratioarray[4:9])),str(sum(so.scorearray[:4,i])+so.scorearray[10,i])]+map(str,list(so.scorerlist[-1].idist[i]))) for  i in range(len(dso.dnlist))]
open('scores','w').write('\n'.join(rl))
print os.getcwd()
pdb.set_trace()
