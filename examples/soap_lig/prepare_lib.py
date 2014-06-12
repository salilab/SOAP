import os
import pdb
import h5py

libpath=os.path.join(os.path.abspath('../../'),'lib')
dbl=open(os.path.join(libpath,'dicls-all.lib')).read()
ac=[i for i in open(os.path.join(libpath,'atmcls-plas.lib')).read().split('\n') if 'HET' in i]
ac=[i.split()[-1][1:-1] for i in ac]
db=[]

for a in ac:
    for b in ac:
        db.append('DBLGRP  "HET:{a}:{b}"\n  DOUBLET  "HET" "{a}" "{b}"'.format(a=a,b=b))

#t=h5py.File('/bell3/gqdong/statpot/ligX/tl400-cs1/ligX.tl400.cs1.hdf5')['mdt']
#pdb.set_trace()

open(os.path.join(libpath,'dicls-lig.lib'),'w').write(dbl+'\n'+'\n'.join(db))

pdb.set_trace()

