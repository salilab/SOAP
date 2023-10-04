from modeller import *
import mdt
import mdt.features
import os
import sys
import pdb

fl=os.listdir(sys.argv[1])
pfl=[f for f in fl if f.endswith('.pdb')]
scriptpath=os.path.dirname(os.path.abspath(__file__))
env=Environ()
env.io.atom_files_directory=[sys.argv[1]]
mlib=mdt.Library(env)
#mlib.atom_classes.read('$(LIB)/atmcls-mf.lib')
mlib.bond_classes.read('${LIB}/bndgrp.lib')
mlib.tuple_classes.read(os.path.join(scriptpath,'dicls-all.lib'))

d=mdt.features.TupleDistance(mlib,bins=mdt.uniform_bins(35, 0, 0.2))
a1=mdt.features.TupleAngle1(mlib,bins=mdt.uniform_bins(6, 0, 30))
a2=mdt.features.TupleAngle2(mlib,bins=mdt.uniform_bins(6, 0, 30))
h=mdt.features.TupleDihedral1(mlib,bins=mdt.uniform_bins(12, -180, 30))
t1=mdt.features.TupleType(mlib)
t2=mdt.features.TupleType(mlib, pos2=True)
m=mdt.Table(mlib,features=(d,a1,a2,h,t1,t2),bin_type=mdt.Float,shape=[-1]*6)
m.read_hdf5(os.path.join(scriptpath,os.path.join(scriptpath,'mdt1.hdf5')))

scorelist=[]
for pf in pfl:
    mdl=Model(env)
    mdl.read(file=pf)
    aln=Alignment(env)
    aln.append_model(mdl,align_codes='p',atom_files=pf)
    sr1=m.open_alignment(aln)
    score1=sr1.sum(chain_span_range=(-9999,0,0,9999),residue_span_range=(-9999,0,0, 9999),bond_span_range=(20,999999999),disulfide=True)
    scorelist.append(pf+': '+str(score1))
    
open(os.path.join(sys.argv[1],'soap_protein_od.score'),'w').write('\n'.join(scorelist))

pdb.set_trace()
