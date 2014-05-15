from modeller import *
import mdt
import mdt.features
import os
import sys

fl=os.listdir(sys.argv[1])
pfl=[f for f in fl if f.endswith('.pdb')]
scriptpath=os.path.dirname(os.path.abspath(__file__))
env=environ()
env.io.atom_files_directory=[sys.argv[1]]
mlib=mdt.Library(env)


mlib.atom_classes.read('$(LIB)/atmcls-mf.lib')
ub=mdt.uniform_bins(30, 0, 0.5)
d=mdt.features.AtomDistance(mlib,bins=ub)
a1=mdt.features.AtomType(mlib)

a2=mdt.features.AtomType(mlib, pos2=True)


m=mdt.Table(mlib,features=(d,a1,a2))
m.read_hdf5(os.path.join(scriptpath,'mdt1.hdf5'))
mlib2=mdt.Library(env)


mlib2.atom_classes.read('$(LIB)/atmcls-mf.lib')
b=mdt.features.FractionalAtomAccessibility(mlib2,bins=mdt.uniform_bins(10, 0, 0.1))
a3=mdt.features.AtomType(mlib2)

m2=mdt.Table(mlib2,features=(b,a3))
m2.read_hdf5(os.path.join(scriptpath,'mdt2.hdf5'))

scorelist=[]
for pf in pfl:
    mdl=model(env)
    mdl.read(file=pf)
    aln=alignment(env)
    aln.append_model(mdl,align_codes='p',atom_files=pf)
    sr1=m.open_alignment(aln)
    score1=sr1.sum(chain_span_range=(-99,-1,1,99),residue_span_range=(9999,0,0, 9999),bond_span_range=(-1,-1))
    sr2=m2.open_alignment(aln)
    score2=sr2.sum()
    scorelist.append(pf+': '+str(score1+score2))
    
open(os.path.join(sys.argv[1],'soap_pp.score'),'w').write('\n'.join(scorelist))