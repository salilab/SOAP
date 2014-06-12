import os
import re
import pickle
import pdb
from modeller import *


bd='/bell2/gqdong/rawdecoyfiles/ppd/ppd4/'

def complete_pdbfile(pdbfile):
    env=environ()
    env.io.hydrogen=True
    env.io.hetatm=False    
    env.libs.topology.read('${LIB}/top_heav.lib')
    env.libs.parameters.read('${LIB}/par.lib')
    m =model(env)
    m.read(pdbfile)
    aln=alignment(env)
    aln.append_model(m,atom_files=pdbfile,align_codes='test')
    m2=model(env)
    m2.generate_topology(aln[0])
    m2.read(pdbfile)
    m2.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES') 
    m2.write(file=pdbfile,no_ter=True)

def replace_end(lpdb):
    fh=open(lpdb)
    fc=fh.read()
    nofe=fc.count('END')
    fh.close()
    if nofe!=1:
        print "more than one END in model"
        pdb.set_trace()
    else:
        fc=fc.replace('END','TER')
    fh=open(lpdb,'w')
    fh.write(fc)
    fh.close()




pdbdir='/bell2/gqdong/rawdecoyfiles/ppd/originaldata/Benchmark4.0/'
for dirpath,dirnames,filenames in os.walk(pdbdir):
    os.chdir(dirpath)
    for f in filenames:
        print f
        #print os.system('cp '+f+' '+pdbdir)
        
fl=os.listdir(bd)
wd={}
for d in fl:
    print 'ppd4 '+d
    sd={}
    wd[d]=sd
    os.chdir(bd+d)
    print os.system('rm '+d+'.pickle')
    print os.system('rm *base.pdb')
    try:
        fh=open('docking.res')
    except:
        pdb.set_trace()
    fc=fh.read()
    fh.close()
    print os.system('cp '+pdbdir+d+'* ./')
    rpdb=re.search('receptorPdb\t\(Str\)\t(.*pdb.HB)\n',fc).group(1)
    lpdb=re.search('ligandPdb\t\(Str\)\t(.*pdb.HB)\n',fc).group(1)
    if rpdb[5]=='r':
        sd['orpdb']=rpdb
        sd['olpdb']=lpdb
        sd['switched']='False'
    elif rpdb[5]=='l':
        sd['orpdb']=lpdb
        sd['olpdb']=rpdb
        sd['switched']='True'
    else:
        pdb.set_trace()
        raise Exception('Do not know how to handle pdbs')     
    env=environ()
    env.io.hetatm=False
    env.io.hydrogen=True
    m1=model(env)
    m2=model(env)
    m1.read(sd['orpdb'])
    m2.read(sd['olpdb'])
    m1cl=''
    try:
        ci=0
        cnl=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n']
        for c1 in m1.chains:
            c1.name=cnl[ci]
            ci+=1
            m1cl=m1cl+c1.name
        m2cl=''
        for c2 in m2.chains:
            c2.name=cnl[ci]
            ci+=1
            m2cl=m2cl+c2.name
    except:
        pdb.set_trace()
    sd['orcl']=m1cl
    sd['olcl']=m2cl    
    m1.write(sd['orpdb'],no_ter=True)
    m2.write(sd['olpdb'],no_ter=True)
    m1=model(env)
    m2=model(env)
    print os.system('cp '+sd['orpdb']+' '+sd['orpdb']+'.AB')
    print os.system('cp '+sd['olpdb']+' '+sd['olpdb']+'.AB')
    sd['orpdbAB']=sd['orpdb']+'.AB'
    sd['olpdbAB']=sd['olpdb']+'.AB'
    print os.system('/bell1/home/gqdong/Download/PatchDock/mainchain.pl '+sd['orpdbAB']+' A')
    print os.system('/bell1/home/gqdong/Download/PatchDock/mainchain.pl '+sd['olpdbAB']+' B')     
    m1.read(sd['orpdbAB'])
    m2.read(sd['olpdbAB'])
    m1.write(sd['orpdbAB'],no_ter=True)
    m2.write(sd['olpdbAB'],no_ter=True)    
    complextype=re.search('receptorSeg\t\(Str\)\t(.*)\n',fc).group(1) #receptorSeg     (Str)   10.0 20.0 1.5 1 1 1 2
    if complextype[-1]=='0':
        ct='other'
    elif complextype[-1]=='2':
        ct='AA'
    elif complextype[-1]=='3':
        ct='EI'
    sd['rpdb']=rpdb
    sd['lpdb']=lpdb
    complete_pdbfile(rpdb)
    complete_pdbfile(lpdb)    
    sd['complextype']=ct
    fh2=open('firedockinput.pickle','w')
    pickle.dump(sd,fh2)
    fh2.close()
    replace_end(lpdb)
    replace_end(rpdb)
    replace_end(sd['orpdbAB'])
    replace_end(sd['olpdbAB'])
    print os.system('scp '+rpdb+' gqdong@baton2:~/decoys/ppd4/'+d+'/')
    print os.system('scp '+lpdb+' gqdong@baton2:~/decoys/ppd4/'+d+'/')
    print os.system('scp '+sd['orpdbAB']+' gqdong@baton2:~/decoys/ppd4/'+d+'/')
    print os.system('scp '+sd['olpdbAB']+' gqdong@baton2:~/decoys/ppd4/'+d+'/')    
    print os.system('scp firedockinput.pickle gqdong@baton2:~/decoys/ppd4/'+d+'/')

for key in wd:
    if wd[key]['complextype']!='AA' and wd[key]['complextype']=='l':
        print key+' '+wd[key]['complextype']+' '+wd[key]['rpdb'][5]

fh3=open('/bell2/gqdong/rawdecoyfiles/ppd/firedocktype.pickle','w')
pickle.dump(wd,fh3)
fh3.close()

#prepare the input files for zdock decoys

bd2='/bell2/gqdong/rawdecoyfiles/ppd/zdock/decoys/'

fl=os.listdir(bd2)
for d in fl:
    print 'zdock '+d
    dd=d.upper()
    sd=wd[dd]
    os.chdir(bd2+d)    
    sd['rpdb']=sd['orpdb']
    sd['lpdb']=sd['olpdb']
    sd['switched']='False'
    fh2=open('firedockinput.pickle','w')
    pickle.dump(sd,fh2)
    fh2.close()
    print os.system('scp '+bd+dd+'/'+sd['rpdb']+' gqdong@baton2:~/decoys/zdock/'+d+'/')
    print os.system('scp '+bd+dd+'/'+sd['lpdb']+' gqdong@baton2:~/decoys/zdock/'+d+'/')
    print os.system('scp '+bd+dd+'/'+sd['orpdbAB']+' gqdong@baton2:~/decoys/zdock/'+d+'/')
    print os.system('scp '+bd+dd+'/'+sd['olpdbAB']+' gqdong@baton2:~/decoys/zdock/'+d+'/')    
    print os.system('scp firedockinput.pickle gqdong@baton2:~/decoys/zdock/'+d+'/')

if 0:
    env=environ()
    env.io.hetatm=True
    env.io.hydrogen=True
    m1=model(env)
    m2=model(env)
    m1.read(sd['orpdb'])
    m2.read(sd['olpdb'])
    m1cl=''
    try:
        ci=0
        cnl=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n']
        for c1 in m1.chains:
            c1.name=cnl[ci]
            ci+=1
            m1cl=m1cl+c1.name
        m2cl=''
        for c2 in m2.chains:
            c2.name=cnl[ci]
            ci+=1
            m2cl=m2cl+c2.name
    except:
        pdb.set_trace()
    sd['orcl']=m1cl
    sd['olcl']=m2cl
    m1.write(sd['orpdb'],no_ter=True)
    m2.write(sd['olpdb'],no_ter=True)