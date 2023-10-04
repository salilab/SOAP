from __future__ import print_function
from SOAP.decoys import *


#define the paths for the decoys
workingPath='/bell3/dina/rawdecoyfiles/ppd/'
originalDecoyPath='/bell3/dina/AAbench/'
preparedDecoyPath=os.path.join(workingPath,'aabench/')
decoySetName='aabench'

#copy the files from originalDecoyPath to target Path

fl=os.listdir(originalDecoyPath)
for f in fl:
    if not os.path.isdir(os.path.join(originalDecoyPath,f)):
        continue
    print(os.system('rm -rf '+os.path.join(preparedDecoyPath,f)))
    os.makedirs(os.path.join(preparedDecoyPath,f))
    os.chdir(os.path.join(preparedDecoyPath,f))
    print(os.system('cp '+os.path.join(originalDecoyPath,f,'docking.res')+' ./'))
    print(os.system('cp '+os.path.join(originalDecoyPath,f,'*pdb')+' ./'))
    print(os.system('touch needtransformation'))

#pre-prepare the decoys

def replace_end(lpdb):
    fh=open(lpdb)
    fc=fh.read()
    nofe=fc.count('END')
    fh.close()
    if nofe!=1:
        print("more than one END in model")
        pdb.set_trace()
    else:
        fc=fc.replace('END','TER')
    fh=open(lpdb,'w')
    fh.write(fc)
    fh.close()

bd=preparedDecoyPath

fl=os.listdir(bd)
wd={}
for d in fl:
    print(' '+d)
    sd={}
    wd[d]=sd
    os.chdir(bd+d)
    print(os.system('rm '+d+'.pickle'))
    print(os.system('rm *base.pdb'))
    try:
        fh=open('docking.res')
    except:
        print(os.system('rm -rf '+bd+d))
        continue
        pdb.set_trace()
    fc=fh.read()
    fh.close()
    rpdb=re.search('receptorPdb\t\(Str\)\t(.*pdb)\n',fc).group(1)
    lpdb=re.search('ligandPdb\t\(Str\)\t(.*pdb)\n',fc).group(1)
    if rpdb=='antigen.pdb':
        sd['orpdb']=rpdb
        sd['olpdb']=lpdb
        sd['switched']='False'
    else:
        pdb.set_trace()
        raise Exception('Do not know how to handle pdbs')
    env=Environ()
    env.io.hetatm=False
    env.io.hydrogen=True
    m1=Model(env)
    m2=Model(env)
    m1.read(sd['orpdb'])
    m2.read(sd['olpdb'])
    m1cl=''
    try:
        ci=0
        cnl=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z','a','b','c','d','e','f','g','h','i','j','k','l','m','n']
        for c1 in m1.chains:
            c1.name='A'
            ci+=1
            m1cl=m1cl+c1.name
        m2cl=''
        for c2 in m2.chains:
            c2.name='B'
            ci+=1
            m2cl=m2cl+c2.name
    except:
        pdb.set_trace()
    sd['orcl']=m1cl
    sd['olcl']=m2cl
    m1.write(sd['orpdb'],no_ter=True)
    m2.write(sd['olpdb'],no_ter=True)
    #m1=Model(env)
    #m2=Model(env)
    print(os.system('cp '+sd['orpdb']+' '+sd['orpdb']+'.AB'))
    print(os.system('cp '+sd['olpdb']+' '+sd['olpdb']+'.AB'))
    sd['orpdbAB']=sd['orpdb']+'.AB'
    sd['olpdbAB']=sd['olpdb']+'.AB'
    #print(os.system('/bell1/home/gqdong/Download/PatchDock/mainchain.pl '+sd['orpdbAB']+' A'))
    #print(os.system('/bell1/home/gqdong/Download/PatchDock/mainchain.pl '+sd['olpdbAB']+' B'))
    #m1.read(sd['orpdbAB'])
    #m2.read(sd['olpdbAB'])
    #m1.write(sd['orpdbAB'],no_ter=True)
    #m2.write(sd['olpdbAB'],no_ter=True)
    complextype=re.search('receptorSeg\t\(Str\)\t(.*)\n',fc).group(1) #receptorSeg     (Str)   10.0 20.0 1.5 1 1 1 2
    if complextype[-1]=='0':
        ct='other'
    elif complextype[-1]=='2':
        ct='AA'
    elif complextype[-1]=='3':
        ct='EI'
    sd['rpdb']=rpdb
    sd['lpdb']=lpdb
    #complete_pdbfile(rpdb)
    #complete_pdbfile(lpdb)
    sd['complextype']=ct
    fh2=open('firedockinput.pickle','wb')
    pickle.dump(sd,fh2)
    fh2.close()
    replace_end(lpdb)
    replace_end(rpdb)
    replace_end(sd['orpdbAB'])
    replace_end(sd['olpdbAB'])

#the decoys are now in the format that SOAP can process, use SOAP.decoys to build the decoyset
do=DecoySet(dsname=decoySetName, sourcedir=preparedDecoyPath)
do.build()
pdb.set_trace()
