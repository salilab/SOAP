"""
   SOAP decoy set module.

"""


from __future__ import print_function
from env import *
import subprocess
from utility import mypickle
from sequences import *
import shutil
import glob
#sys.modules['sp']=sys.modules[__name__]

def cat_files(infiles, outfile):
    """Concatenate `infiles` to the file handle `outfile`"""
    for f in infiles:
        for line in open(f):
            outfile.write(line)

class Decoys4Single(object):
    """Single set of decoys, serving as the building block for :class:`DecoySet`, corresponding to a single subfolder in the :class:`DecoySet`'s folder. All the decoy structure files should present in this subfolder.

    This class contains method for converting the input decoy strctures in the subfolder into formats necessary for using the decoy structures for statistical potential training, and for copying the prepared decoys into the correponding folders in SGE cluster.


    :Parameters:
      - `code`: the name for this single set of decoys, which is also the subfolder name.
      - `dsname`: the :class:`DecoySet` name.
      - `path`: the path to the :class:`DecoySet`, the place to create the subfoler for this set of decoys.
      - `dtype`: list of numpy data types defining the strctured array (you should not change this)
      - `extrapar`: numpy data types for additonal propertites of each decoy structure is define here
        specially for atom class features (see :class:`features.AtomType`)"""

    def __init__(self,code='',dsname='',path=''):
        self.code=code
        self.dsname=dsname
        self.dir=path
        self.path = os.path.join(self.dir, self.code)
        self.score=[]
        self.dtype=[('decoyset','i1'),('decoyseqposn','i2'),('decoyposn','i4'),('rmsd','f4'),('score','f4'),('score2','f4')]
        self.extrapar=[]

    def rename(self):
        #rename decoy names for easy retrieve
        os.chdir(self.dir)
        fl=os.listdir('./')
        for file in fl:
            if file[-3:]!='pdb' and file[-6:]!='pdb.gz':
                continue
            if file[0:(2+len(self.code)+len(self.dsname))]==self.dsname+'.'+self.code+'.':
                continue
            elif file[0:len(self.code)+1]==self.code+'.':
                shutil.move(file, self.dsname+'.'+file)
            else:
                shutil.move(file, self.dsname+'.'+self.code+'.'+file)
        #print(os.system('for file in *pdb; do mv $file '+self.dsname+'.'+self.code+'.$file; done'))

    def filter_pdb(self):
        #filter the decoys to eliminate bad ones, used for filter HR decoys. Ideally, the input decoys need not be filtered.
        os.chdir(self.dir)
        rmsdlist=[]
        flist=os.listdir(self.dir)
        nativepattern='native.*pdb$'
        print('running'+self.code)
        nativef=[elem for elem in flist if re.search(nativepattern,elem)][0]
        fh=open(nativef)
        fhc=fh.read()
        atomlines=re.findall('\\n(ATOM.*)','\n'+fhc+'\n')
        atomlist=[]
        seqnuml=[]
        seqnamel=[]
        oldseqnum=0
        atomcoord={}
        atomlinelist=[]
        atomnamedict={}
        for atomline in atomlines:
            atomcode=atomline[12:16].strip()
            if (atomline[16]==' ' or atomline[16]=='A') and atomcode!='OXT' and (not atomcode[0] in ['H','Q']) and (not atomcode[0:2] in ['1H', '2H','3H']):
                atomlinelist.append(atomline)
                atomlist.append(atomline[22:26].strip()+atomline[17:20].strip()+atomline[12:16].strip())
                if int(atomline[22:26])!=oldseqnum:
                    seqnuml.append(int(atomline[22:26]))
                    seqnamel.append(atomline[17:20])
                    oldseqnum=int(atomline[22:26])
        atomlistlen=len(atomlist)
        atomcoord[nativef]=atomlinelist
        atomnamedict[nativef]=atomlist
        #fist filter the decoy structures to eliminate outliers
        otherf=[elem for elem in flist if re.search('pdb$', elem)]
        otherf.remove(nativef)
        posdiffl=[]
        atomdict={}
        otherffiltered=[]
        for otherprotein in otherf:
            fs=''
            try:
                fht=open(otherprotein,'r')
                fhtc=fht.read()
                tatomlines=re.findall('\\n(ATOM.*)','\n'+fhtc+'\n')
                oldseqnum=0
                tseqnuml=[]
                tseqnamel=[]
                for tatomline in tatomlines:
                    if int(tatomline[22:26])!=oldseqnum:
                        tseqnuml.append(int(tatomline[22:26]))
                        tseqnamel.append(tatomline[17:20])
                        oldseqnum=int(tatomline[22:26])
            except:
                print(tatomline)
                print('structure problem '+otherprotein)
                shutil.move(otherprotein, otherprotein[:-4]+'.atombad')
                continue
            posdiff=self.align_seq(seqnamel,seqnuml,tseqnamel,tseqnuml)
            if posdiff==-9999:
                #if the decoy structure can not be aligned with the native structure, mark it as a bad structure
                print('can not align with native '+otherprotein)
                shutil.move(otherprotein, otherprotein[:-4]+'.alnbad')
                continue
            posdiffl.append(posdiff)
            tatomlist=[]
            atomlinelist=[]
            for tatomline in tatomlines:
                atomnamestr=tatomline[12:16].strip()
                if atomnamestr=='OT1':
                    atomnamestr='O'
                    tatomline.replace('OT1','O  ')
                atomlinelist.append(tatomline)
                atomstr=str(int(tatomline[22:26])+posdiff)+tatomline[17:20].strip()+atomnamestr
                tatomlist.append(atomstr)
                if atomstr in atomdict:
                    atomdict[atomstr]=atomdict[atomstr]+1
                else:
                    atomdict[atomstr]=1
            atomcoord[otherprotein]=atomlinelist
            atomnamedict[otherprotein]=tatomlist
            otherffiltered.append(otherprotein)
        otherf=otherffiltered
        #calculate the concensus atom list of decoys
        tatomlist=[]
        minnum=min(int(len(otherf)*0.7),len(otherf)-3)
        for key in atomdict:
            if atomdict[key]>minnum:
                tatomlist.append(key)
        #build the final atom list based on both decoys and native structures
        newatomlist=[]
        for item in atomlist:
            if item in tatomlist:
                newatomlist.append(item)
            else:
                print('native structure and decoys uncommon atom '+item)
        if len(newatomlist)< len(atomlist)*0.90:
            print(nativef)
            print('too many missing atoms')
            pdb.set_trace()
        atomlist=newatomlist
        #now filter both the native and the decoy structure based on the list
        otherf.append(nativef)
        posdiffl.append(0)
        k=0
        for otherprotein in otherf:
            fs=''
            atomlines=atomcoord[otherprotein]
            atomnames=atomnamedict[otherprotein]
            for atom in atomlist:
                try:
                    index=atomnames.index(atom)
                except:
                    print(otherprotein+' missing atoms '+atom)
                    shutil.move(otherprotein, otherprotein[:-3]+atom+'.bad')
                    break
                fs=fs+atomlines[index]+'\n'
            with open(otherprotein,'w') as fht:
                fht.write(fs)
                fht.flush()

    def align_seq(self,seqnamel,seqnuml,tseqnamel,tseqnuml):
        k=0
        i=len(seqnamel)/3
        fail=1
        while i<len(seqnamel)-10:
            try:
                while sum(seqnuml[i:i+10])!=(seqnuml[i]*10+45):
                    i=i+1
                seq10=''.join(seqnamel[i:i+10])
                seqpos=seqnuml[i]
                tseqnames=[''.join(tseqnamel[j:j+10]) for j in range(0,len(tseqnamel)-10)]
                index=tseqnames.index(seq10)
                tseqpos=tseqnuml[index]
                posdiff=seqpos-tseqpos
                fail=0
                break
            except Exception as e:
                if k==0:
                    i=k
                i=i+1
                k=k+1
        if fail:
            posdiff=-9999
        return posdiff

    def copy2ddb(self):
        #copy decoys to jobserver
        tdir = os.path.join(runenv.serverUserPath, 'decoys', self.dsname,
                            self.code)
        subprocess.check_call(['gzip'] + glob.glob('*.pdb'))
        subprocess.check_call(['gunzip'] + glob.glob('*base*.gz'))
        pdblist=" ".join([item+'.pdb.gz' for item in self.dnlist])
        subprocess.check_call(['ssh', runenv.jobserver, 'mkdir', '-p', tdir])
        subprocess.check_call(['scp', '-q'] + glob.glob('*.pdb.gz') \
                              + glob.glob('*base*') \
                              + [runenv.jobserver+':'+tdir])

    def calc_rmsd(self,nativepattern='*native*pdb$'):
        """Calculate the RMSDs for the decoys
        It takes either the rosetta sc file, rmsd.pickle, or raw decoys files with a native structure (named \*native\*pdb).
        rmsd.pickle
        {decoyname:rmsdvalue} or {decoyname:[rmsdvalue,extra properties defined by extrapar]}
        """
        #change the native name to be native.pdb
        fl=os.listdir('./')
        sfl=[item for item in fl if item.endswith('sc')]
        if len(sfl)==1:
            return self.read_rmsd_rosettadock(sfl[0])
        elif len(sfl)>1:
            raise Exception('more than 1 score file, do not know how to handle')
        flist=os.listdir(self.dir)
        rmsdfl=[elem for elem in flist if re.search('rmsd',elem) and not re.search('matlab',elem)]
        print(rmsdfl)
        if len(rmsdfl)==0:
            self.calc_rmsd_fromfiles()
            return 0
        if rmsdfl[0]=='rmsd.pickle':
            return self.load_rmsd(rmsdfl[0])
        rmsdf=rmsdfl[0]
        with open(rmsdf) as rmsdfh:
            rmsdfc=rmsdfh.read()
        nativef=[elem for elem in flist if re.search(nativepattern,elem)][0]
        flist.remove(nativef)
        otherf=[elem for elem in flist if re.search('pdb$', elem)]
        rmsdlist=[]
        rmsdlist.append([nativef[:-4],0.0])
        for file in otherf:
            filenamelist=file.split('.')
            if self.dsname=='HR':
                rer=re.search('\s*'+filenamelist[-2]+'\s+([0-9\.]+)\s*\\n',rmsdfc)
            else:
                rer=re.search('.'.join(filenamelist[-3:])+'\s+is\s+([0-9\.]+)\s*\\n',rmsdfc)
                if not rer:
                    rer=re.search('.'.join(filenamelist[-2:])+'\s+is\s+([0-9\.]+)\s*\\n',rmsdfc)
            try:
                rmsdlist.append([file[:-4],float(rer.group(1))])
            except:
                pdb.set_trace()
        rmsdlist.sort(key=itemgetter(1))
        self.dnlist=[item[0] for item in rmsdlist]
        self.build_sa(rmsdlist)

    def build_sa(self,rmsdlist):
        print("building sa")
        self.sa=np.zeros(len(self.dnlist), dtype=self.dtype+self.extrapar)
        for i in range(0,len(self.dnlist)):
            self.sa['rmsd'][i]=rmsdlist[i][1]
        k=1
        for item in self.extrapar:
            for i in range(0,len(self.dnlist)):
                self.sa[item[0]][i]=rmsdlist[i][k+1]
            k=k+1
        if len(self.score)>0:
            self.sa['score']=self.score
        print("building sa finished")

    def read_rmsd_rosettadock(self,rmsdf):
        fh=open(rmsdf)
        rmsdlist=[]
        k=-1
        fl=os.listdir('./')
        pdbfl=[item for item in fl if item.endswith('.pdb')]
        for line in fh:
            if k==-1:
                k=k+1
                continue
            ls=line.split()
            if k==0:
                if ls[0].startswith(self.dsname+'.'+self.code):
                    prefix=''
                elif ls[0].startswith(self.code):
                    prefix=self.dsname+'.'
                else:
                    prefix=self.dsname+'.'+self.code+'.'
            if not prefix+ls[0] in pdbfl:
                continue
            rmsdlist.append([prefix+ls[0].strip()[:-4],float(ls[2]),float(ls[5]),float(ls[1])])
            k=k+1
        self.extrapar=[('strmad','f4')]
        rmsdlist.sort(key=itemgetter(1))
        for item in rmsdlist:
            self.score.append(item[-1])
        self.dnlist=[item[0] for item in rmsdlist]
        self.build_sa(rmsdlist)

    def load_rmsd(self,rmsdf):
        with open(rmsdf,'rb') as fh:
            rmsddict=pickle.load(fh)
        rmsdlist=[]
        k=0
        dopescore=False
        if os.path.isfile('score.pickle'):
            dopescore=True
            with open('score.pickle','rb') as fh:
                dopedict=pickle.load(fh)
        if os.path.isfile('extrapar.pickle'):
            self.extrapar=pickle.load(open('extrapar.pickle', 'rb'))
        for key in rmsddict:
            if dopescore:
                if not (key in dopedict):
                    continue
                else:
                    self.score.append(dopedict[key])
            if k==0:
                if key.startswith(self.dsname+'.'+self.code+'.'):
                    prefix=''
                elif key.startswith(self.code+'.'):
                    prefix=self.dsname+'.'
                else:
                    prefix=self.dsname+'.'+self.code+'.'
            if isinstance(rmsddict[key],list):
                if not os.path.isfile('extrapar.pickle'):
                    raise Exception('no extrapar specified')
                rmsdlist.append([prefix+key]+rmsddict[key]+[key])
            else:
                rmsdlist.append([prefix+key,rmsddict[key],key])
        rmsdlist.sort(key=itemgetter(1))
        self.dnlist=[item[0] for item in rmsdlist]
        self.build_sa(rmsdlist)

    def calc_rmsd_fromfiles(self,nativepattern='native.*pdb$'):
        #calculate rmsd for files in a single directory#
        os.chdir(self.dir)
        rmsdlist=[]
        flist=os.listdir(self.dir)
        rmsdfile=open(self.path+'.pir','w')
        nativef=[elem for elem in flist if re.search(nativepattern,elem)][0]
        rmsdfile.write('cRMSD between '+nativef+' and '+nativef+' is   '+str(0.0)+'\n')
        rmsdlist.append([nativef[:-4],0.0])
        flist.remove(nativef)
        otherf=[elem for elem in flist if re.search('pdb$', elem)]
        for i in otherf:
            try:
                rmsdv=self.rmsd_2file(self.dir+'/'+nativef,self.dir+'/'+i)
                rmsdfile.write('RMSD between '+nativef+' and '+i+' is   '+str(rmsdv)+'\n')
                rmsdlist.append([i[:-4],rmsdv])
            except Exception as e:
                traceback.print_exc()
        rmsdfile.close()
        rmsdlist.sort(key=itemgetter(1))
        self.dnlist=[item[0] for item in rmsdlist]
        self.build_sa(rmsdlist)

    def build(self,filter=False):
        """
        Build the Decoys4Single object from the files in the directory, and copy the built decoys into SGE cluster.

        The directory should have the following files:
            Decoy structure files
                * structure pdbs
                * needtransformation: a base pdb file *base* and a transformation file *res*, see :meth:`Decoys4Single.build_dockingdecoys_withtransformation()`
                * needattachtobase: a base pdb file *base* and files with the remaining parts of decoy structures, the final decoy structure files are assembled on the fly by attaching the decoy file to the base file.
                * needcombinewithbase: a base pdb file *base* and files with the remaining parts of decoy structures. The base file should contains "replacethis" str, which will be replaced by the decoy file on the fly.
            .. note::
                Example: *base* can be the receptor file shared by all decoys/ligands/peptides
            rmsd.pickle | rosetta sc file |  a native structure \*native\*, see :meth:`Decoys4Single.calc_rmsd`
        """
        #pdb.set_trace()
        os.chdir(self.dir)
        if os.path.isfile("ma"):
            return self.build_dockingdecoys_ma()
        elif os.path.isfile(self.code+'.zd3.0.zr.out'):
            return self.build_dockingdecoys_zdock()
        elif os.path.isfile("needtransformation"):
            return self.build_dockingdecoys_withtransformation()
        else:
            self.rename()
            if filter:
                self.filter_pdb()
            self.calc_rmsd()
            if os.path.isfile('base.pir'):# or os.path.isfile('needattachtobase'):
                self.make_basepir()
                self.replicate_pir()
            else:
                ndnlist=pir(path=self.path+'.pir').make_pir_fromlist(self.dnlist)
                if len(ndnlist)!=len(self.dnlist):
                    nindex=[]
                    for i in range(len(self.dnlist)):
                        if self.dnlist[i] in ndnlist:
                            nindex.append(i)
                    self.dnlist=ndnlist
                    self.sa=self.sa[nindex]
            self.save()
            self.copy2ddb()

    def make_basepir(self,path=''):
        fl=os.listdir('./')
        bf=[f for f in fl if 'base' in f and (f.endswith('pdb') or f.endswith('.pdb.gz') or f.endswith('gz'))]
        if path:
            nativefile=path
        elif os.path.isfile("needtransformation"):
            pass
        elif os.path.isfile('basebase.pir'):
            return
        elif len(bf)>0 and os.path.isfile('needcombinewithbase'):
            basefile=bf[0]
            fh=gzip.open(basefile,'r')
            fc=fh.read()
            fh.close()
            fh=gzip.open(self.dsname+'.'+self.code+'.native.pdb.gz','r')
            dfc=fh.read()
            fh.close()
            fc=fc.replace('replacethis',dfc)
            with open('nativepdb','w') as fh:
                fh.write(fc)
            nativefile='nativepdb'
        elif len(bf)>0 and os.path.isfile('needattachtobase'):
            nativeligand=[f for f in fl if self.dnlist[0] in f and (f.endswith('pdb.gz') or f.endswith('pdb'))][0]
            if 'gz' in nativeligand:
                lf=gzip.open(nativeligand).read()
            else:
                lf=open(nativeligand).read()
            if 'gz' in bf[0]:
                open('temp','w').write(gzip.open(bf[0]).read()+'\n'+lf)
            else:
                open('temp','w').write(open(bf[0]).read()+'\n'+lf)
            nativefile='temp'
        env=Environ()
        mdl=Model(env,file=nativefile)
        aln=Alignment(env)
        aln.append_model(mdl,align_codes='base')
        aln[0].name='tbrname'
        aln.write('base.pir')

    def build_dockingdecoys_ma(self):
        print('building docking decoys')
        fl=os.listdir('./')
        #renaming the chains
        subprocess.check_call(['PatchDock/mainchain.pl',
                               self.code+'_r_u.pdb', 'A'])#use unbound
        subprocess.check_call(['PatchDock/mainchain.pl',
                               self.code+'_l_u.pdb', 'B'])
        #combine the two files
        with open(self.code+'_r_b.pdb','a') as fh:
            fh.write('\n')
        cat_files((self.code+'_r_u.pdb', self.code+'_l_u.pdb'),
                  open(self.dsname+'.'+self.code+'.base.pdb', 'w'))
        #generate pir file for the base structure
        self.make_basepir(self.dsname+'.'+self.code+'.base.pdb')
        basepir=open('base.pir').read()
        totalpir=basepir.replace('base',self.dsname+'.'+self.code+'.base')+'\n'
        with open(self.code+'.pir','w') as fh:
            fh.write(totalpir)
        #Read in the rmad, transformation matrix
        with open('ma') as fh:
            fc=fh.read()
        fcl=fc.split(',')
        va=[float(item) for item in fcl]
        rmsdlist=[[self.dsname+'.'+self.code+'.base',0]+va]
        self.dnlist=[item[0] for item in rmsdlist]
        self.extrapar=[('dg','f4'),('dasa','f4'),('irmsd','f4')]
        self.build_sa(rmsdlist)
        #copy to cluster
        tdir = os.path.join(runenv.serverUserPath, 'decoys', self.dsname,
                            self.code)
        subprocess.check_call(['ssh', runenv.jobserver, 'mkdir', '-p', tdir])
        subprocess.check_call(['scp', '*base*', runenv.jobserver+':'+tdir])
        subprocess.check_call(['scp', 'need*', runenv.jobserver+':'+tdir])
        self.save()

    def build_dockingdecoys_zdock(self):
        print('building docking decoys')
        fl=os.listdir('./')
        pdbfilelist=[item for item in fl if item.endswith('pdb')]
        #renaming the chains
        for g in glob.glob(os.path.join(runenv.ddbdir, 'ppd', 'benchmark4',
                                        self.code, '*')):
            shutil.copy(g, '.')
        subprocess.check_call(['PatchDock/mainchain.pl',
                               self.code.upper()+'_r_u.pdb', 'A'])
        subprocess.check_call(['PatchDock/mainchain.pl',
                               self.code.upper()+'_l_u.pdb', 'B'])
        #combine the two files
        with open(self.code.upper()+'_r_u.pdb','a') as fh:
            fh.write('\n')
        cat_files((self.code.upper()+'_r_u.pdb', self.code.upper()+'_l_u.pdb'),
                  open(self.dsname+'.'+self.code+'.base.pdb', 'w'))
        #generate pir file for the base structure
        open('needtransformation', 'w')
        self.make_basepir(self.dsname+'.'+self.code+'.base.pdb')
        #Read in the rmad, transformation matrix
        fh=open('rmsd.res')
        collect=False
        rmsdlist=[]
        transmatrixlist=[]
        score=[]
        for line in fh:
            #if re.search('# | rmsd',line):
            #    collect=True
            #    continue
            if True:#collect:
                try:
                    lines=line.split('|')
                    if len(lines)<2:
                        break
                    srmsd=float(lines[1].split('(')[0])
                    irmsd=float(lines[1].split('(')[1].split(')')[0])
                    rmsdlist.append([lines[-1][:-1],self.dsname+'.'+self.code+'.'+lines[0].strip(),srmsd,irmsd])
                except:
                    pass
                    #pdb.set_trace()
        k=0
        rfh=open(self.code+'.zd3.0.zr.out.rmsds')
        for line in rfh:
            rmsdlist[k].append(float(line.split()[1]))
            k=k+1
        k=0
        sfh=open(self.code+'.zd3.0.zr.out.scores.txt')
        for line in sfh:
            rmsdlist[k].append(float(line.split()[1]))
            k=k+1
        rmsdlist.sort(key=itemgetter(2))
        transmatrixlist=[item[0] for item in rmsdlist]
        rmsdlist=[item[1:] for item in rmsdlist]
        scorelist=[item[-1] for item in rmsdlist]
        rmsdlist=[item[:-1] for item in rmsdlist]
        self.dnlist=[item[0] for item in rmsdlist]
        self.extrapar=[('irmsd','f4'),('birmsd','f4')]
        self.build_sa(rmsdlist)
        self.sa['score']=np.array(scorelist)
        #replicate the pir file
        basepir=open('base.pir').read()
        bpl=basepir.split('-1.00\n')
        basepir=bpl[0]+'-1.00\nACS*\n'
        totalpir=''
        for k in range(0,len(self.dnlist)):
            totalpir=totalpir+basepir.replace('base',self.dnlist[k]).replace('tbrname',transmatrixlist[k])+'\n'
        with open(self.code+'.pir','w') as fh:
            fh.write(totalpir)
        #copy files to cluster
        tdir = os.path.join(runenv.serverUserPath, 'decoys', self.dsname,
                            self.code)
        subprocess.check_call(['ssh', runenv.jobserver, 'mkdir', '-p', tdir])
        subprocess.check_call(['scp', '*base*', runenv.jobserver+':'+tdir])
        subprocess.check_call(['scp', 'need*', runenv.jobserver+':'+tdir])
        self.save()

    def build_dockingdecoys_withtransformation(self):
        """Build decoys that need transformation, for pathdock only?"""
        with open('firedockinput.pickle', 'rb') as fh:
            sd=pickle.load(fh)
        print('building docking decoys')
        fl=os.listdir('./')
        #pf=[item for item in fl if re.search('_r_',item)][0]
        #lf=[item for item in fl if re.search('_l_',item)][0]
        resfile=[item for item in fl if item[-3:]=='res'][0]
        #renaming the chains
        #print(os.system('PatchDock/mainchain.pl '+sd['rpdb'][:-3]+' A'))
        #print(os.system('PatchDock/mainchain.pl '+sd['lpdb'][:-3]+' B'))
        #combine the two files
        cat_files((sd['rpdb'], sd['lpdb']),
                  open(self.dsname+'.'+self.code+'.base.pdb', 'w'))
        #generate pir file for the base structure
        self.make_basepir(self.dsname+'.'+self.code+'.base.pdb')
        #Read in the rmad, transformation matrix
        fh=open(resfile)
        collect=False
        rmsdlist=[]
        transmatrixlist=[]
        print("Getting list")
        for line in fh:
            if re.search('# | rmsd',line):
                collect=True
                continue
            if collect:
                try:
                    lines=line.split('|')
                    if len(lines)<2:
                        break
                    srmsd=float(lines[1].split('(')[0])
                    irmsd=float(lines[1].split('(')[1].split(')')[0])
                    rmsdlist.append([self.dsname+'.'+self.code+'.'+lines[0].strip(),srmsd,irmsd,lines[-1][:-1]])
                except:
                    pass
                    #pdb.set_trace()
        print("Finished getting list")
        rmsdlist.sort(key=itemgetter(1))
        print("sort finished")
        self.dnlist=[item[0] for item in rmsdlist]
        transmatrixlist=[item[-1] for item in rmsdlist]
        rmsdlist=[item[:-1] for item in rmsdlist]
        self.extrapar=[('irmsd','f4')]
        self.build_sa(rmsdlist)
        #replicate the pir file
        basepir=open('base.pir').read()
        bpl=basepir.split('-1.00\n')
        sl=int(len(''.join(bpl[1:]))/50+0.5)
        residuecodes=''.join(['A' for i in range(0,sl)])
        basepir=bpl[0]+'-1.00\n'+residuecodes+'*\n\n'
        totalpir=''
        print("generating pir")
        basepirlist=basepir.split('base')
        basepirlist=basepirlist[:2]+basepirlist[-1].split('tbrname')
        totalpir=''.join([basepirlist[0]+self.dnlist[k]+basepirlist[1]+self.dnlist[k]+basepirlist[2]+transmatrixlist[k]+basepirlist[3] for k in range(0,len(self.dnlist))])
        print("generating pir finished")
        with open(self.code+'.pir','w') as fh:
            fh.write(totalpir)
        tdir = os.path.join(runenv.serverUserPath, 'decoys', self.dsname,
                            self.code)
        subprocess.check_call(['ssh', runenv.jobserver, 'mkdir', '-p', tdir])
        subprocess.check_call(['scp', '*base*', runenv.jobserver+':'+tdir])
        subprocess.check_call(['scp', 'need*', runenv.jobserver+':'+tdir])
        subprocess.check_call(['scp', sd['lpdb'], runenv.jobserver+':'+tdir])
        subprocess.check_call(['scp', sd['rpdb'], runenv.jobserver+':'+tdir])
        subprocess.check_call(['scp', 'firedockinput.pickle',
                               runenv.jobserver+':'+tdir])
        self.save()

    def update_ppdt_pir(self):
        self.make_basepir(self.dsname+'.'+self.code+'.base.pdb')
        basepir=open('base.pir').read()
        bpl=basepir.split('-1.00\n')
        sl=int(len(''.join(bpl[1:]))/50+0.5)
        basepir=''.join(['A' for i in range(0,sl)])
        print(basepir)
        with open(self.code+'.pir','r') as fh:
            fc=fh.read()
        fc=fc.replace('ACS',basepir)
        with open(self.code+'.pir','w') as fh:
            fh.write(fc)
        return fc

    def get_sorted_rmsdlist(self):
        fh=open(runenv.decoysbasedir+'ddbrmsds/'+self.dsname+'/'+self.code+'.rmsds')
        codere=re.compile('and\s*(\S*)\s*is\s*([.0-9]*)')
        rmsddict={}
        rmsdlist=[]
        for f in fh:
            codeg=codere.search(f)
            if codeg.group(1)[:2]!=self.code[:2]:
                coden=self.code+'.'+codeg.group(1)
            else:
                coden=codeg.group(1)
            coden=coden[:-4]
            rmsddict[coden]=codeg.group(2)
            rmsdlist.append([coden,float(codeg.group(2))])
            rmsdlist.sort(key=itemgetter(1))
        self.rmsdlist=rmsdlist
        return rmsdlist

    def filter_rmsdlist(self,dnlist):
        newrmsdlist=[]
        for i in range(0,len(self.rmsdlist)):
            if dnlist.count(self.rmsdlist[i][0]):
                newrmsdlist.append(self.rmsdlist[i])
        self.rmsdlist=newrmsdlist

    def get_sorted_pir_str(self):
        print(0)

    def rmsd_2file(self,file1,file2,rmsdtype='all',looprmsd=False, nofit=True):
        e=runenv.env
        aln=Alignment(runenv.env)
        m1=Mymodel(e,file=file1)
        print('start calculating rmsd')
        print(file1)
        aln.append_model(m1, align_codes='c1', atom_files=file1)
        print(file2)
        m2=Mymodel(e,file=file2)
        aln.append_model(m2, align_codes='c2', atom_files=file2)
        if looprmsd:
            atmsel=m1.select_loop_atoms(self.loops)
        else:
            atmsel=Selection(m1)
        if rmsdtype=='all':
            pass
        elif rmsdtype=='mcrmsd':
            atmsel = atmsel.only_mainchain()
        elif rmsdtype=='scrmsd':
            atmsel = atmsel.only_sidechain()
        else:
            raise Exception('Unknown rmsdtype: '+str(rmsdtype))
        if not nofit:
            r = atmsel.superpose(m2, aln)
        else:
            r=atmsel.superpose(m2, aln,fit=False,refine_local=False)
        # We can now use the calculated RMS, DRMS, etc. from the returned 'r' object:
        rms = r.rms
        drms = r.drms
        print(rms)
        print('finished calculating rmsd')
        return rms

    def save(self):
        with open(self.path+'.pickle','wb') as fh:
            cPickle.dump(self,fh)

    def load(self):
        with open(self.path+'.pickle','rb') as fh:
            no=cPickle.load(fh)
        return no

    def replicate_pir(self):
        basepir=open('base.pir').read()
        totalpir=''
        for dn in self.dnlist:
            totalpir=totalpir+basepir.replace('base',dn)+'\n'
        with open(self.code+'.pir','w') as fh:
            fh.write(totalpir)

    def get(self):
        if os.path.isfile(self.path+'.pickle'):
            return self.load()
        else:
            self.build()
            return self

    def add_property_tosa(self,sa,dnlist,property='mcrmsd',ddbdir=''):
        os.chdir(ddbdir)
        if not (dnlist[0].endswith('.native.pdb') or dnlist[0].endswith('.native.pdb.gz') or dnlist[0].endswith('.native')):
            pdb.set_trace()
            raise Exception('The first structure in dnlist is not native structure')
        nativefile=dnlist[0]
        self.get_loop_region(nativefile)
        if len(self.loops)>0:
            looprmsd=True
        else:
            looprmsd=False
        e = Environ()
        e.io.atom_files_directory = ['./']
        for i in range(len(dnlist)):
            print(dnlist[i] + '  '+str(sa[i]['rmsd']))
            sa[i][property]=self.rmsd_2file(e,nativefile+'.pdb',dnlist[i]+'.pdb',rmsdtype=property,looprmsd=looprmsd)

    def get_loop_region(self,dname):
        if not dname.startswith(('ml','mgl')):
            self.loops=[]
            return 0
        ns=dname.split('.')
        lcode=ns[1][4:]
        lcode=lcode.replace('-',' ')
        rer=re.search('([0-9]{1,4})([a-zA-Z\s]{1,2})([0-9]{1,4})',lcode)
        self.loops=[[rer.group(2),int(rer.group(1)),int(rer.group(3)),int(rer.group(1))+int(rer.group(3))-1]]

class DecoySet(object):
    """
    A single set of decoys
    """
    def __init__(self,dsname='',sourcedir=''):
        self.dsname=dsname
        self.dir=os.path.join(runenv.decoysbasedir, dsname)
        self.path=self.dir+dsname
        self.pirpath=self.dir+dsname+'.pir'
        self.sourcedir=sourcedir

    def mkdir(self):
        if os.path.isdir(self.dir):
            return 0
        os.mkdir(self.dir)

    def clean_pir(self):
        pir(self.pirpath).get_pir_fromlist(self.dnlist,self.pirpath)

    def bulid_rmsd_array(self,sal):
        #build the rmsd score array based on all the known parameters
        self.dtype=sal[0].dtype+sal[0].extrapar
        itemlist=[]
        for item in sal[0].dtype+sal[0].extrapar:
            itemlist.append(item[0])
        self.sa=np.zeros(len(self.dnlist), dtype=sal[0].sa.dtype)#,('decoyn','a50')
        for i in range(0,len(self.pos)-1):
            self.sa['decoyseqposn'][self.pos[i]:self.pos[i+1]]=i
            self.sa['decoyposn'][self.pos[i]:self.pos[i+1]]=range(self.pos[i],self.pos[i+1])
            for item in itemlist[3:]:
                self.sa[item][self.pos[i]:self.pos[i+1]]=sal[i].sa[item]

    def get_codelist(self):
        #get the list of decoy sequences
        rmsdfiles=os.listdir(runenv.decoysbasedir+'ddbrmsds/'+self.dsname+'/')
        self.codelist=[item[:-6] for item in rmsdfiles]
        self.codelist.sort()

    def pre_build(self):
        os.chdir(self.sourcedir)
        folderlist=os.listdir('./')
        for folder in folderlist:
            if not os.path.isdir(self.sourcedir+folder):
                continue
            os.chdir(self.sourcedir+folder)
            subfolderlist=os.listdir('./')
            for subfolder in subfolderlist:
                #if subfolder[-10:]=='.coord.pdb':
                #    if os.path.isfile(subfolder[-10:]+'.pdb'):
                #        print(os.system('rm '+subfolder))
                #        continue
                #    print(os.system('python ../ergePDB.py '+subfolder[:-10]+' '+folder+' '+subfolder[:-10]))
                #    print(os.system('rm '+subfolder))
                #print(os.system('mv  ./'+subfolder+'/'+subfolder+'.pdb ./'+folder+'/'+subfolder+'/'+subfolder+'.native.pdb'))
                shutil.move(subfolder, '..')
            shutil.rmtree(os.path.join('..', folder))

    def build(self):
        """
        Build the decoy set.

        The pre-prepared decoys directory should look like the following::

          root/
              decoySet1/ could use the PDB code of the receptor here
                  Decoy files
                  rmsd file
              decoySet2/
              ...

        .. seealso::
            :meth:`Decoys4Single.build`
        """
        dir=self.sourcedir
        self.mkdir()
        codelist=[]
        os.chdir(dir)
        fl=os.listdir('./')
        fll=[]
        for f in fl:
            if os.path.isdir(f):
                fll.append(f)
        self.codelist=fll
        sdl=[]
        tdir = os.path.join(runenv.serverUserPath, 'decoys', self.dsname)
        subprocess.check_call(['ssh', runenv.jobserver, 'mkdir', '-p', tdir])
        sal=self.cat_singles()
        self.bulid_rmsd_array(sal)
        print("building array finished")
        self.save()

    def filter2good(self,bnors=2000,scorefile=''):
        if not scorefile:
            rsp=rawsp(pdbset='X_2.2A',features='a158a158d30'#,'a158a158d30','a158a158d200','a158a158d300']
                 ,genmethod='ss1')
            ssp=scaledsp(rspo=rsp,pm='npend')
            sps=spdsscore(ssp=ssp,dslist=[self.dsname],bm='ss1')
            sps.test=True
            sps.get_score()
            fh=open(runenv.basedir+'X_2.2A/a158a158d30-ss3/X_2.2A.a158a158d30.ss1.npend.'+self.dsname+'.ss1.pickle','rb')
        else:
            fh=open(scorefile)
        spscore=cPickle.load(fh)
        fh.close()
        ndnlist=set(spscore)
        odnlist=self.dnlist
        for i in range(0,len(odnlist)):
            if len(odnlist[i])>40:
                odnlist[i]=odnlist[i][0:40]
        dnlist=[]
        #rmsdlist=[]
        codelist=[]
        pos=[0]
        ind=[]
        dpn=0
        for i in range(0,len(self.codelist)):
            print("filtering "+self.codelist[i])
            k=0
            for j in range(self.pos[i],self.pos[i+1]):
                if (odnlist[j] in ndnlist) and odnlist[j][0:len(self.dsname)]==self.dsname:
                    #rmsdlist.append(self.rmsds[j])
                    dnlist.append(self.dnlist[j])
                    k=k+1
                    ind.append(j)
            if k>0:
                codelist.append(self.codelist[i])
            else:
                continue
            self.sa['decoyseqposn'][self.pos[i]:self.pos[i+1]]=dpn
            pos.append(pos[-1]+k)
            dpn=dpn+1
        if len(dnlist)<0.8*len(odnlist):
            print('Error: too many decoys (>20%) are filtered out using the filter, please check code for bugs')
            sys.exit('Error: too many decoys (>20%) are filtered out using the filter, please check code for bugs')
        self.pos=pos
        #self.rmsds=rmsdlist
        self.dnlist=dnlist
        self.codelist=codelist
        self.sa=self.sa[ind]
        self.sa['decoyposn']=np.arange(0,self.sa.shape[0])
        self.clean_pir()
        self.save()

    def cat_singles(self):
        sdl=[]
        pos=[0]
        dnlist=[]
        sal=[]
        pirlist=[]
        print('##############################################building decoys '+self.dsname)
        #runtask=task()
        #runtask.parallel_local(self.build_single,inputlist=self.codelist,nors=10)
        #pdb.set_trace()
        for code in self.codelist:
            print("***********building "+code)
            self.build_single(code)
            sd=Decoys4Single(code,self.dsname,
                             os.path.join(self.sourcedir, code))
            sd.get()
            with open(os.path.join(self.sourcedir, code, code+'.pickle'),'rb') as fh:
                sd=cPickle.load(fh)
            dnlist=dnlist+sd.dnlist
            sal.append(sd)
            pos.append(pos[-1]+len(sd.dnlist))
            if len(sd.dnlist)==1:
                raise Bugs('Only 1 decoys in the decoy set, probably something wrong with it')
            pirlist.append(sd.path+'.pir')
        self.pos=pos
        self.dnlist=dnlist
        cat_files(pirlist, open(self.pirpath, 'w'))
        return sal

    def build_single(self,code):
        sd=Decoys4Single(code,self.dsname, os.path.join(self.sourcedir, code))
        print('generating decoy for '+code)
        sd=sd.get()
        print('generating finished '+code)
        return sd

    def save(self):
        os.chdir(self.dir)
        try:
            mypickle().dump(self,self.dsname)
        except Exception as e:
            traceback.print_exc()
            with open(self.dsname, 'rb') as fh:
                pickle.dump(fh)
            pdb.set_trace()
            mypickle().dump(self,self.dsname)

    def load(self):
        os.chdir(self.dir)
        return mypickle().load(self.dsname)

    def get(self,rebuild=False):
        if os.path.isfile(self.path+'.pickle') and not rebuild:
            return self.load()
        else:
            self.build()
            return self

    def check_decoys(self):
        for pos in self.pos:
            code1=self.dnlist[pos]
            code2=self.dnlist[pos+1]
            fh1=open(runenv.ddbdir+code1+'.pdb')
            fh2=open(runenv.ddbdir+code2+'.pdb')
            fh1c=fh1.read()
            fh2c=fh2.read()
            l1=len(re.findall('\\nATOM',fh1c))
            l2=len(re.findall('\\nATOM',fh2c))
            if l1!=l2:
                print(code1+' '+str(l1)+' '+str(l2))

    def update_ppdt_pir(self):
        pirlist=[]
        ss=self.load()
        for code in ss.codelist:
            os.chdir(self.sourcedir+code)
            sd=Decoys4Single(code,self.dsname,self.sourcedir+code+'/')
            sd.update_ppdt_pir()
            pirlist.append(sd.path+'.pir')
        cat_files(pirlist, open(self.pirpath, 'w'))

    def get_loopdict(self,dp):
        ld={}
        for code in a.codelist:
            rer=re.search('([0-9]+)([A-Z\-]{1})([0-9]+)',code[4:])
            cl=int(rer.group(1))
            cs=int(rer.group(3))
            cc=rer.group(2)
            if cc=='-':
                cc=''
            if code[0:4] in ld:
                ld[code[0:4]].append([cc,cl,cs,cs+cl-1])
            else:
                ld[code[0:4]]=[[cc,cl,cs,cs+cl-1]]
            with open(dp,'wb') as fh:
                pickle.dump(ld,fh)

    def build_index_list(self):
        self.indexlist=[]
        for i in range(0,len(self.pos)-1):
            self.indexlist.append([self.pos[i],self.pos[i+1]])

    def get_codepos_fordifferentset(self):
        sposindex=1
        codepos=[0]
        for i in range(1,len(self.indexlist)):
            if self.indexlist[i][1]>self.spos[sposindex] and self.indexlist[i-1][1]<=self.spos[sposindex]:
                codepos.append(i)
                sposindex+=1
        codepos.append(len(self.indexlist))
        self.setpos=codepos

    def rename_codes(self):
        if not 'dslist' in self.__dict__.keys():
            for i in range(len(self.codelist)):
                self.codelist[i]=self.dsname+self.codelist[i]
        else:
            for i in range(0,len(self.setpos)-1):
                for j in range(self.setpos[i],self.setpos[i+1]):
                    self.codelist[j]=self.dslist[i]+self.codelist[j]

    def add_property_tosa(self,sa, property='mcrmsd',ddbdir=''):
        self=self.load()
        self.rebuild_codelist()
        if re.search(property, str(self.sa.dtype)):
            sa[...]=self.sa
            return 0
        self.sa=sa
        for i in range(len(self.codelist)):
            try:
                singlecode=self.dnlist[self.pos[i]].split('.')[1]
                Decoys4Single(singlecode).add_property_tosa(self.sa[self.pos[i]:self.pos[i+1]],self.dnlist[self.pos[i]:self.pos[i+1]],property,ddbdir=ddbdir+singlecode+'/')
            except:
                pdb.set_trace()
        self.save()

    def rebuild_codelist(self):
        codelist=[]
        for i in range(len(self.pos)-1):
            codelist.append(self.dnlist[self.pos[i]].split('.')[1])
        self.codelist=codelist

class DecoySets(DecoySet):
    def __init__(self,dslist=[],sourcedir=''):
        self.dslist=dslist
        DecoySet.__init__(self,'_'.join(dslist))
        self.sourcedir=sourcedir

    def get(self,rebuild=False):
        self.mkdir()
        dss=[]
        self.dnlist=[]
        self.rmsds=[]
        self.pos=[]
        self.codelist=[]
        self.spos=[0]
        for dsname in self.dslist:
            ds=DecoySet(dsname,sourcedir=self.sourcedir+dsname+'/')
            ds=ds.get(rebuild)
            dss.append(ds)
            self.dnlist=self.dnlist+ds.dnlist
            #self.rmsds=self.rmsds+ds.rmsds
            self.codelist=self.codelist+ds.codelist
            self.pos=self.pos[:-1]+[item+self.spos[-1] for item in ds.pos]
            self.spos.append(self.spos[-1]+len(ds.dnlist))
        self.cat_rmsd_array(dss)
        self.gen_pir()
        return self

    def build(self,rebuild=False):
        self.get(rebuild)
        self.save()

    def cat_rmsd_array(self,dss):
        #build the rmsd score array based on all the known parameters
        itemlist=[]
        self.dtype=dss[0].sa.dtype.descr
        for item in self.dtype:
            itemlist.append(item[0])
        self.sa=np.zeros(len(self.dnlist), dtype=dss[0].sa.dtype)#,('decoyn','a50')
        for i in range(0,len(self.spos)-1):
            self.sa['decoyset'][self.spos[i]:self.spos[i+1]]=i
            self.sa['decoyseqposn'][self.spos[i]:self.spos[i+1]]=dss[i].sa['decoyseqposn']
            self.sa['decoyposn'][self.spos[i]:self.spos[i+1]]=range(self.spos[i],self.spos[i+1])
            for item in itemlist[3:]:
                self.sa[item][self.spos[i]:self.spos[i+1]]=dss[i].sa[item]

    def filter2good(self):
        for dsname in self.dslist:
            DecoySet(dsname=dsnmae).filter2good()
        self.build()

    def load_oldlist(self):
        #use this function to load the rmsdlist generated by old codes
        fh=open(runenv.basedir+'_'.join(['DecoysRUS','HR','casp58'])+'/'+'_'.join(['DecoysRUS','HR','casp58'])+'.rmsdlist')
        return cPickle.load(fh)

    def calc_nativelist(self):
        rmsdlist=self.load_oldlist()
        sa=rmsdlist[0]
        #sas=np.sort(sa,order=['decoyset','decoyseqposn','rmsd'])
        dnlist=rmsdlist[4]
        nativelist=[dnlist[sa['decoyposn'][i]] for i in range(0, len(dnlist)) if (sa['rmsd'][i]==0 and sa['decoyset'][i]==0)]
        self.gen_pir()
    #    pirfile=pir(self.pirpath)
    #    pirfile.get_pir_fromlist(nativelist,self.pirpath[:-3]+'native.pir')
        with open(self.dir+'nativelist','w') as nlfh:
            nlfh.write(','.join(nativelist))
        return nativelist

    def load_nativelist(self):
        nlfh=open(self.dir+'nativelist')
        nlc=nlfh.read()
        return nlc.split(',')

    def get_nativelist(self):
        if os.path.isfile(self.dir+'nativelist'):
            return self.load_nativelist()
        else:
            return self.calc_nativelist()

    def gen_pir(self):
        for i in range(0, len(self.dslist)):
            shutil.copy(DecoySet(self.dslist[i]).pirpath,
                        self.dir+'set'+str(i)+'.pir')
        cat_files(glob.glob(os.path.join(self.dir, 'set*pir')),
                  open(self.pirpath, 'w'))
        for g in glob.glob(os.path.join(self.dir, 'set*pir')):
            os.unlink(g)

    def get_ds_part(self,sample):
        sset=set(sample)
        codelist=[]
        codelist2=[]
        indexlist=[]
        indexlist2=[]
        for i in range(0,len(self.codelist)):
            code=self.codelist[i]
            if not code in sset:
                codelist2.append(code)
                indexlist2.append(self.indexlist[i])
            else:
                codelist.append(code)
                indexlist.append(self.indexlist[i])
        nds=copy.copy(self)
        nds.indexlist=indexlist
        nds.codelist=codelist
        nds.get_codepos_fordifferentset()
        nds2=copy.copy(self)
        nds2.indexlist=indexlist2
        nds2.codelist=codelist2
        nds2.get_codepos_fordifferentset()
        return [nds,nds2]

    def get_ds_list(self):
        ndslist=[copy.copy(self) for i in range(len(self.codelist))]
        for i in range(0,len(self.codelist)):
            ndslist[i].codelist=[self.codelist[i]]
            ndslist[i].indexlist=[self.indexlist[i]]
            ndslist[i].get_codepos_fordifferentset()
        return ndslist

    def add_property_tosa(self,property='mcrmsd',ddbdir=''):
        self=self.load()
        self.rebuild_codelist()
        from numpy.lib.recfunctions import append_fields
        self.sa=append_fields(self.sa,data=np.zeros(len(self.sa)),names=property,dtypes='f4',usemask=False)
        for i in range(len(self.dslist)):
            DecoySet(dsname=self.dslist[i]).add_property_tosa(self.sa[self.spos[i]:self.spos[i+1]],property,ddbdir=ddbdir+self.dslist[i]+'/')
        self.save()

    def load(self,mode='dss'):
        if mode=='dss':
            os.chdir(self.dir)
            if os.path.isfile('dss.pickle'):
                return mypickle().load('dss')
            else:
                if len(self.dslist)==1:
                    ndss=DecoySet.load(self)
                else:
                    ndss=self.get()
                ndss.rebuild_codelist()
                if len(self.dslist)==1:
                    ndss.spos=[0,ndss.pos[-1]]
                ndss.build_index_list()
                ndss.get_codepos_fordifferentset()
                ndss.rename_codes()
                nds=DecoySets(self.dslist)
                del ndss.dnlist
                nds.pos=copy.copy(ndss.pos)
                #self.ds.spos=copy.copy(nds.spos)
                nds.codelist=copy.copy(ndss.codelist)
                nds.sa=np.copy(ndss.sa)
                nds.spos=copy.copy(ndss.spos)
                nds.indexlist=copy.copy(ndss.indexlist)
                nds.setpos=copy.copy(ndss.setpos)
                del ndss
                os.chdir(self.dir)
                mypickle().dump(nds,'dss')
                return nds
        else:
            if len(self.dslist)==1:
                ndss=DecoySet.load(self)
            else:
                ndss=self.get()
            ndss.rebuild_codelist()
            if len(self.dslist)==1:
                ndss.spos=[0,ndss.pos[-1]]
            ndss.build_index_list()
            ndss.get_codepos_fordifferentset()
            ndss.rename_codes()
            return ndss

    def get_rmsd_stats(self,num=1,rmsdtype='rmsd'):
        self.rmsdarray=np.zeros([len(self.indexlist),1])
        poslist=self.indexlist
        for i in range(0,len(poslist)):
            if self.withnative:
                rmsds=self.sa[poslist[i][0]+1:poslist[i][1]]
            else:
                rmsds=self.sa[poslist[i][0]:poslist[i][1]]
            self.rmsdarray[i]=rmsds[rmsdtype][0:num].mean()
        rmsdl={}
        for i in range(len(self.setpos)-1):
            nod.append(self.setpos[i+1]-self.setpos[i])
            srmsd=self.rmsdarray[self.ds.setpos[i]:self.ds.setpos[i+1]]
            rmsdl[self.ds.dslist[i]]=[srmsd.mean(), srmsd.std(),srmsd.min(),srmsd.max()]
        print(rmsdl)
        return rmsdl

    def filter(self,fa):
        #self.codelist=self.codelist #
        #self.pos #pos in sa
        #self.spos #setpos in sa
        #self.indexlist=copy.copy(nds.indexlist)
        #self.setpos #setpos in codelist #
        #pdb.set_trace()
        #if len(self.indexlist) != len(self.pos)-1:
        #    raise Exception('DecoySets.filter() can only be run on full decoy set, indexlist must cover all decoys')
        self.sa=self.sa[fa]
        pos=[0]
        spos=[0]
        indexlist=[]
        codelist=[]
        setpos=[0]
        dslist=[]
        sk=-1
        k=-1
        for sp in range(len(self.setpos)-1):
            sk=sk+1
            for index in self.indexlist[self.setpos[sp]:self.setpos[sp+1]]:
                k=k+1
                if fa[index[0]:index[1]].sum()>0:
                    pos.append(pos[-1]+fa[index[0]:index[1]].sum())
                    indexlist.append([pos[-2],pos[-1]])
                    codelist.append(self.codelist[k])
            if not spos[-1]==pos[-1]:
                spos.append(pos[-1])
                dslist.append(self.dslist[sk])
                setpos.append(len(codelist))
        self.pos=pos
        self.spos=spos
        self.indexlist=indexlist
        if self.indexlist[-1][-1]!=self.sa.shape[0]:
            raise Exception('Error in filtering...')
        self.setpos=setpos
        self.codelist=codelist
        self.dslist=dslist
