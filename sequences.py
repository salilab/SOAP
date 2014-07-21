"""
   A module for processing protein sequences in PIR format.
"""

from env import *

class pir(object):
    """
    Process PIR files for use by MDT.

    MDT calculates statistics from structures by looping through
    a PIR file containing all sequences for the structures.

    The last two fields in the PIR header (X-ray resolution and the r factor)
    can be repurposed to define the error bar parameters associated
    with each structure, used by MDT to calculate the error bar
    on the position of each atom.

    :Parameters:
      - `path`: path to the PIR file.
    """
    def __init__(self, path=''):
        self.pirpath=path

    def get_pir(self,pdbset):#X-xray, U-biological unit, T-transmembrane protein, M-membrane protein,S-single chain, W-non-membrane protein,1-single chain pir or none, 2-multiple chain pir, H-helical protein, B-beta protein
        import pp
        pirfiledirectory=runenv.pdbpirdir
        pdbfile='pdb_'+pdbset+'.pir'
        listoffiles=os.listdir(pirfiledirectory)
        if pdbfile in listoffiles:
            return pdbfile
        pirlist=pdbset.split('_')
        #if pirlist[-1].startswith('')
        if re.search('_[0-9]{2}',pdbset[-3:]):
            self.property_filter(pdbset[:-3])
            os.chdir(pirfiledirectory)
            pp.pir_sif('pdb_'+pdbset[:-3]+'.pir', pdbfile,si=int(pdbset[-2:]))
        else:
            self.property_filter(pdbset)
        return pdbfile

    def get_property_dict(self,pdbset):
        #pdb.set_trace()
        fc=pdbset.split('_')
        propertydict={}
        k=-1
        for fi in fc:
            k=k+1
            if k==0:
                propertydict['type']=fi #X-xray, S-single chain, M-Membrane, W-non-membrane, T-transmembrane, U-biological units
            elif fi.find('ph')!=-1:
                fic=fi.split('ph')
                propertydict['ph']=[float(fic[0]),  float(fic[1])]
                continue
            elif fi[-1]=='A':
                propertydict['resolution']= float(fi[:-1])
                continue
            elif fi[-1]=='r':
                propertydict['r']= float(fi[:-1])
                continue
            elif fi[-5:]=='rfree':
                propertydict['rfree']= float(fi[:-5])
                continue
            elif fi[-1]=='p':
                propertydict['similar']=[fi[:-3],float(fi[-3:-1])]
                continue
            else:
                propertydict[fi]=fi
                continue
        return propertydict

    def property_filter(self,pdbset):
        pirfiledirectory=runenv.pdbpirdir
        pdbfile='pdb_'+pdbset+'.pir'
        listoffiles=os.listdir(pirfiledirectory)
        if pdbfile in listoffiles:
            return pdbfile
        pirstring=self.get_pir_from_pdbdict(pdbset)
        fh=open(runenv.pdbpirdir+'pdb_'+pdbset+'.pir','w')
        fh.write(pirstring)
        fh.close()
        pir(path=runenv.pdbpirdir+'pdb_'+pdbset+'.pir').duplicate_check()

    def get_pir_from_pdbdict(self,pdbset):
        import pp
        if pdbset.startswith('THAC'):
            ep=pp.tmhpdb()
            ppdict=self.get_property_dict(pdbset[4:])
        elif pdbset.startswith('THMC'):
            ep=pp.tmhmcpdb()
            ppdict=self.get_property_dict(pdbset[4:])
        elif pdbset.startswith('BIOU'):
            ep=pp.bioupdb()
            ppdict=self.get_property_dict(pdbset[4:])
        else:
            ppdict=self.get_property_dict(pdbset)
            ep=pp.entirepdb()
        return ep.get_pir(ppdict)

    def make_pir_single_forautomodel(self):
        code=self.pirpath
        e=environ()
        e.io.atom_files_directory=runenv.pdbdir
        aln=alignment(e)
        mdl=model(e)
        mdl.read(file=code)
        aln.append_model(mdl,align_codes=code,atom_files=code)
        aln.append_model(mdl,align_codes=code+'decoys')
        aln.write(file= code+'.pir',alignment_format='PIR')

    def make_pir_fromlist(self,filelist):
        code=self.pirpath
        e=environ()
        e.io.atom_files_directory=['./']
        output = modfile.File(self.pirpath, 'w')
        ndnlist=[]
        fl=os.listdir('./')
        bf=[f for f in fl if 'base' in f and (f.endswith('pdb') or f.endswith('.pdb.gz') or f.endswith('gz'))]

        for f in filelist:
            try:
            #if 1:
                mdl=model(e)
                if os.path.isfile('needattachtobase'):
                    if os.path.isfile(f+'.pdb.gz'):
                        lf=gzip.open(f+'.pdb.gz').read()
                    else:
                        lf=open(f+'.pdb').read()
                    if 'gz' in bf[0]:
                        open('temp.pdb','w').write(gzip.open(bf[0]).read()+'\n'+lf)
                    else:
                        open('temp.pdb','w').write(open(bf[0]).read()+'\n'+lf)
                    mdl.read('temp.pdb')
                else:
                    mdl.read(f)
                aln=alignment(e)
                aln.append_model(mdl,align_codes=f,atom_files=f)
                aln.write(output,alignment_format='PIR')
                ndnlist.append(f)
            except Exception, error:
                print os.listdir('./')
                print error
                print f
                traceback.print_exc()
        return ndnlist

    def count_pir(self):
        input=file(self.pirpath)
        n=0;
        pirhead=re.compile(">P1;(.*)$")
        for line in input:
            if pirhead.match(line):
                n=n+1
        return n

    def build_dict(self):
        pirhead=re.compile(">P1;(.+)$")
        fh=open(self.pirpath,'r')
        pirdict={}
        code=''
        codepir=''
        for line in fh:
            if pirhead.match(line):
                pirdict[code]=codepir
                code=pirhead.match(line).group(1)
                codepir=''
            codepir=codepir+line
        pirdict[code]=codepir
        self.pirdict=pirdict
        fh.close()

    def get_codelist(self):
        input=file(self.pirpath)
        n=0;
        pirhead=re.compile(">P1;(.*)$")
        codelist=[]
        for line in input:
            if pirhead.match(line):
                rer=pirhead.match(line)
                codelist.append(rer.group(1))
        return codelist

    def count_residue2(self,residuepower=2):
        input=file(self.pirpath)
        n=0;
        pirhead=re.compile(">P1;(.*)$")
        linematch=re.compile("\S+.*")
        sly=0
        county=0
        for line in input:
            if county:
                if linematch.match(line):
                    numofr=numofr+len(line)-1
                else:
                    n=n+numofr**residuepower
                    county=0
            if sly==1:
                sly=0
                county=1
            if pirhead.match(line):
                sly=1;
                numofr=0;
                #print line
        return n

    def sep_pir(self, tardir, numoff,permin=100,residuepower=2):
        pirfile=self.pirpath
        bdir=tardir # input
        totalcodenum=self.count_pir()
        if totalcodenum<permin:
            print os.system('cp '+pirfile+' '+tardir+'/pdbs1')
            return 0
        if float(totalcodenum)/numoff<permin:
            numoff=int(float(totalcodenum)/permin)
        input=file(pirfile)
        tn=self.count_residue2(residuepower)
        maxn=tn/numoff;#input
        k=1;#input
        tc=0;
        n=0;
        os.chdir(tardir)
        out=file('pdbs'+str(k),"w")
        pirhead=re.compile(">P1;(.*)$")
        linematch=re.compile("\S+.*")
        sly=0
        county=0
        for line in input:
            if county:
                if linematch.match(line):
                    numofr=numofr+len(line)-1
                else:
                    n=n+numofr**residuepower
                    county=0
            if sly==1:
                sly=0
                county=1
            if pirhead.match(line):
                if n > maxn:
                    n=0
                    k=k+1
                    out.flush()
                    out.close()
                    out=file('pdbs'+str(k),"w")
                    #print 'seperating'+ str(k)
                sly=1
                numofr=0;
                tc=tc+1
            out.write(line)
        out.flush()
        out.close()

    def gen_pir(self, pirfiledirectory,pdbset):
        print 'Generating pdb file'
        env = environ()
        aln = alignment(env)
        pdbfile='pdb_'+pdbset+'.pir'
        listoffiles=os.listdir(pirfiledirectory)
        for filen in listoffiles:
            if pdbfile==filen:
                return pdbfile
        for filen in listoffiles:
            if re.match(filen[:-4],pdbfile):
                if re.match('[0-9.]*A',pdbfile[len(filen)-3:-4]):
                    print filen
                    input = modfile.File(pirfiledirectory+filen, 'r')
                    output = modfile.File(pirfiledirectory+pdbfile, 'w')
                    while aln.read_one(input, alignment_format='PIR'):
                        print "Read code %s" % aln[0].code
                        if aln[0].resolution < float(pdbfile[len(filen)-3:-5]):
                            aln.write(output, alignment_format='PIR')
                    input.close()
                    output.close()
                    return pdbfile
        print 'The program do not know how to generate the pdbpir file. Please generate it manually.'
        return 0

    def filter_pir(self,inputpirfile,outputpirfile,structuretype='structureX',structureresolution=99.0):
        if (os.access(outputpirfile,os.F_OK)):
            return 0
        env = environ()
        aln = alignment(env)
        input = modfile.File(inputpirfile, 'r')
        output = modfile.File(outputpirfile, 'w')
        while aln.read_one(input, alignment_format='PIR'):
            if aln[0].resolution < structureresolution and aln[0].prottyp==structuretype:
                aln.write(output, alignment_format='PIR')
        input.close()
        output.close()

    def sif_pir(self,inputpirfile,outputpirfile,si=30):
        if (os.access(outputpirfile,os.F_OK)):
            return 0
        #log.verbose()
        env = environ()
        sdb = sequence_db(env, seq_database_file=inputpirfile,
                      seq_database_format='PIR',
                      chains_list='ALL', minmax_db_seq_len=(30, 3000),
                      clean_sequences=True)
        sdb.filter(rr_file='${LIB}/blosum62.sim.mat', gap_penalties_1d=(-500, -50),
               matrix_offset = -450, max_diff_res=30, seqid_cut=si,
               output_grp_file=outputpirfile+'.grp', output_cod_file=outputpirfile+'.cod')
        out = file(outputpirfile, "w")
        codes = [line.rstrip('\r\n') for line in file(outputpirfile+'.cod')]
        codes = dict.fromkeys(codes)
        pirhead = re.compile(">P1;(.*)$")
        printline = False
        for line in file(inputpirfile):
            m = pirhead.match(line)
            if m:
                printline = m.group(1) in codes
            if printline:
                out.write(line)

    def gen_pir_new(self,pirfiledirectory,pdbset):
        print 'Generating pdb file ' + pdbset
        env = environ()
        aln = alignment(env)
        pdbfile='pdb_'+pdbset+'.pir'
        listoffiles=os.listdir(pirfiledirectory)
        for filen in listoffiles:
            if pdbfile==filen:
                return pdbfile
        if pdbset[0]=='X':
            sttype='structureX'
        elif pdbset[0]=='N':
            sttype='structureN'
        else:
            print 'unsupported pdb database file type (example X_2.2A_60): '+pdbset
        if re.search('[0-9.]*A',pdbset):
            stresolution=float(pdbset[2:5])
            print stresolution
        else:
            stresolution=99.0
        if not re.search('_[0-9]*$',pdbset):
            pdb1name=pdbfile
        else:
            pdb1name=pdbfile[-7:]+'.pir'
        filter_pir(pirfiledirectory+'pdball.pir',pirfiledirectory+pdbfile[0:10]+'.pir', structuretype=sttype,structureresolution=stresolution)
        if re.search('_[0-9]*$',pdbset):
            sic=int(pdbset[-2:])
            sif_pir(pirfiledirectory+pdbfile[0:10]+'.pir', pdbfile,si=sic)
        return pdbfile

    def get_pir_fromlist(self,list,targetpath):
        pirstring=''
        self.build_dict()
        k=0
        for fn in list:
            k=k+1
            print str(k)+' outof '+str(len(list))
            pirstring=pirstring+self.pirdict[fn]
        out=open(targetpath,'w')
        out.write(pirstring)
        out.close()

    def get_first_n(self,n,targetpath):
        pirstring=''
        self.build_dict()
        k=0
        for fn in self.pirdict.keys():
            k=k+1
            pirstring=pirstring+self.pirdict[fn]
            if k>=n:
                break
        out=open(targetpath,'w')
        out.write(pirstring)
        out.close()

    def filter_pir(self,xray=2,rfactor=0.25,phrange=[6.5,7.5],sid=60):
        tp=self.pirpath[:-4]+'_'+str(xray)+'_'+str(rfactor)+'_'+str(phrange[0])+'-'+str(phrange[1])+'.pir'
        if os.path.isfile(tp[:-4]+'_'+str(sid)+'.pir'):
            return 0
        self.build_dict()
        pirstring=''
        env=environ()
        env.io.atom_files_directory=[runenv.opdbdir]
        for key in self.pirdict:
            if len(key)<4:
                continue
            print key
            mdl=model(env)
            mdl.read(file=key[0:4])
            if mdl.rfactor>rfactor or mdl.resolution>xray or mdl.rfactor==-1 or mdl.resolution==-1:
                continue
            ph=get_ph(key[0:4])
            if ph<phrange[0] or ph>phrange[1]:
                print 'ph not in range '+str(ph)
                continue
            pirstring=pirstring+self.pirdict[key]
        out=open(tp,'w')
        out.write(pirstring)
        out.close
        import pp
        pp.pir_sif(tp, tp[:-4]+'_'+str(sid)+'.pir',si=sid)

    def clean_pdb(self):
        env=environ()
        os.chdir(runenv.basedir+'pdb/')
        env.io.atom_files_directory =['/salilab/park2/database/pdb/divided/']
        codelist=self.get_codelist()
        for code in codelist:
            mdl=model(env)
            mdl.read(code[0:4])
            mdl.write('pdb'+code[0:4]+'.ent')
            print os.system('scp pdb'+code[0:4]+'.ent '+runenv.jobserver+':~/pdb/')

    def duplicate_check(self):
        codelist=self.get_codelist()
        if len(set(codelist))<len(codelist):
            raise Exception('Duplicate entries in pir file possbily due to bugs in codes to generate the pir file: '+self.pirpath)

    def exclude(self,excludepirfile,resultpir):
        po=pir(path=excludepirfile)
        codelist=po.get_codelist()
        codelist=[code.split('.')[1][:4].lower() for code in codelist]
        codeset=set(codelist)
        self.build_dict()
        pirlist=[]
        for code in self.pirdict:
            if not code[:4].lower() in codeset:
                pirlist.append(self.pirdict[code])
        open(resultpir,'w').write('\n'.join(pirlist))
