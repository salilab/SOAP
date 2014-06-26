#functions to process pdb lists

from modeller import *
import re
import os
import time
import urllib
from modeller.automodel import *
from modeller.automodel.autosched import *
# Load in optimizer and schedule support
from modeller import optimizers
from modeller.schedule import schedule, step
from modeller.optimizers import actions
import sys
import pdb
import numpy as np
import cPickle
import pickle
from operator import itemgetter, attrgetter
import copy
import SOAP
import traceback
from SOAP import *

basedir='/bell2/gqdong/statpot/'
pfd=basedir+'pdbpir/'
resolre=re.compile("\s{2}RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*([0-9\.]+)")
completenessre=re.compile("\s{2}COMPLETENESS.+\(%\)\s*:\s*([0-9\.]+)")
rfre=re.compile("R VALUE.{1,40}:\s*(0\.\d+)")
rfactorre=re.compile("FREE R VALUE\s*\(WORKING SET\)\s*:\s*(0\.\d+)")
frfactorre=re.compile("\s{2}FREE R VALUE\s{14}.{1,20}:\s*(0\.\d+)")
norre=re.compile("NUMBER OF REFLECTIONS.{1,30}:\s*([0-9\.]+)")
atomre=re.compile("\\nATOM.{52}1.00\s*|\\nHETATM.{50}1.00\s*")
bre=re.compile("\\nATOM.{57}\s?([0-9\.]+)\s*")
abfactorre=re.compile("\s{2}MEAN B VALUE\s*\(OVERALL, A\*\*2\)\s*:\s*([0-9\.]+)")
originaldir='/bell2/gqdong/pdb/pdb1/'
flippeddir='/bell2/gqdong/pdb/pdb2/'
finaldir='/bell2/gqdong/pdb/pdb3/'
decoysdir={'ppd4s':'/bell2/gqdong/rawdecoyfiles/ppd/ppd4s/'}
debug=False

class entirepdb(object):
    #key:code
        #['code', 'medianCbeta', 'aveB', 'MolProbityScore', 'numRama', 'nucacids', 'pct_badangles', 'pir',
        #'completeness', 'maxCbeta', 'pct_badbonds', 'aveerr', 'numRota', 'meanCbeta', 'n_samples', 'ramaOutlier',
        #'maxresol', 'ph', 'type', 'Mol_pct_rank', 'errorscale', 'pct_rank40', 'residues', 'clashscoreB<40',
        #'clashscore', 'rota<1%', 'r', 'ramaAllowed', 'chains', 'pct_rank', 'minresol', 'numCbeta', 'ramaFavored',
        #'rfree', 'atomnum', 'numofreflections', 'resolution', 'cbeta>0.25']
            #chains:chainid
                #['missingresnum', 'length', 'missingatomnum', 'atomnum', 'type', 'pir']
    def __init__(self,originalpdbdir='',update=False):
        self.originalpdbdir=originalpdbdir
        self.xslist=[]
        self.update=update
        if update:
            fh=open('xslist.pickle')
            self.oxslist=cPickle.load(fh)
            fh.close()
        else:
            self.oxslist=[]
            self.nxslist=self.xslist
        self.errfh=open('/bell2/gqdong/pdb/errors','a')
        self.minchainlen=30
        self.pdbdict={}
        self.pdbdictname='pdbdict.pickle'

    def get_xray_structure_list(self,path=originaldir):
        os.chdir(path)
        print os.system('scp gqdong@baton2:/netapp/database/pdb/remediated/pdball.pir '+basedir+'pdbpir/pdball.pir')
        fh=open(basedir+'pdbpir/pdball.pir')
        fc=fh.read()
        pl=fc.split('>P1;')
        pl=pl[1:]
        for p in pl:
            code=p[0:4]
            psl=p.split('\n')
            if psl[1].startswith('structureX'):
                self.xslist.append(code)
        self.xslist=list(set(self.xslist))
        if self.update:
            fh=open('xslist.pickle')
            self.oxslist=cPickle.load(fh)
            fh.close()
            self.nxslist=[item for item in self.xslist if not item in self.oxslist]
        else:
            self.nxslist=self.xslist
            self.oxslist=[]
        #pdb.set_trace()
        fh=open('xslist.pickle','w')
        cPickle.dump(self.xslist,fh)
        fh.close()

    def copy_xray_structures(self,path=originaldir):
        os.chdir(path)
        for code in self.nxslist:
            print os.system('cp /salilab/park2/database/pdb/divided/'+code[1:3]+'/pdb'+code+'.ent.gz ./')
            print os.system('gunzip pdb'+code+'.ent.gz')
            print os.system('mv '+'pdb'+code+'.ent '+'pdb'+code+'.pdb')

    def load_xslist(self):
        os.chdir(originaldir)
        fh=open('xslist.pickle')
        self.xslist=cPickle.load(fh)
        self.nxslist=self.xslist
        fh.close()

    def preprocess(self):
        import sp
        if debug:
            for code in self.nxslist:
                self.preprocess_one(code)
        else:
            sp.task().parallel_local(self.preprocess_one,self.nxslist,10)

    def preprocess_one(self,code):
        if not os.path.isfile(originaldir+'pdb'+code+'.pdb'):
            return 0
        if os.path.isfile(finaldir+code+'.pickle'):# or code in ['3lmf','3alp','3ajn','3bcr','2c28','2c5h','2c5b','3e3n','9ga2','9ga1']:#or code in ['3a01','3a6p','2a79','3a9c','1abw']'3aej','3aem','1abw','3aen','2ajh','3ajn','3alp'
            return 0
        if debug:
            ss=SingleStructure(code)
            ss.preprocess()
        else:
            try:
                ss=SingleStructure(code)
                ss.preprocess()
            except Exception,e:
                print code
                print e
                self.errfh.write(code+'\n'+str(e))

    def collect_dict(self):
        os.chdir(finaldir)
        fl=os.listdir('./')
        fl=[f for f in fl if f[5:]=='pickle']
        pdbdict={}
        for f in fl:
            print f
            fh=open(f)
            sd=cPickle.load(fh)
            fh.close
            if 1:
                ss=SingleStructure(f[0:4])
                ss.pdict=sd
                ss.write_xray_rfactor()
                sd['type']='X'
                ss.save_dict()
            pdbdict[f[0:4]]=sd
        os.chdir('..')
        self.pdbdict=pdbdict
        self.save_dict()

    def save_dict(self):
        os.chdir('/bell2/gqdong/pdb/')
        fh=open(self.pdbdictname,'w')
        cPickle.dump(self.pdbdict,fh)
        fh.close()

    def load_dict(self):
        os.chdir('/bell2/gqdong/pdb/')
        fh=open(self.pdbdictname)
        self.pdbdict=cPickle.load(fh)
        fh.close()

    def get_seqlist(dsname):
        dsdir=decoysdir[dsname]
        fl=os.listdir(dsdir)
        env=environ()
        aln=alignment(env)
        asl=[]
        for f in fl:
            cl,sl=mymodel(f+'_r_u.pdb.HB').get_structure_chain_sequence()
            asl=asl+sl
            cl,sl=mymodel(f+'_l_u.pdb.HB').get_structure_chain_sequence()
            asl=asl+sl
        return asl

    def get_pir(self,ppdict={}):#X-xray, U-biological unit, T-transmembrane protein, M-membrane protein,S-single chain, W-non-membrane protein,1-single chain pir or none, 2-multiple chain pir, H-helical protein, B-beta protein
        #pdb.set_trace()
        if len(self.pdbdict)==0:
            self.load_dict()
        if 'similar' in ppdict:
            similarlist=ppdict['similar']
            similarfilter=True
            decoysseqlist=self.get_seqlist(similarlist[0])
            del ppdict['similar']
        else:
            similarfilter=False
        pirlist=[]
        cp=['T','H','B'] #all possible chain types
        ctype=''
        #pdb.set_trace()
        for cpi in cp:
            if cpi in ppdict['type']:
                ctype=ctype+cpi
        for code in self.pdbdict:
            codedict=self.pdbdict[code]
            collect=True
            complexpir=0
            ctype=''
            cnom=False#the chain/s can not have missing atoms/residues
            for key in ppdict:
                if key=='type': #X-xray, U-biological unit, T-transmembrane protein, M-membrane protein,S-single chain, W-non-membrane protein,1-single chain pir or none, 2-multiple chain pir, H-helical protein, B-beta protein
                    for item in ppdict[key]:
                        if item=='S':
                            if len(codedict['chains'].keys())>1:
                                collect=False
                                break
                        elif item in ['X','U','T','M','S','W','H','B']:
                            if not item in codedict[key]:
                                collect=False
                                break
                        elif item in ['2','3','4']:
                            complexpir=int(item)
                elif key=='ph':
                    #pdb.set_trace()
                    if codedict[key]<ppdict[key][0] or codedict[key]>ppdict[key][1]:
                        collect=False
                        break
                elif key=='cnom':#the single chain has no missing atoms or residues
                    cnom=True
                elif key in ['r','rfree','resolution']:
                    if codedict[key]>ppdict[key]:
                        collect=False
                        break
                else:
                    raise Exception('Do not know the property: '+key)
            if not collect:
                continue
            print code
            if complexpir==2:
                if len(codedict['chains'])>1:
                    pirlist.append((codedict['pir'],codedict['resolution'],codedict['rfree']))
                else:
                    continue
            elif complexpir==3:
                if len(codedict['chains'])==2:
                    pirlist.append((codedict['pir'],codedict['resolution'],codedict['rfree']))
                else:
                    continue
            elif complexpir==4:
                if len(codedict['chains'])>1:
                    for cc in codedict['chains']:
                        cd=codedict['chains'][cc]
                        pirlist.append((cd['pir'],codedict['resolution'],codedict['rfree']))
                else:
                    continue
            else:
                for cc in codedict['chains']:
                    cd=codedict['chains'][cc]
                    ctc=True
                    for item in ctype: #check whether the chain type matches the type dict
                        if not item in cd['type']:
                            ctc=False
                            break
                    if not ctc:
                        continue
                    print cd['missingatomnum']
                    if cnom and cd['missingatomnum']>0:
                        continue
                    if (cd['length']-cd['missingresnum'])>self.minchainlen:
                        pirlist.append((cd['pir'],codedict['resolution'],codedict['rfree']))
        pirlist.sort(key=itemgetter(1,2))
        pirstring='\n'.join([pir[0] for pir in pirlist])
        return pirstring

    def takethischain(self,seqid,decoysseqlist,csl):
        env
        for cs in csl:
            for dss in decoysseqlist:
                aln=alignment(env)
                aln.append_sequence(cs)
                aln.append_sequence(dss)
                sid=aln.id_table()
                if sid>seqid:
                    return False
        return True

    def update_S2C(self):
        pass

    def get_alltm_codes(self):
        paths=['http://pdbtm.enzim.hu/data/pdbtm_all.list','http://pdbtm.enzim.hu/data/pdbtm_alpha.list','http://pdbtm.enzim.hu/data/pdbtm_beta.list']
        keynames=['T','H','B']
        dictlist={}
        for i in range(0,len(paths)):
            path=paths[i]
            ld={}
            mlf=urllib.urlopen(path)
            mlfc=mlf.read()
            mll=mlfc.split('\n')
            mpl=[item for item in mll if len(item)>3]
            for item in mpl:
                if item[0:4] in ld:
                    ld[item[0:4]]=ld[item[0:4]]+item[-1]
                else:
                    ld[item[0:4]]=item[-1]
            dictlist[keynames[i]]=ld
        self.THBdict=dictlist

    def update_THB(self):
        self.get_alltm_codes()
        self.load_dict()
        for code in self.pdbdict:
            for dictkey in ['T','H','B']:
                if code in self.THBdict[dictkey]:
                    print code
                    self.pdbdict[code]['type']=self.pdbdict[code]['type']+dictkey
                    for cc in self.pdbdict[code]['chains']:
                        if cc in self.THBdict[dictkey][code]:
                            if 'type' in self.pdbdict[code]['chains'][cc]:
                                self.pdbdict[code]['chains'][cc]['type']=self.pdbdict[code]['chains'][cc]['type']+dictkey
                            else:
                                self.pdbdict[code]['chains'][cc]['type']=dictkey
        self.save_dict()

    def update_U(self):
        pass

    def update_ph(self):
        os.chdir(originaldir)
        env=environ()
        self.load_dict()
        for code in self.pdbdict:
            self.pdbdict[code]['ph']=sp.get_ph(code)
        self.save_dict()

class tmhpdb(entirepdb):
    def __init__(self):
        self.minchainlen=1
        self.pdbdictname='tmhpdbdict.pickle'
        self.pdbdict={}

    def convert_tochains(self,path=''):
        self.path=path
        env=environ()
        os.chdir(path)
        fl=os.listdir('./')
        fl=[item[:-4] for item in fl if item.endswith('pdb')]
        self.fl=fl
        ep=entirepdb()
        ep.load_dict()
        for f in fl:
            os.chdir(path)
            if f[0:4] not in ep.pdbdict:
                print 'skipping because of not in pdbdict'+f
                continue
            if  not os.path.isfile(path+f+'.log'):
                print 'skipping because of no log'+f
                continue
            #mdlo=mymodel(env,code=f+'.pdb')
            #mdlo.convert_tmh_tochains(path,f)
            mdl=mymodel(env,code='THAC'+f+'.pdb')
            key='THAC'+f
            self.pdbdict[key]=copy.deepcopy(ep.pdbdict[f[0:4]])
            self.pdbdict[key]['code']=key
            self.pdbdict[key]['pir']=mdl.get_structure_sequence()
            cd=mdl.build_chain_dict('TH')
            self.pdbdict[key]['chains']=cd
            mdl=[]
        self.save_dict()

    def copy_pdb(self):
        os.chdir(self.path)
        os.system('cp THAC*pdb '+finaldir)
        os.system('scp THAC*pdb gqdong@baton2:~/pdb/')

class tmhmcpdb(entirepdb):
    def __init__(self):
        self.minchainlen=1
        self.pdbdictname='tmhmcpdbdict.pickle'
        self.pdbdict={}

    def convert_tochains(self,path=''):
        self.path=path
        env=environ()
        os.chdir(path)
        fl=os.listdir('./')
        fl=[item for item in fl if item.startswith('THAC') and item.endswith('pdb')]
        fd={}
        for f in fl:
            if f[4:8] in fd:
                fd[f[4:8]].append(f)
            else:
                fd[f[4:8]]=[f]
        self.fl=fl
        ep=entirepdb()
        ep.load_dict()
        for pdbkey in fd:
            os.chdir(path)
            if pdbkey not in ep.pdbdict:
                print 'skipping because of not in pdbdict'+f
                continue
            combine_models(fd[pdbkey],'THMC'+pdbkey+'.pdb')
            mdl=mymodel(env,code='THMC'+pdbkey+'.pdb')
            key='THMC'+pdbkey
            self.pdbdict[key]=copy.deepcopy(ep.pdbdict[pdbkey])
            self.pdbdict[key]['code']=key
            self.pdbdict[key]['pir']=mdl.get_structure_sequence()
            cd=mdl.build_chain_dict('TH')
            self.pdbdict[key]['chains']=cd
            mdl=[]
        self.save_dict()

    def copy_pdb(self):
        os.chdir(self.path)
        os.system('cp THMC*pdb '+finaldir)
        os.system('scp THMC*pdb gqdong@baton2:~/pdb/')

class bioupdb(entirepdb):
    #Biounit dictionary
    def __init__(self):
        self.minchainlen=30
        self.pdbdictname='bioupdbdict.pickle'
        self.pdbdict={}

    def build_pdbdict(self,opath,tpath):
        env=environ()
        ep=entirepdb()
        ep.load_dict()
        fl=os.listdir(opath)
        for f1 in fl:
            if len(f1)!=2:
                continue
            cdir=opath+f1+'/'
            os.chdir(cdir)
            print os.system('gunzip *gz')
            print os.system('rm *r')
            print os.system('rm *.model*')
            fl2=os.listdir(cdir)
            for f2 in fl2:
                print f2
                os.chdir(cdir)
                pdbcode=f2[0:4]
                mycode='BIOU'+pdbcode+f2[-1]
                if pdbcode not in ep.pdbdict or pdbcode in ['1smv']:
                    print 'skipping because of not in pdbdict'+f2
                    continue
                mdl=model(env)
                try:
                    if not os.path.isfile(tpath+mycode+'.pdb'):
                        try:
                            mdl.read(f2)
                        except:
                            continue
                        if len(mdl.chains)<=1:
                            print "single chain structure, skipping...."
                            continue
                        self.rebuild_model(f2)
                        print os.system('cp '+f2+'r '+tpath+mycode+'.pdb')
                except Exception,e:
                    if re.search('Too many chains',str(e)):
                        print os.system('touch '+opath+f2)
                        continue
                    else:
                        traceback.print_exc()
                        pdb.set_trace()
                os.chdir(tpath)
                mdl2=mymodel(env,code=mycode)
                self.pdbdict[mycode]=copy.deepcopy(ep.pdbdict[pdbcode])
                self.pdbdict[mycode]['code']=mycode
                self.pdbdict[mycode]['pir']=mdl2.get_structure_sequence()
                cd=mdl2.build_chain_dict('')
                self.pdbdict[mycode]['chains']=cd
        self.save_dict()

    def rebuild_model(self,fp):
        ml=seperate_models(fp)
        if ml:
            combine_models(ml,fp+'r')
        else:
            print os.system('cp '+fp+' '+fp+'r')

    def copy_pdb(self,tpath):
        os.chdir(self.tpath)
        os.system('scp *pdb gqdong@baton2:~/biounits/')

    def get_dockgroundlist(self): #label,G
        lp='http://dockground.bioinformatics.ku.edu/BOUND/REPRESENTATIVES/result/represent.easy.dataset/represent_res_cluster_id_30_dimer.txt'
        pass

class loop(object): #loopchain,looplength,loopstart,loopend,looptype
    def __init__(self,env=env()):
        self.env=env
        self.nofloops=3000
        self.refpot=['$(LIB)/atmcls-melo.lib','$(LIB)/melo1-dist.lib']
        self.bm='slow'
        self.j=[]
        self.pn=8
        self.calc_rmsds=111
        self.input=''
        self.pe='local'

    def read_dssp_for_allpdb(self,dir='',filelist=[]):
        if dir:
            os.chdir(dir)
        else:
            os.chdir('/bell2/gqdong/S2C/')
        dssplist=os.listdir('./')
        if filelist:
            dssplist=filelist
        dsspdict={}
        try:
            for file in dssplist:
                print file
                if file[-4:]!='dssp':
                    continue
                try:
                    fh=open(file,'r')
                except:
                    print 'Open file problem '+file
                    continue
                dssplist=[]
                startrecord=False
                for line in fh:
                    if line[0:25]=='  #  RESIDUE AA STRUCTURE':
                        startrecord=True
                        continue
                    if startrecord:
                        try:
                            resnum=int(line[6:10])
                            chain=line[11]
                            resname=line[13]
                            ssc=line[16]
                        except Exception,e:
                            traceback.print_exc()
                            print line
                            continue
                        dssplist.append([chain,resnum,resname,ssc])
                dsspdict[file[:-5]]=dssplist
        except Exception,e:
            traceback.print_exc()
        return dsspdict

    def define_loops(self,dslist=[],filelist=[]):
        dsspdict=self.read_dssp_for_allpdb(filelist=filelist)
        if dslist:
            dslist=set(dslist)
        else:
            dslist=set([' ','T','S','B'])
        loopdict={}
        for key in dsspdict:
            dssplist=dsspdict[key]
            inloop=False
            looplist=[]
            loopchain=0
            loopstart=0
            loopend=0
            looplength=0
            looptype=0# good loop:0, broken loops:1, sterminal loops:10, end terminal loop:100
            for i in range(0,len(dssplist)):
                item=dssplist[i]
                if ((not item[3] in dslist) or (item[0]!=dssplist[i-1][0])) and inloop:
                    inloop=False
                    loopend=dssplist[i-1][1]
                    if i==1 or item[0]!=dssplist[i-1][0]:
                        looptype=looptype+100
                    if looplength!=loopend-loopstart+1:
                        looptype=looptype+1
                    looplist.append([loopchain,looplength,loopstart,loopend,looptype])
                if item[3] in dslist:
                    if not inloop:
                        inloop=True
                        loopchain=item[0]
                        loopstart=item[1]
                        if i==0 or i==(len(dssplist)-1) or item[0]!=dssplist[i-1][0]:
                            looptype=10
                        else:
                            looptype=0
                        looplength=1
                    else:
                        looplength=looplength+1
                        if i==(len(dssplist)-1):
                            loopend=dssplist[i][1]
                            looptype=looptype+100
                            if looplength!=loopend-loopstart+1:
                                looptype=looptype+1
                            looplist.append([loopchain,looplength,loopstart,loopend,looptype])
            loopdict[key]=looplist
        dslistname=''.join(dslist).strip()
        fh=open('/bell3/gqdong/statpot/loop/'+dslistname+'loop.pickle','wb')
        cPickle.dump(loopdict,fh)
        fh.close()

    def filter_loops(self,pdbset='pdb_SX.pir',Xray=2.0,Rfactor=0.25,looplength=[4,25],nonstdres=False,
                     terminal=False,broken=False,occ=1,accessibilitycutoff=[5,60],phrange=[6.5,7.5], sid=60,
                     missingatoms=False,isotempfactorstd=2,minhetamdist=1,metaliondist=3.5,maxcacadist=3.7,dslist=[' ','T','S','B']):
        os.chdir('/bell3/gqdong/statpot/loop/')
        dslistname=''.join(dslist).strip()
        fh=open('/bell3/gqdong/statpot/loop/SBTloop.pickle','rb')
        loopdict=cPickle.load(fh)
        fh.close()
        selectedloopdict={}
        env = environ()
        log.minimal()
        env.io.atom_files_directory=['/salilab/park2/database/pdb/divided/']
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        env.io.hetatm=True
        aln=alignment(env)
        from sp import *
        po=pir(self.env,runenv.pdbpirdir+pdbset)
        po.filter_pir(xray=Xray,rfactor=Rfactor,phrange=phrange,sid=sid)
        po=pir(self.env,po.pirpath[:-4]+'_'+str(Xray)+'_'+str(Rfactor)+'_'+str(phrange[0])+'-'+str(phrange[1])+'_'+str(sid)+'.pir')
        codelist=po.get_codelist()
        resatoms=get_residues_atoms()
        nofloops=0
        aacc=[]
        numofloops=0
        for code in codelist:
            if code[0:4] in loopdict:
                loops=loopdict[code[0:4]]
            else:
                continue
            selectedloops=[]
            mdl=model(env)
            mdl.read(file=code[0:4])
            accmdl=model(env)
            accmdl.read(file=code[0:4])
            accmdl.write_data('PSA',accessibility_type=2,file=None)
            asel=selection(mdl)
            nonstds=selection(mdl)-asel.only_std_residues()
            meanbiso=0
            bisomean,bstd=self.calc_model_biso(mdl)
            bisothresd=bisomean+isotempfactorstd*bstd
            clash,contact=Mymodel(env).clash_contact(code[0:4])
            self.clash=clash
            self.contact=contact
            for loop in loops:
                print loop
                sacc=self.aveacc(loop,accmdl)
                aacc.append(str(sacc))
                #try:
                if 1:
                    if not (terminal or broken):
                        if loop[4]>0:
                            continue
                    if loop[1]<looplength[0] or loop[1]>looplength[1]:
                        print 'loop length '+str(loop[1])
                        continue
                    if not (nonstdres or missingatoms):
                        try:
                            s=selection(mdl.residue_range(str(loop[2])+':'+loop[0],str(loop[3])+':'+loop[0]))
                        except:
                            continue
                        s=s.only_std_residues()
                        residuedict={}
                        for atom in s:
                            if atom.residue.name in ['ASX','GLX']:
                                continue
                            if not atom.residue.num in residuedict:
                                residuedict[atom.residue.num]=[atom.residue.name,[atom.name],atom.biso]
                            else:
                                residuedict[atom.residue.num][1].append(atom.name)
                                residuedict[atom.residue.num][2]=residuedict[atom.residue.num][2]+atom.biso
                        if len(residuedict)!=loop[1]:
                            print 'missing residue'
                            continue
                        select=True
                        for key in residuedict:
                            if (not (set(residuedict[key][1])>=set(resatoms['resatomdict'][residuedict[key][0]]))):
                                select=False
                                print 'missing atoms'
                                print residuedict[key][1]
                                print resatoms['resatomdict'][residuedict[key][0]]
                                break
                            if (residuedict[key][2]/len(residuedict[key][1]))>bisothresd:
                                print 'disorder'
                                select=False
                                break
                        if select==False:
                            continue
                    else:
                        print 'please update filter_loops function to do this, non implemented yet'
                    if sacc>accessibilitycutoff[1] or sacc<accessibilitycutoff[0]:
                        print 'Accessibility too high '+str(sacc)
                        continue
                    if self.clash_filter(loop):
                        continue
                    if self.contact_filter(loop,mdl,mindist=minhetamdist,ionmindist=metaliondist):
                        continue
                    if self.ca_distance_filter(loop,mdl,maxcacadist=maxcacadist):
                        continue
                #except Exception,e:
                #    print loop
                #    traceback.print_exc()
                print 'selected:', loop
                selectedloops.append(loop)
            if selectedloops:
                print "#######################################"+code+'    '+str(mdl.rfactor)+'   '+str(mdl.resolution)
                selectedloopdict[code[0:4]]=selectedloops
                numofloops=numofloops+len(selectedloops)
        print 'num of selected loops: '+str(numofloops)
        print 'numof structures: '+str(len(selectedloopdict))
        fh=open('/bell3/gqdong/statpot/acc','w')
        fh.write(','.join(aacc))
        fh.close()
        fh=open('/bell3/gqdong/statpot/loop/selectedloop.pickle','wb')
        cPickle.dump(selectedloopdict,fh)
        fh.close()

    def filter_bs(self,pdbset='pdb_SX.pir',Xray=2.0,Rfactor=0.25,looplength=[4,25],nonstdres=False,
                     terminal=False,broken=False,occ=1,accessibilitycutoff=[5,60],phrange=[6.5,7.5], sid=60,
                     missingatoms=False,isotempfactorstd=2,minhetamdist=1,metaliondist=3.5,maxcacadist=3.7,dslist=[' ','T','S','B']):
        os.chdir('/bell3/gqdong/statpot/loop/')
        dslistname=''.join(dslist).strip()
        fh=open('/bell3/gqdong/statpot/loop/SBTloop.pickle','rb')
        loopdict=cPickle.load(fh)
        fh.close()
        selectedloopdict={}
        env = environ()
        log.minimal()
        env.io.atom_files_directory=['/salilab/park2/database/pdb/divided/']
        env.libs.topology.read(file='$(LIB)/top_heav.lib')
        env.libs.parameters.read(file='$(LIB)/par.lib')
        env.io.hetatm=True
        aln=alignment(env)
        from sp import *
        po=pir(self.env,runenv.pdbpirdir+pdbset)
        po.filter_pir(xray=Xray,rfactor=Rfactor,phrange=phrange,sid=sid)
        po=pir(self.env,po.pirpath[:-4]+'_'+str(Xray)+'_'+str(Rfactor)+'_'+str(phrange[0])+'-'+str(phrange[1])+'_'+str(sid)+'.pir')
        codelist=po.get_codelist()
        resatoms=get_residues_atoms()
        nofloops=0
        aacc=[]
        numofloops=0
        for code in codelist:
            if code[0:4] in loopdict:
                loops=loopdict[code[0:4]]
            else:
                continue
            selectedloops=[]
            mdl=model(env)
            mdl.read(file=code[0:4])
            accmdl=model(env)
            accmdl.read(file=code[0:4])
            accmdl.write_data('PSA',accessibility_type=2,file=None)
            asel=selection(mdl)
            nonstds=selection(mdl)-asel.only_std_residues()
            meanbiso=0
            bisomean,bstd=self.calc_model_biso(mdl)
            bisothresd=bisomean+isotempfactorstd*bstd
            clash,contact=Mymodel(env).clash_contact(code[0:4])
            self.clash=clash
            self.contact=contact
            for loop in loops:
                print loop
                sacc=self.aveacc(loop,accmdl)
                aacc.append(str(sacc))
                #try:
                if 1:
                    if not (terminal or broken):
                        if loop[4]>0:
                            continue
                    if loop[1]<looplength[0] or loop[1]>looplength[1]:
                        print 'loop length '+str(loop[1])
                        continue
                    if not (nonstdres or missingatoms):
                        try:
                            s=selection(mdl.residue_range(str(loop[2])+':'+loop[0],str(loop[3])+':'+loop[0]))
                        except:
                            continue
                        s=s.only_std_residues()
                        residuedict={}
                        for atom in s:
                            if atom.residue.name in ['ASX','GLX']:
                                continue
                            if not atom.residue.num in residuedict:
                                residuedict[atom.residue.num]=[atom.residue.name,[atom.name],atom.biso]
                            else:
                                residuedict[atom.residue.num][1].append(atom.name)
                                residuedict[atom.residue.num][2]=residuedict[atom.residue.num][2]+atom.biso
                        if len(residuedict)!=loop[1]:
                            print 'missing residue'
                            continue
                        select=True
                        for key in residuedict:
                            if (not (set(residuedict[key][1])>=set(resatoms['resatomdict'][residuedict[key][0]]))):
                                select=False
                                print 'missing atoms'
                                print residuedict[key][1]
                                print resatoms['resatomdict'][residuedict[key][0]]
                                break
                            if (residuedict[key][2]/len(residuedict[key][1]))>bisothresd:
                                print 'disorder'
                                select=False
                                break
                        if select==False:
                            continue
                    else:
                        print 'please update filter_loops function to do this, non implemented yet'
                    if sacc>accessibilitycutoff[1] or sacc<accessibilitycutoff[0]:
                        print 'Accessibility too high '+str(sacc)
                        continue
                    if self.clash_filter(loop):
                        continue
                    if self.contact_filter(loop,mdl,mindist=minhetamdist,ionmindist=metaliondist):
                        continue
                    if self.ca_distance_filter(loop,mdl,maxcacadist=maxcacadist):
                        continue
                #except Exception,e:
                #    print loop
                #    traceback.print_exc()
                print 'selected:', loop
                selectedloops.append(loop)
            if selectedloops:
                print "#######################################"+code+'    '+str(mdl.rfactor)+'   '+str(mdl.resolution)
                selectedloopdict[code[0:4]]=selectedloops
                numofloops=numofloops+len(selectedloops)
        print 'num of selected loops: '+str(numofloops)
        print 'numof structures: '+str(len(selectedloopdict))
        fh=open('/bell3/gqdong/statpot/acc','w')
        fh.write(','.join(aacc))
        fh.close()
        fh=open('/bell3/gqdong/statpot/loop/selectedloop.pickle','wb')
        cPickle.dump(selectedloopdict,fh)
        fh.close()



    def clash_filter(self,loop):
        loopres=set([loop[0]+str(lnum) for lnum in range(loop[2],loop[3]+1)])
        clashlist=[]
        for item in self.clash:
            if item[0][0]=='0':
                clashlist.append(item[0][1]+item[0][2])
            if item[1][0]=='0':
                clashlist.append(item[1][1]+item[1][2])
        for item in clashlist:
            if item in loopres:
                print 'clash with '+','.join(item)+'  '
                return True
        return False

    def contact_filter(self,loop,mdl,mindist,ionmindist):
        loopres=[loop[0]+str(lnum) for lnum in range(loop[2],loop[3]+1)]
        metalset=set(metalions())
        for item in self.contact:
            try:
                if not (item[0][1]+item[0][2]) in loopres:
                    continue
                if item[1][3]=='WAT' or item[1][3]=='HOH':
                    continue
                else:
                    atomname=item[1][4].split('.')[0]
                    if item[1][3]=='MSE' and atomname=='SE':
                        element='SE'
                    elif atomname[0]=='H' or atomname[0:2] in ['1H', '2H','3H']:
                        element='H'
                    else:
                        element=mdl.residues[item[1][2]+':'+item[1][1]].atoms[atomname].element
                    if element in metalset:
                        if -float(item[2])<ionmindist:
                            print 'metal ion too close, ',item
                            return True
                    else:
                        if -float(item[2])<mindist:
                            print 'other atoms too close, ',item
                            return True
            except:
                pdb.set_trace()
        return False

    def ca_distance_filter(self,loop,mdl,maxcacadist):
        atom1=mdl.residues[str(loop[2])+':'+loop[0]].atoms['CA']
        atom2=mdl.residues[str(loop[3])+':'+loop[0]].atoms['CA']
        dist=np.sqrt((atom1.x-atom2.x)**2+(atom1.y-atom2.y)**2+(atom1.z-atom2.z)**2)
        if dist>maxcacadist*(loop[1]-1):
            print 'distance between ca atoms are too large: '+str(dist)
            return True
        return False

    def hetatm_filter(self,loopsel,haa,mindist=4): #not used anymore
        laa=np.zeros([len(loopsel),3])
        i=0
        for atom in loopsel:
            laa[i,0]=atom.x
            laa[i,1]=atom.y
            laa[i,2]=atom.z
            i=i+1
        for i in range(0,haa.shape[0]):
            for j in range(0,len(loopsel)):
                if np.sqrt(np.sum((haa[i]-laa[j])**2))<mindist:
                    print 'hetatm too close'
                    print haa[i]
                    print laa[j]
                    return False
        return True

    def aveacc(self,loop,mdl):
        for i in range(0,len(mdl.residues)):
            try:
                if int(mdl.residues[i].num)==loop[2]:
                    break
            except Exception,e:
                print e
        aveacc=0
        for j in range(i,i+loop[1]):
            aveacc=aveacc+mdl.residues[j].atoms[0].biso
        return aveacc/loop[1]

    def calc_model_biso(self,mdl):
        biso=np.zeros(len(mdl.atoms))
        for i in range(0,len(mdl.atoms)):
            biso[i]=mdl.atoms[i].biso
        return (biso.mean(), biso.std())

    def get_loop_subset(self,inputloop='selectedloop.pickle', outputloop='m5l4t20.pickle',length=range(4,21),maxn=20):
        fh=open(inputloop)
        als=pickle.load(fh)
        fh.close()
        nld={}
        length.reverse()
        for ll in length:
            k=0
            for code in als:
                if code in nld:
                    continue
                loops=als[code]
                nll=[]
                for loop in loops:
                    if loop[1]==ll:
                        nll.append(loop)
                        break
                if len(nll)>0:
                    nld[code]=nll
                    k=k+1
                    if k>=maxn:
                        break
        pdb.set_trace()
        fh=open(outputloop,'w')
        pickle.dump(nld,fh)
        fh.close()

def list2pir(pdblist):
    pdb.set_trace()

class mymodel(model):
    def __init__(self,env,code):
        if code.endswith('pdb'):
            model.__init__(self,env,file=code)
            self.code=code[:-4]
        elif code.startswith('BIOU'):
            model.__init__(self,env,file=code+'.pdb')
            self.code=code
        else:
            model.__init__(self,env,file=code)#'pdb'+code+'.pdb'
            self.code=code
        self.type='structureX'

    def read_tm_region(self,fp=''):
        fh=open(fp)
        fc=fh.read()
        fh.close()
        rer=re.findall('TMH\s+[0-9]+\s+residue\s+([0-9]+)\s+([0-9]+)',fc)
        tmrl=[]
        for item in rer:
            tmrl.append([int(item[0]),int(item[1])])
        return rer

    def select_loop_atoms(self,loops):
        s=selection()
        for loop in loops:
            try:
                s.add(selection(self.residue_range(str(loop[2])+':'+loop[0],str(loop[3])+':'+loop[0])))
            except:
                lind=self.residues[str(loop[2])+':'+loop[0]].index
                s.add(selection(self.residues[lind:(lind+loop[1])]))
        return s

    def convert_tmh_tochains(self,path='',fn=''):
        #convert different transmembrane helix to different chains
        os.chdir(path)
        tmrl=self.read_tm_region(path+fn+'.log')
        self.read(fn+'.pdb')
        cl=['A','B','C','D','E','F','G','H','I','J','K','L','M','N','O','P','Q','R','S','T']
        k=-1
        for tmr in tmrl:
            k=k+1
            self.chains[0].name=cl[k]
            s=self.select_loop_atoms([[self.chains[0].name,-1,tmr[0],tmr[1]]])
            s.write('tmregiontemp'+cl[k])
        print os.system('cat tmregiontemp* > THAC'+fn+'.pdb')
        print os.system('rm tmregiontemp*')

    def get_ph(self):
        return sp.get_ph(self.code)

    def get_structure_sequence(self):
        code=self.code
        env=environ()
        m=model(env,file=code)
        aln=alignment(env)
        aln.append_model(self,align_codes=code,atom_files=code)
        aln.write(code+'.s.pir')
        fh=open(code+'.s.pir')
        fc=fh.read()
        fh.close()
        return replace_structure_type(fc,self.type)

    def get_structure_chain_sequence(self):
        e=environ()
        cc=[]
        cs=[]
        for c in self.chains:
            filename=self.code
            print "Wrote out ",' ',filename, self.seq_id
            atom_file, align_code = c.atom_file_and_code(filename)
            #print os.system('touch '+code+'.pir')
            c.write('temp.pir', atom_file, align_code, format='PIR',
                chop_nonstd_termini=False)
            cc.append(align_code[-1])
            fh=open('temp.pir')
            fc=fh.read()
            fh.close()
            cs.append(replace_structure_type(fc,self.type))
        return [cc,cs]

    def build_chain_dict(self,ctype=''):
        cc,cs=self.get_structure_chain_sequence()
        cd={}
        for i in range(0,len(cc)):
            scd={}
            scd['type']=ctype
            scd['pir']=cs[i]
            scd['length']=len(self.chains[i].residues)
            scd['missingresnum']=self.count_missing_residues(i)
            scd['atomnum']=len(self.chains[i].atoms)
            scd['missingatomnum']=self.count_missing_atoms(i)
            cd[cc[i]]=scd
        return cd

    def count_missing_atoms(self,cn): #count atoms with coordinates -999.0
        c=self.chains[cn]
        k=0
        for atom in c.atoms:
            if atom.x==-999.0:
                k=k+1
        return k

    def count_missing_residues(self,cn):#count resiudes with coordiantes -999.0
        c=self.chains[cn]
        k=0
        for res in c.residues:
            am=1
            for atom in res.atoms:
                if atom.x!=-999.0:
                    am=0
                    break
            if am:
                k=k+1
        return k

    def build_residue_nonhet_set(self):
        nonhetdict={}
        for chain in self.chains:
            hetset=[]
            chaincode=chain.name
            for res in chain.residues:
                if not res.hetatm:
                    hetset.append(chaincode+res.num+res.pdb_name)
            nonhetdict[chaincode]=hetset
        return nonhetdict

    def build_residue_list(self):
        nonhetdict={}
        sequencedict={}
        resnumdict={}
        for chain in self.chains:
            hetset=[]
            resnumlist=[]
            sequencel=''
            chaincode=chain.name
            for res in chain.residues:
                hetset.append(chaincode+res.num+res.pdb_name)
                sequencel=sequencel+res.code
                resnumlist.append(res.num)
            if chaincode in resnumdict:
                continue
            resnumdict[chaincode]=resnumlist
            nonhetdict[chaincode]=hetset
            sequencedict[chaincode]=sequencel
        #if len(sequencedict)!=len(self.chains):
        #    print 0
        #    raise Exception('Modeller model has multiple chains with the same chain code '+self.code)
        return [nonhetdict,sequencedict,resnumdict]

class SingleStructure(object):
    def __init__(self,code):
        self.code=code
        self.type='X'

    def preprocess(self):
        self.errdict=self.precalc_errorscale()
        self.load_original_model()
        self.flip_sidechain()
        self.molprobity_analysis()
        self.pdict=self.load_molprobity()
        combine_dict1d(self.pdict,self.errdict)
        self.read_full_sequence()
        self.merge_sequence()
        self.build_alignment()
        self.complete_pdb()
        self.load_final_model()
        self.pdict['chains']=self.get_chain_pir()
        self.pdict['pir']=self.finalm.get_structure_sequence()
        self.calc_errorscale()
        self.write_xray_rfactor()
        os.chdir(finaldir)
        self.save_dict()

    def save_dict(self):
        ph=open(self.code+'.pickle','wb')
        cPickle.dump(self.pdict,ph)
        ph.close()

    def precalc_errorscale(self):
        os.chdir(originaldir)
        f='pdb'+self.code+'.pdb'
        if not os.path.isfile('pdb'+self.code+'.pdb'):
            raise Exception('Pdb file does not exist '+self.code)
        n=0
        pno='NULL'
        rwrfactor=0;
        pnn=self.code
        print self.code
        cdict={}
        collect=1
        # read the data from file and calculate the scale factor
        pdbfile=open(f)
        pline=pdbfile.read()
        #free rfactor
        if frfactorre.search(pline):
            frfactorrer=frfactorre.search(pline)
            frfactor=float(frfactorrer.group(1))
        elif rfactorre.search(pline):
            rfactorrer=rfactorre.search(pline)
            frfactor=float(rfactorrer.group(1))
        else:
            frfactor=0
        cdict['rfree']=frfactor
        #r factor
        if rfre.search(pline):
            rfrer=rfre.search(pline)
            rfactor=float(rfrer.group(1))
        else:
            rfactor=0
        cdict['r']=rfactor
        #resolution
        if resolre.search(pline):
            resolrer=resolre.search(pline)
            resol=float(resolrer.group(1))
        else:
            resol=None
        cdict['resolution']=resol
        #number of reflections
        if norre.search(pline):
            norrer=norre.search(pline)
            nor=np.abs(float(norrer.group(1)))
        else:
            nor=None
        cdict['numofreflections']=nor
        #average B factor
        abfactor=0
        if 1:
            brer=bre.findall(pline)
            bl=0
            #print brer
            for j in range(0,len(brer)):
                #print brer[j][0]
                if brer[j]:
                    bl=bl+float(brer[j])
    #          if len(brer)>0:
            abfactor=bl/len(brer)
        cdict['aveB']=abfactor
        print  'abfactor: '+str(abfactor)
        #completeness
        completeness=0
        if completenessre.search(pline):
            completenessrer=completenessre.search(pline)
            completeness=float(completenessrer.group(1))
        if completeness==0:
            completeness=100
        if completeness<50:
            completeness=100-completeness
        cdict['completeness']=completeness
        #print  'completeness: '+str(completeness)
        #print  'number of reflection: '+str(nor)
        atomnrer=atomre.findall(pline)
        atomn=len(atomnrer)
        cdict['atomnum']=atomn
        return cdict

    def load_original_model(self):
        os.chdir(originaldir)
        if not os.path.isfile('pdb'+self.code+'.pdb'):
            raise Exception('Pdb file does not exist '+self.code)
        env=environ()
        env.io.hetatm=False
        self.originalm=mymodel(env,self.code)

    def flip_sidechain(self):
        code=self.code
        os.chdir(originaldir)
        if not os.path.isdir(self.code):
            os.mkdir(self.code)
        file='pdb'+code+'.pdb'
        print os.system('cp '+file+' ./'+self.code+'/')
        inputdir=originaldir+self.code+'/'
        outputdir=flippeddir
        print os.system('/bell1/home/gqdong/software/molprobity3/cmdline/reduce-build '+inputdir+' '+outputdir)
        print os.system('rm -r '+self.code)

    def molprobity_analysis(self):
        os.chdir(flippeddir)
        file='pdb'+self.code+'.pdb'
        if not os.path.isdir(file[:-4]):
            os.mkdir(file[:-4])
        print os.system('cp '+file+' ./'+file[:-4]+'/' )
        print os.system('/bell1/home/gqdong/software/molprobity3/cmdline/oneline-analysis '+file[:-4]+' >'+file[:-4]+'.molprobity')
        print os.system('rm -r '+file[:-4])

    def load_molprobity(self):
        f='pdb'+self.code+'.molprobity'
        f=open(f)
        fc=f.readline()
        fc=f.readline()
        fc=f.readline()
        f.close()
        fd={}
        if fc[-1]=='\n':
            fc=fc[:-1]
        fl=fc.split(':')
        fd['code']=fl[0][3:-4]
        nl=[ 'chains', 'residues', 'nucacids', 'resolution', 'r', 'rfree', 'clashscore', 'clashscoreB<40', 'minresol', 'maxresol', 'n_samples', 'pct_rank', 'pct_rank40', 'cbeta>0.25', 'numCbeta', 'maxCbeta', 'medianCbeta', 'meanCbeta', 'rota<1%', 'numRota', 'ramaOutlier', 'ramaAllowed', 'ramaFavored', 'numRama', 'pct_badbonds', 'pct_badangles', 'MolProbityScore', 'Mol_pct_rank']
        for i in range(0,len(nl)):
            if fl[i+1]:
                fd[nl[i]]=float(fl[i+1])
            else:
                fd[nl[i]]=0
        return fd

    def read_full_sequence(self):
        code=self.code
        bdir='/bell2/gqdong/S2C/'
        mlf=open(bdir+code+'.sc')
        mlfc=mlf.read()
        mlf.close()
        sc=re.findall('(SEQCRD.*)\\n',mlfc)
        sc=sc[1:]
        sl={}
        cc=''
        for item in sc:
            if item[15:18]=='MSE':
                item=item[0:15]+'MET'+item[18:]
            if item[7]!=cc:
                if cc:
                    sl[cc]={'ss':ss,'fs':fs,'reslist':reslist,'restype':restype}
                cc=item[7]
                fs=''
                ss=''
                restype=[]
                reslist=[]
                psn=-999999
            csn=int(item[19:24])
            if csn==psn:
                restype[-1].append(item[11:14])
                restype[-1].append(item[15:18])
                continue
            fs=fs+item[9]
            reslist.append(item[7]+item[25:31].strip()+item[15:18])
            restype.append([item[11:14],item[15:18]])
            if item[15:18]=='---':
                ss=ss+'-'
            else:
                ss=ss+item[9]
            psn=csn
        sl[cc]={'ss':ss,'fs':fs,'reslist':reslist,'restype':restype}
        self.fullseq=sl

    def merge_sequence(self):
        nreslist,nresseq,resnumdict=self.originalm.build_residue_list()
        sl={}
        for cc in nreslist:
            if cc not in self.fullseq:
                sl[cc]={'ss':nresseq[cc],'fs':nresseq[cc]}
                continue
            ssm=''
            fsm=''
            cnrl=nreslist[cc]
            crl=self.fullseq[cc]['reslist']
            ss=self.fullseq[cc]['ss']
            fs=self.fullseq[cc]['fs']
            rli=[-1 for i in range(0,len(cnrl))]
            allmissing=True
            for i in range(0,len(cnrl)):
                if cnrl[i] in crl:
                    rli[i]=crl.index(cnrl[i])
                    if rli[i]!=-1:
                        allmissing=False
            if allmissing:
                sl[cc]={'ss':nresseq[cc],'fs':nresseq[cc]}
                continue
            #identify all the gaps in sequence
            gl=[]
            gs=0
            for i in range(0,len(cnrl)):
                if rli[i]==-1:
                    if not gs:
                        gl.append([i])
                    gs=1
                else:
                    if gs:
                        gl[-1].append(i)
                    gs=0
            if gs:
                gl[-1].append(i+1)
            #go through gaps and try to fill in the gaps using the residue types
            #pdb.set_trace()
            for g in gl:
                ngl=g[1]-g[0]
                resolved=False
                if g[0]==0:
                    fsep=rli[g[1]]
                    fssp=0
                elif g[1]==len(cnrl):
                    fssp=rli[g[0]-1]+1
                    fsep=len(crl)
                else:
                    fssp=rli[g[0]-1]+1
                    fsep=rli[g[1]]
                fsgl=fsep-fssp
                if ngl==fsgl:
                    fss=fs[fssp:fsep].replace('X','\w')
                    if re.match(fss,nresseq[cc][g[0]:g[1]]):
                        rli[g[0]:g[1]]=range(fssp,fsep)
                        resolved=True
                elif fsgl==0:
                    rli[g[0]:g[1]]=[-9 for i in range(g[0],g[1])]#a valid gap to be inserted here
                    resolved=True
                elif fsgl>ngl:
                    #pdb.set_trace()
                    fsgs=fs[fssp:fsep]
                    fss=fs[fssp:fsep].replace('X','\w')
                    ngs=nresseq[cc][g[0]:g[1]]
                    for gresi in range(g[0],g[1]):
                        if g[0]==0:
                            try:
                                rpte=int(resnumdict[cc][g[1]])-int(resnumdict[cc][gresi])
                                rpts=len(fsgs)-rpte+1
                            except:
                                rer=re.search('([0-9]+)',resnumdict[cc][g[1]])
                                rer2=re.search('([0-9]+)',resnumdict[cc][gresi])
                                #pdb.set_trace()
                                rpte=int(rer.group(1))-int(rer2.group(1))
                                if rpte==0:
                                    rpte=1
                                rpts=len(fsgs)-rpte+1
                        else:
                            try:
                                rpts=int(resnumdict[cc][gresi])-int(resnumdict[cc][g[0]-1])
                            except:
                                rer=re.search('([0-9]+)',resnumdict[cc][gresi])
                                rer2=re.search('([0-9]+)',resnumdict[cc][g[0]-1])
                                #pdb.set_trace()
                                rpts=int(rer.group(1))-int(rer2.group(1))
                                if rpts==0:
                                    rpts=1
                        if not re.match(fsgs[rpts-1],nresseq[cc][gresi]):
                            resolved=False
                            break
                        else:
                            resolved=True
                            rli[gresi]=rli[g[0]-1]+rpts
                if not resolved:
                    raise Exception('can not resolve the gap between native sequence and the S2C sequence at gap pos'+nreslist[cc][fssp]+' '+nreslist[cc][fsep]+'  code '+self.code)
            #based on the resolved gaps we will rebuild the full sequencences
            for i in range(0,len(cnrl)):
                if (i==(len(cnrl)-1)) and rli[i]<(len(crl)-1):
                    ssm=ssm+nresseq[cc][i]+''.join(['-' for j in range(rli[i]+1,len(crl))])
                    fsm=fsm+nresseq[cc][i]+fs[(rli[i]+1):]
                    break
                if i==0 and rli[i]>0:
                    ssm=ssm+''.join(['-' for j in range(0,rli[i])])
                    fsm=fsm+fs[0:rli[i]]
                if (i<(len(cnrl)-1)) and (rli[i+1]-rli[i])>1 and rli[i]!=-9 and rli[i+1]!=-9:
                    ssm=ssm+nresseq[cc][i]+''.join(['-' for j in range(rli[i],rli[i+1]-1)])
                    fsm=fsm+nresseq[cc][i]+fs[(rli[i]+1):rli[i+1]]
                else:
                    ssm=ssm+nresseq[cc][i]
                    fsm=fsm+nresseq[cc][i]
            sl[cc]={'ss':ssm,'fs':fsm}
            #pdb.set_trace()
        self.fullseq=sl
        #pdb.set_trace()

    def read_full_sequence_nohetatm(self):
        code=self.code
        bdir='/bell2/gqdong/S2C/'
        nonhetdict=self.originalm.build_residue_nonhet_set()
        od=copy.deepcopy(nonhetdict)
        mlf=open(bdir+code+'.sc')
        mlfc=mlf.read()
        mlf.close()
        sc=re.findall('(SEQCRD.*)\\n',mlfc)
        sc=sc[1:]
        sl={}
        cc=''
        k=0
        al=0
        resindexdict={}
        resnamedict={}
        for item in sc:
            if item[7]!=cc:
                if cc:
                    sl[cc]={'ss':ss,'fs':fs}
                    al=al+len(od[cc])
                ki=0
                cc=item[7]
                resnamedict[cc]=[]
                if cc in nonhetdict:
                    cd=nonhetdict[cc]
                else:
                    nonhetdict[cc]=[]
                    cd=nonhetdict[cc]
                fs=''
                ss=''
            fs=fs+item[9]
            if item[15:18]=='MSE':
                item=item[0:15]+'MET'+item[18:]
            if ((item[7]+item[25:31].strip()+item[15:18]) in cd) and item[15:18]!='---':
                cd.remove(item[7]+item[25:31].strip()+item[15:18])
                ss=ss+item[9]
                resindexdict[item[7]+item[25:31].strip()+item[15:18]]=ki
                k=k+1
            else:
                ss=ss+'-' #item[15:18]=='---'
            resnamedict[cc].append([item[11:14],item[15:18]])
            ki=ki+1
        sl[cc]={'ss':ss,'fs':fs}
        al=al+len(od[cc])
        #pdb.set_trace()
        #deal with the case when the residue is in the full sequence but not in the pdb sequence of the alignment file
        for cc in sl:
            if len(nonhetdict[cc])>0:
                for res in nonhetdict[cc]:
                    ri=od[cc].index(res)
                    if ri>0 and (od[cc][ri-1] in resindexdict):
                        rci=resindexdict[od[cc][ri-1]]+1
                        #pdb.set_trace()
                        if res[-3:] in resnamedict[cc][rci]:
                            sl[cc]['ss']=sl[cc]['ss'][0:rci]+sl[cc]['fs'][rci]+sl[cc]['ss'][rci+1:]
                            k=k+1
                            nonhetdict[cc].remove(res)
                    elif (len(od[cc])>ri+1) and (od[cc][ri+1] in resindexdict):
                        rci=resindexdict[od[cc][ri+1]]-1
                        if res[-3:] in resnamedict[cc][rci]:
                            sl[cc]['ss']=sl[cc]['ss'][0:rci]+sl[cc]['fs'][rci]+sl[cc]['ss'][rci+1:]
                            k=k+1
                            nonhetdict[cc].remove(res)
        #dealwith the case when the residue in not even in the full sequence

        #pdb.set_trace()
        if al!=k:
            raise Exception('can not build correct alignment file from the full sequence file '+self.code)
        dl=0
        for key in nonhetdict:
            dl=dl+len(nonhetdict[key])
        if dl>0:
            raise Exception('Atoms in pdb file are missing in alignment files, probably due to DNA or RNA sequence'+self.code)
        self.fullseq=sl

    def build_alignment(self):
        os.chdir(flippeddir)
        code=self.code
        sl=self.fullseq
        cc,cs=self.originalm.get_structure_chain_sequence()
        k=0
        lc=0
        cl=[]
        for i in range(0,len(cc)):
            if cc[i] in sl:
                cl.append(cc[i])
                if k==0:
                    ae=cs[i].split(':')
                    sp=ae[2:4]
                k=k+1
                lc=i
        ae=cs[lc].split(':')
        ep=ae[4:6]
        ns='>P1;'+code+'8\nsequence::     : :     : :::-1.00:-1.00\n'
        nss='>P1;'+code+'\nstructure'+self.type+':'+code+':'+sp[0]+':'+sp[1]+':'+ep[0]+':'+ep[1]+':::-1.00:-1.00\n'
        mrc=[]
        for c in cl:
            ns=ns+sl[c]['fs']+'/'
            nss=nss+sl[c]['ss']+'/'
            mrc.append(sl[c]['ss'].count('-'))
        ns=ns[:-1]+'*\n\n'
        nss=nss[:-1]+'*\n\n'
        rs=[]
        fh=open(code+'.ali','w')
        fh.write(nss+ns)
        fh.close()
        self.missingresnum=mrc
        #pdb.set_trace()

    def complete_pdb(self):
        path=flippeddir
        tarpath=finaldir
        code=self.code
        os.chdir(path)
        e=environ()
        #io=e.io
        #e.io.hetatm=False
        e.libs.topology.read(file='$(LIB)/top_heav.lib')
        e.libs.parameters.read(file='$(LIB)/par.lib')
        aln=alignment(e)
        aln.append(code+'.ali')
        mdl=model(e)
        mdl.generate_topology(aln[-1])
        mdl.transfer_xyz(aln)
        mdl.write(tarpath+'pdb'+code+'.pdb')

    def load_final_model(self):
        os.chdir(finaldir)
        env=environ()
        self.finalm=mymodel(env,self.code)

    def get_chain_pir(self):
        cc,cs=self.finalm.get_structure_chain_sequence()
        cd={}
        for i in range(0,len(cc)):
            scd={}
            scd['pir']=cs[i]
            scd['length']=len(self.finalm.chains[i].residues)
            scd['missingresnum']=self.finalm.count_missing_residues(i)
            scd['atomnum']=len(self.finalm.chains[i].atoms)
            scd['missingatomnum']=self.finalm.count_missing_atoms(i)
            cd[cc[i]]=scd
        return cd

    def calc_errorscale(self):
        codedict=self.pdict
        collect=True
        if not (codedict['atomnum'] and codedict['numofreflections'] and codedict['completeness'] and codedict['resolution'] and (codedict['r']+codedict['rfree'])):
            codedict['aveerr']=0
            codedict['errorscale']=0
            return 0
        if codedict['rfree']:
            aveerr=((codedict['atomnum']/codedict['numofreflections'])**0.5)*((codedict['completeness']/100)**(-0.33333333))*codedict['rfree']*codedict['resolution']
            #print  'number of atom: '+str(atomn)
        else:
            aveerr=((codedict['atomnum']/codedict['numofreflections'])**0.5)*((codedict['completeness']/100)**(-0.33333333))*codedict['r']*codedict['resolution']
        codedict['aveerr']=aveerr
        if codedict['aveB']:
            codedict['errorscale']=aveerr/codedict['aveB']
        else:
            codedict['errorscale']=0
        #line=line[:-6]+scalefactor[0:5]

    def write_xray_rfactor(self):
        self.pdict['pir']=replace_xray_rfactor(self.pdict['pir'],self.pdict['resolution'],self.pdict['errorscale']*1000)
        for c in self.pdict['chains']:
            cd=self.pdict['chains'][c]
            cd['pir']=replace_xray_rfactor(cd['pir'],self.pdict['resolution'],self.pdict['errorscale']*1000)

def combine_models(modellist,fmn):
    if len(modellist)==1:
        print os.system('cp '+modellist[0]+' '+fmn)
        return 0
    env=environ()
    env.io.hetatm=False
    cl=['A','B','C','D','E','F','G','H','I','J','K',
        'L','M','N','O','P','Q','R','S','T','U','V','W','X','Y','Z',
        '0','1','2','3','4','5','6','7','8','9','a','b','c','d','e',
        'f','g','h','i','j','k','l','m','n','o','p','q','r','s','t',
        'u','v','w','x','y','z','/','[',']','{','}','(',')','&','%',
        '#','@','!','<','>']
    k=0
    ci=0
    mnl=''
    for mn in modellist:
        mdl=model(env,file=mn)
        for i in range(0,len(mdl.chains)):
            chain=mdl.chains[i]
            chain.name=cl[ci]
            ci=ci+1
            if ci>=76:
                raise Exception('Too many chains')
        mdl.write(mn)
        mnl=mnl+' '+mn
        k=k+1
    print os.system('cat '+mnl+' >'+fmn)

def seperate_models(fp):
    fh=open(fp)
    fc=fh.read()
    fh.close()
    nm=re.findall('MODEL',fc)
    if len(nm)<2:
        return []
    fl=fc.split('ENDMDL')
    k=0
    mnl=[]
    for mf in fl:
        if len(mf)<400:
            continue
        mn=fp+'.model'+str(k)
        fh=open(mn,'w')
        fh.write(mf)
        fh.close()
        mnl.append(mn)
        k=k+1
    return mnl

def complete_pdb_dir(path,tarpath):
    os.chdir(path)
    fl=os.listdir('./')
    fl=[item for item in fl if item[-3:]=='pdb']
    inputlist=[]
    for f in fl:
        code=f[3:-4]
        inputlist.append([code,path,tarpath])
    import sp
    sp.task().parallel_local(complete_pdb_one,inputlist)

def replace_structure_type(s,type='structureX'):
    sls=s.split('\n')
    k=0
    for i in range(0,len(sls)):
        if re.search('.*:.*:.*:',sls[i]):
            break
    tl=sls[i].split(':')
    tl[0]=type
    sls[i]=':'.join(tl)
    s='\n'.join(sls)
    return s

def replace_xray_rfactor(s,xray,rfactor):
    sls=s.split('\n')
    k=0
    for i in range(0,len(sls)):
        if re.search('.*:.*:.*:',sls[i]):
            break
    tl=sls[i].split(':')
    tl[-1]='{:2.2f}'.format(rfactor)
    tl[-2]='{:2.2f}'.format(xray)
    sls[i]=':'.join(tl)
    s='\n'.join(sls)
    return replace_structure_type(s)

def precalc_errorscale(pdbdir):
    os.chdir(pdbdir)
    fl=os.listdir('./')
    fl=[item for item in fl if item[-3:]=='pdb']
    n=0
    pno='NULL'
    rwrfactor=0;
    collect=1
    pdbdict={}
    for f in fl:
        pnn=f[3:-4]
        print pnn
        pdbdict[pnn]={}
        collect=1
        # read the data from file and calculate the scale factor
        pdbfile=open(f)
        pline=pdbfile.read()
        #free rfactor
        if frfactorre.search(pline):
            frfactorrer=frfactorre.search(pline)
            frfactor=float(frfactorrer.group(1))
        elif rfactorre.search(pline):
            rfactorrer=rfactorre.search(pline)
            frfactor=float(rfactorrer.group(1))
        else:
            frfactor=None
        pdbdict[pnn]['rfree']=frfactor
        #r factor
        if rfre.search(pline):
            rfrer=rfre.search(pline)
            rfactor=float(rfrer.group(1))
        else:
            rfactor=None
        pdbdict[pnn]['r']=rfactor
        #resolution
        if resolre.search(pline):
            resolrer=resolre.search(pline)
            resol=float(resolrer.group(1))
        else:
            resol=None
        pdbdict[pnn]['resolution']=resol
        #number of reflections
        if norre.search(pline):
            norrer=norre.search(pline)
            nor=np.abs(float(norrer.group(1)))
        else:
            nor=None
        pdbdict[pnn]['numofreflections']=nor
        #average B factor
        abfactor=0
        if 1:
            brer=bre.findall(pline)
            bl=0
            #print brer
            for j in range(0,len(brer)):
            #print brer[j][0]
                if brer[j]:
                    bl=bl+float(brer[j])
  #          if len(brer)>0:
            abfactor=bl/len(brer)
        pdbdict[pnn]['aveB']=abfactor
        print  'abfactor: '+str(abfactor)
        #completeness
        completeness=0
        if completenessre.search(pline):
            completenessrer=completenessre.search(pline)
            completeness=float(completenessrer.group(1))
        if completeness==0:
            completeness=100
        if completeness<50:
            completeness=100-completeness
        pdbdict[pnn]['completeness']=completeness
        #print  'completeness: '+str(completeness)
        #print  'number of reflection: '+str(nor)
        atomnrer=atomre.findall(pline)
        atomn=len(atomnrer)
        pdbdict[pnn]['atomnum']=atomn
    out=file(pdbdir+'pdbdict2.pickle',"wb")
    ccPickle.dump(pdbdict,out)
    out.close()

def load_molprobity_dir(path):
    os.chdir(path)
    fl=os.listdir('./')
    fl=[item for item in fl if len(item)>11 and item[-11:]=='.molprobity']
    pdbdict={}
    for f in fl:
        fd=load_molprobity_singlefile(f)
        pdbdict[fd['code']]=fd
    out=file('pdbdict1.pickle',"wb")
    ccPickle.dump(pdbdict,out)
    out.close()
    return pdbdict

def combine_dict2d(dict1,dict2):
    for key in dict1:
        if key in dict2:
            combine_dict1d(dict1[key],dict2[key])
        else:
            del dict1[key]

def combine_dict1d(dict1,dict2):
    for key in dict2:
        if (key in dict1) and dict1[key]:
            continue
        else:
            dict1[key]=dict2[key]

def update_pdball():
    print os.system('scp gqdong@baton2:/netapp/database/pdb/remediated/pdball.pir '+basedir+'pdbpir/pdball.pir')
    pir_err(pfd+'pdball.pir')

def gen_pdb_A():
    e = environ()
    e.io.atom_files_directory = ['/salilab/park2/database/pdb/divided/']
    k=0
    pdball=modfile.File(pfd+'pdb_A_temp.pir','w')
    for root, dir, file in os.walk('/salilab/park2/database/pdb/divided/'):
        for fname in file:
            fpath=root+'/'+fname
            try:
                m = model(e, file=fpath)
                for c in m.chains:
                    if c.filter(minimal_chain_length=1, minimal_resolution=99.0,
                        minimal_stdres=1, chop_nonstd_termini=True,
                        structure_types='structureX structureN'):
                        filename = re.findall('pdb\w{4}',fname)[0]
                        print "Wrote out ",str(k),' ',filename, m.seq_id
                        atom_file, align_code = c.atom_file_and_code(filename)
                        c.write(pdball, atom_file, align_code, format='PIR',
                            chop_nonstd_termini=True)
                        k+=1
            except:
                print e
    pir_cc(pfd+'pdb_A_temp.pir',pfd+'pdb_A.pir')
    pir_err(pfd+'pdb_A.pir')

def opm_mpl():
    #Function to read the complete list of transmembrane proteins from OPM
    mlf=urllib.urlopen('http://opm.phar.umich.edu/proteins.php')
    mlfc=mlf.read()
    mpl=re.findall('Text\[\d*\]=\[\"([0-9a-z]{4})',mlfc)
    for pn in mpl:
        print pn
        time.sleep(1)
        mlf=urllib.urlopen('http://opm.phar.umich.edu/protein.php?pdbid='+pn)
        mlfc=mlf.read()
        opt=re.findall('Related PDB Sum entries</b></td>\\r\\n.*\\r\\n',mlfc)
        mpl=mpl+re.findall('>([0-9a-z]{4})</a>',opt[0])
    fh=open(basedir+'pdbpir/mpl','w')
    fh.write(' '.join(mpl))
    fh.close()
    return mpl

def opm_tmcl():
    sfn=open(basedir+'pdbpir/tmplts','r')
    fc=sfn.read()
    cfn=open(basedir+'pdbpir/tmcl','w')
    mpl=re.findall('>([0-9a-z]{4})([^>]*)>',fc)
    tmcl=''
    for pn in mpl:
        print pn[0]
        opt=re.findall('\\n([A-Z0-9]{1})\s{1}Tilt',pn[1])
        for item in opt:
            tmcl=tmcl+pn[0]+item+" "
    cfn.write(tmcl)
    cfn.close()
    sfn.close()
    return tmcl

def pcl(of):
    sfn=open(of,'r')
    cfn=open(basedir+'pdbpir/pcl','w')
    pcll=''
    for fc in sfn:
        if re.match('\s*([0-9a-z]{4})\s*;',fc):
            rer=re.match('\s*([0-9a-z]{4})\s*;',fc)
            pcll=pcll+rer.group(1)+" "
    cfn.write(pcll)
    cfn.close()
    sfn.close()
    return pcll

def opm_tmpl():
        #Function to read in all the of transmembrane proteins from OPM, with the transmembrane segments and tilt angles.
    mlf=urllib.urlopen('http://opm.phar.umich.edu/classes.php?type=1')
    mlfc=mlf.read()
    mpl=re.findall('Text\[\d*\]=\[\"([0-9a-z]{4})',mlfc)
    fh=open(basedir+'pdbpir/tmpl','w') # this file contains the codes of unique transmembrane proteins defined by OPM
    fh.write(' '.join(mpl))
    fh.close()
    sfn=open(basedir+'pdbpir/tmplts','w') # this file contains the tilt angles and the transmembrane segments of unique transmembrane proteins defined by OPM
    for pn in mpl:
        time.sleep(1)
        print pn
        sfn.write('>'+pn+'\n')
        mlf=urllib.urlopen('http://opm.phar.umich.edu/protein.php?pdbid='+pn)
        mlfc=mlf.read()
        opt=re.findall('Related PDB Sum entries</b></td>\\r\\n.*\\r\\n',mlfc)
        mpl=mpl+re.findall('>([0-9a-z]{4})</a>',opt[0])
        mplts=re.findall('<b>([A-Z]{1})</b>\s*-\s*(Tilt:\s*[0-9]{1,3})&#176;\s*-\s*(Segments.*)</td>',mlfc)
        for mpltsi in mplts:
            sfn.write(' '.join(mpltsi)+';\n')
    sfn.close()
    fh=open(basedir+'pdbpir/tmpla','w') #this file contains the codes of all transmembrane proteins from OPM
    fh.write(' '.join(mpl))
    fh.close()
    return mpl

def gen_pir(pdbfiles):
    pfd=basedir+'pdbpir/'
    print 'Generating pdb file ' + pdbfiles
    env = environ()
    aln = alignment(env)
    pdbfile='pdb_'+pdbfiles+'.pir'
    if re.search('([0-9]{1}\.[0-9]{1})A',pdbfiles):
        rer=re.search('([0-9]{1}\.[0-9]{1})A',pdbfiles)
        rsvs=rer.group(1)
        stresolution=float(rsvs)
        print stresolution
    else:
        stresolution=99.0
    if not re.search('_[0-9]{2}$',pdbfiles):
        pdb1name=pdbfile
    else:
        pdb1name=pdbfile[:-7]+'.pir'
    if pdbfiles[0]=='X':
        pir_xrf(pfd+'pdball.pir',pfd+pdb1name, structuretype='structureX',structureresolution=stresolution)
    elif pdbfiles[0]=='N':
        pir_xrf(pfd+'pdball.pir',pfd+pdb1name, structuretype='structureN',structureresolution=stresolution)
    elif pdbfiles[0]=='M':
        pir_ief(pfd+'pdb_X'+pdb1name[5:],pfd+pdb1name,pfd+'mpl')
    elif pdbfiles[0]=='W':
        pir_ief(pfd+'pdb_X'+pdb1name[5:],pfd+pdb1name,pfd+'mpl',inc=0)
    elif pdbfiles[0:2]=='TC':
        pir_ief(pfd+'pdb_X'+pdb1name[6:],pfd+pdb1name,pfd+'tmcl',fnl=5)
    elif pdbfiles[0]=='T':
        pir_ief(pfd+'pdb_X'+pdb1name[5:],pfd+pdb1name,pfd+'tmpl')
    if pdbfiles[0:2]=='AX':
        pir_xrf(pfd+'pdb_A.pir',pfd+pdb1name, structuretype='structureX',structureresolution=stresolution)
    elif pdbfiles[0:2]=='AM':
        pir_ief(pfd+'pdb_AX'+pdb1name[6:],pfd+pdb1name,pfd+'mpl')
    elif pdbfiles[0:2]=='AW':
        pir_ief(pfd+'pdb_AX'+pdb1name[6:],pfd+pdb1name,pfd+'mpl',inc=0)
    elif pdbfiles[0:2]=='AT':
        pir_ief(pfd+'pdb_AX'+pdb1name[6:],pfd+pdb1name,pfd+'tmpl')
    elif pdbfiles[0:2]=='AC':
        pir_ief(pfd+'pdb_AX'+pdb1name[6:],pfd+pdb1name,pfd+'pcl')
    elif pdbfiles[0]=='A':
        gen_pdb_A()
    elif pdbfiles[0:2]=='SX':
        pir_xrf(pfd+'pdb_S.pir',pfd+pdb1name, structuretype='structureX',structureresolution=stresolution)
    elif pdbfiles[0:2]=='SN':
        pir_xrf(pfd+'pdb_S.pir',pfd+pdb1name, structuretype='structureN',structureresolution=stresolution)
    elif pdbfiles[0:2]=='SM':
        pir_ief(pfd+'pdb_SX'+pdb1name[6:],pfd+pdb1name,pfd+'mpl')
    elif pdbfiles[0:2]=='SW':
        pir_ief(pfd+'pdb_SX'+pdb1name[6:],pfd+pdb1name,pfd+'mpl',inc=0)
    elif pdbfiles[0:3]=='STC':
        pir_ief(pfd+'pdb_SX'+pdb1name[7:],pfd+pdb1name,pfd+'tmcl',fnl=5)
    elif pdbfiles[0:2]=='ST':
        pir_ief(pfd+'pdb_SX'+pdb1name[6:],pfd+pdb1name,pfd+'tmpl')
    else:
        print 'pp.py unsupported pdb database file type (example X_2.2A_60): '+pdbfiles
    if re.search('_[0-9]*$',pdbfiles):
        sic=int(pdbfiles[-2:])
        pir_sif(pfd+pdb1name, pfd+pdbfile,si=sic)
    return pdbfile

def pir_ief(inputpirfile,outputpirfile,mlist,inc=1,fnl=4):
    print 'converting '+inputpirfile+' to '+outputpirfile
    if not (os.access(inputpirfile,os.F_OK)):
        print "please generate "+inputpirfile+" first"
    f=open(mlist,'r' )
    fc=f.read()
    out=file(outputpirfile,"w")
    pirhead=re.compile(">P1;(\S*)$")
    inputfc=open(inputpirfile)
    n=0
    pno=0
    for line in inputfc:
        if pirhead.match(line):
            pno=0
            mg=pirhead.match(line)
            pn=mg.group(1)
            pnn=pn[0:fnl]
            if bool(inc)==bool(re.search(pnn,fc,re.IGNORECASE)):
                print pnn
                pno=1
                n=n+1
        if pno:
            out.write(line)
    print n

def pir_xrf(inputpirfile,outputpirfile,structuretype='structureX',structureresolution=99.0):
    print 'converting '+inputpirfile+' to '+outputpirfile
    env = environ()
    env.io.atom_files_directory=['/salilab/park2/database/pdb/divided/']
    aln = alignment(env)
    input = modfile.File(inputpirfile, 'r')
    output = modfile.File(outputpirfile, 'w')
#    rs=1
#    b=1
#    while rs and b:
#        try:
#            b=aln.read_one(input, allow_alternates=True)
#            print a[0].code
#            rs=0
#        except Exception:
#            rs=1
    while aln.read_one(input, allow_alternates=True):
        if aln[0].resolution < structureresolution and aln[0].prottyp==structuretype:
            aln.write(output, alignment_format='PIR')
#        while rs and b:
#            try:
#                b=aln.read_one(input, allow_alternates=True)
#                print a[0].code
#                rs=0
#            except Exception:
#                rs=1
    input.close()
    output.close()

def pir_sif(inputpirfile,outputpirfile,si=30):
    print 'converting '+inputpirfile+' to '+outputpirfile
    log.verbose()
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

def pir_err(inputpirfile):
    print os.system('cp '+inputpirfile+' '+inputpirfile+'.temp')
    pirfile=inputpirfile+'.temp'
    pdbdir='/salilab/park2/database/pdb/divided/'
    os.chdir('/bell2/gqdong/temp/')
    out=file(inputpirfile,"w")
    error=file('/bell2/gqdong/statpot/pdbpir/errorlist',"w")
    pirhead=re.compile(">P1;(\S*)$")
    resolre=re.compile("\s{2}RESOLUTION RANGE HIGH \(ANGSTROMS\)\s*:\s*([0-9\.]+)")
    completenessre=re.compile("\s{2}COMPLETENESS.+\(%\)\s*:\s*([0-9\.]+)")
    rfre=re.compile("R VALUE.{1,40}:\s*(0\.\d+)")
    rfactorre=re.compile("FREE R VALUE\s*\(WORKING SET\)\s*:\s*(0\.\d+)")
    frfactorre=re.compile("\s{2}FREE R VALUE\s{14}.{1,20}:\s*(0\.\d+)")
    norre=re.compile("NUMBER OF REFLECTIONS.{1,30}:\s*([0-9\.]+)")
    atomre=re.compile("\\nATOM.{52}1.00\s*|\\nHETATM.{50}1.00\s*")
    bre=re.compile("\\nATOM.{57}\s?([0-9\.]+)\s*")
    abfactorre=re.compile("\s{2}MEAN B VALUE\s*\(OVERALL, A\*\*2\)\s*:\s*([0-9\.]+)")
    n=0
    input=open(pirfile)
    pno='NULL'
    rwrfactor=0;
    collect=1
    for line in input:
        if pirhead.match(line):
            mg=pirhead.match(line)
            pn=mg.group(1)
            pnn=pn[0:4]
            print pnn
            collect=1
            # read the data from file and calculate the scale facbasedir+'pdbpirtor
            ppathgz=pdbdir+pnn[1:3]+'/pdb'+pnn+'.ent.gz'
            if not (os.access(ppathgz,os.F_OK)):
                error.write(ppathgz)
                print ppathgz
                continue
            print os.system('cp '+ppathgz+' /bell2/gqdong/temp/')
            if not (os.access('/bell2/gqdong/temp/pdb'+pnn+'.ent',os.F_OK)):
                print os.system('gunzip '+' /bell2/gqdong/temp/pdb'+pnn+'.ent.gz')
            pdbfile=open('/bell2/gqdong/temp/pdb'+pnn+'.ent')
            pline=pdbfile.read()
            if frfactorre.search(pline):
                frfactorrer=frfactorre.search(pline)
                frfactor=float(frfactorrer.group(1))
            elif rfactorre.search(pline):
                rfactorrer=rfactorre.search(pline)
                frfactor=float(rfactorrer.group(1))
            elif rfre.search(pline):
                rfrer=rfre.search(pline)
                frfactor=float(rfrer.group(1))
            else:
                collect=0
            #print  'free R factor: '+str(frfactor)
            if resolre.search(pline):
                resolrer=resolre.search(pline)
                resol=float(resolrer.group(1))
            else:
                collect=0
            if norre.search(pline):
                norrer=norre.search(pline)
                nor=float(norrer.group(1))
                if nor<=0:
                    collect=0
            else:
                collect=0
            #print  'resolution: '+str(resol)
            abfactor=0
            if 1:
                brer=bre.findall(pline)
                bl=0
                #print brer
                for j in range(0,len(brer)):
                    #print brer[j][0]
                    if brer[j]:
                        bl=bl+float(brer[j])
  #              if len(brer)>0:
                abfactor=bl/len(brer)
            if abfactor==0:
                collect=0
            print  'abfactor: '+str(abfactor)
            completeness=0
            if completenessre.search(pline):
                completenessrer=completenessre.search(pline)
                completeness=float(completenessrer.group(1))
            if completeness==0:
                completeness=100
            if completeness<50:
                completeness=100-completeness
            #print  'completeness: '+str(completeness)
            #print  'number of reflection: '+str(nor)
            atomnrer=atomre.findall(pline)
            atomn=len(atomnrer)
            #print  'number of atom: '+str(atomn)
            print n
            #print sf
            if collect==0:
                if abfactor==0:
                    sf=0.001
                else:
                    sf=0.1/abfactor
                error.write(line)
                collect=1
            else:
                sf=((atomn/nor)**0.5)*((completeness/100)**(-0.33333333))*frfactor*resol/abfactor
            scalefactor=str(sf*1000)
            print 'Mean error: '+str(sf*abfactor)+'   Scale factor:'+str(scalefactor)
            rwrfactor=1
            out.write(line)
            continue
        if rwrfactor:
            rwrfactor=0
            n=n+1
            line=line[:-6]+scalefactor[0:5]
            print os.system('rm '+' /bell2/gqdong/temp/pdb'+pnn+'.ent*')
            out.write(line+'\n')
            continue
        if collect:
            out.write(line)
   # if n>1:
   #     break;
    input.close()
    print os.system('rm  '+inputpirfile+'.temp')

#combine chains in the same proteins
def pir_cc(inputfile,outputfile):
    print 'converting '+inputfile+' to '+outputfile
    pdball=open(inputfile,'r')
    out=open(outputfile,'w')
    pirhead=re.compile(">P1;(\S*)$")
    sequence=re.compile("[a-zA-Z0-9]*:[a-zA-Z0-9_\.]*:[a-zA-Z0-9\s\-]*:[a-zA-Z0-9\s]*:([a-zA-Z0-9\s\-]*:[a-zA-Z0-9\s]*):")
    n=0
    pno='NULL'
    pnno=''
    wl=''
    srer=''
    lastline=''
    rwrfactor=0;
    line=pdball.readline()
    while line:
        if pirhead.match(line):
            mg=pirhead.match(line)
            pn=mg.group(1)
            pnn=pn[0:4]
            if pnn!=pnno:
                wl=re.sub('1234:1234',srer,wl)
                out.write(wl+lastline+'\n')
                name=line
                wl=line
                b=pdball.readline()
                sre=sequence.match(b)
                srer=sre.group(1)
                c=re.sub(srer,'1234:1234',b)
                wl=wl+c
                line=pdball.readline()
                while not re.search('\*', line):
                    wl=wl+line
                    line=pdball.readline()
                lastline=line
            else:
                d=re.sub('\*','/',lastline)
                wl=wl+d
                line=pdball.readline()
                sre=sequence.match(line)
                srer=sre.group(1)
                line=pdball.readline()
                while not re.search('\*', line):
                    wl=wl+line
                    line=pdball.readline()
                lastline=line
                n=n+1
#               if n>2:
#                   break
            pnno=pnn
        line=pdball.readline()

def gen_sc():
    e = environ()
    e.io.atom_files_directory = ['/salilab/park2/database/pdb/divided/']
    k=0
    pdball=modfile.File(pfd+'pdb_S.pir','w')
    for root, dir, file in os.walk('/salilab/park2/database/pdb/divided/'):
        for fname in file:
            fpath=root+'/'+fname
            try:
                m = model(e, file=fpath)
                if len(m.chains) > 1:
                    continue
                for c in m.chains:
                    if c.filter(minimal_chain_length=30, minimal_resolution=99.0,
                            minimal_stdres=30, chop_nonstd_termini=True,
                            structure_types='structureX structureN'):
                        filename = re.findall('pdb\w{4}',fname)[0]
                        print "Wrote out ",str(k),' ',filename, m.seq_id
                        atom_file, align_code = c.atom_file_and_code(filename)
                        c.write(pdball, atom_file, align_code, format='PIR',
                        chop_nonstd_termini=True)
                        k+=1
            except:
                print e

def decoys_pir(ddir):
    e = environ()
    e.io.atom_files_directory = [ddir]
    k=0
    pirfile=modfile.File(ddir+'all.pir','w')
    files=os.listdir(ddir)
    for name in files:
        try:
            m = model(e, file=ddir+name)
            for c in m.chains:
                print name
                c.write(pirfile, name[:-4], name[:-4], format='PIR',chop_nonstd_termini=True)
                k+=1
                print k;
        except:
            print e

def pir2fsa(inputfile):
    input=file(pfd+inputfile)
    output=file(inputfile[:-3]+'fsa','w')
    n=0;
    pirhead=re.compile(">P1;(.*)$")
    linematch=re.compile("C; P")
    sly=0
    county=0
    for line in input:
        if county:
            if linematch.match(line):
                county=0
                output.write('\n')
            else:
                output.write(line)
        if sly==1:
            sly=0
            county=1
        if pirhead.match(line):
            sly=1;
            output.write(line)
            #print line
    return 0
