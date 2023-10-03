from __future__ import print_function
import os
import pickle
from collections import defaultdict
import pdb

loopsetname='fiser'
loopdir='/bell3/gqdong/statpot/loop/'+loopsetname
rloopdir='/netapp/sali/gqdong/loop/'+loopsetname
if not os.path.isdir(loopdir):
    os.mkdir(loopdir)
print(os.system('ssh gqdong@baton2 mkdir '+rloopdir))
tl=open('fiser-targets.txt').read().split('\n')
ll=[]
ld=defaultdict(list)
for t in tl[1:]:
    if len(t)<7:
        continue
    ts=t.split()
    ld[ts[2]].append((ts[0],[['', int(ts[2]), int(ts[1]), int(ts[1])+int(ts[2])-1]]))
    ll.append((ts[0],[['', int(ts[2]), int(ts[1]), int(ts[1])+int(ts[2])-1]]))
    #print(os.system('cp fiser-structures/'+ts[0]+'.sali.pdb '+os.path.join(loopdir,'pdb'+ts[0]+'.ent')))
    #print(os.system('scp fiser-structures/'+ts[0]+'.sali.pdb gqdong@baton2:'+os.path.join(rloopdir,'pdb'+ts[0]+'.ent')))


pickle.dump(dict(ll),open(loopdir+'.pickle','wb'))
for l in ld:
    pickle.dump(dict(ld[l]),open(loopdir+l+'.pickle','wb'))
    loopdir='/bell3/gqdong/statpot/loop/'+loopsetname+l
    rloopdir='/netapp/sali/gqdong/loop/'+loopsetname+l
    if not os.path.isdir(loopdir):
        os.mkdir(loopdir)
    print(os.system('ssh gqdong@baton2 mkdir '+rloopdir))
    for sl in ld[l]:
        key=sl[0]
        print(os.system('cp fiser-structures/'+key+'.sali.pdb '+os.path.join(loopdir,'pdb'+key+'.ent')))
        print(os.system('scp fiser-structures/'+key+'.sali.pdb gqdong@baton2:'+os.path.join(rloopdir,'pdb'+key+'.ent')))
