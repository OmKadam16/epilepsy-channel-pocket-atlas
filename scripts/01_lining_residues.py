import sys, numpy as np, collections
f=sys.argv[1]
atoms=[];sph=collections.defaultdict(list)
for l in open(f):
    if l.startswith('ATOM') and l[12:16].strip()[0]!='H':
        atoms.append((l[21],int(l[22:26]),l[17:20],np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])))
    elif l.startswith('HETATM'):
        sph[int(l[22:26])].append(np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]))
X=np.array([a[3] for a in atoms])
out={}
for p,s in sorted(sph.items()):
    S=np.array(s); d=np.linalg.norm(X[:,None,:]-S[None,:,:],axis=2).min(axis=1)
    res=sorted({(atoms[i][0],atoms[i][1]) for i in np.where(d<=4.0)[0]})
    out[p]=res
import json
json.dump({p:res for p,res in out.items()},open('lining.json','w'))
for p in [23,33,42,67,36,66,72,79,47,60,63,73]:
    r=out[p]; chains=sorted({c for c,_ in r})
    print(p,len(sph[p]),'spheres',len(r),'res chains',chains,'resnums',sorted({n for _,n in r})[:40])
