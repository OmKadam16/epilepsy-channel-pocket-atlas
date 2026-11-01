import sys, numpy as np, collections
f,ids=sys.argv[1],[int(x) for x in sys.argv[2].split(',')]
atoms=[];sph=collections.defaultdict(list)
for l in open(f):
    if l.startswith('ATOM') and l[12:16].strip()[0]!='H':
        atoms.append((int(l[22:26]),np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])])))
    elif l.startswith('HETATM'):
        sph[int(l[22:26])].append(np.array([float(l[30:38]),float(l[38:46]),float(l[46:54])]))
X=np.array([a[1] for a in atoms]);print('pocket ids',min(sph),max(sph),len(sph))
for p in ids:
    if p not in sph: print(p,'missing');continue
    S=np.array(sph[p]);d=np.linalg.norm(X[:,None]-S[None],axis=2).min(1)
    r=sorted({atoms[i][0] for i in np.where(d<=4)[0]});print(p,len(S),'spheres',len(r),'res:',r)
