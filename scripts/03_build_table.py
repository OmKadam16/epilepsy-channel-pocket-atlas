import json,numpy as np,collections
from Bio.PDB import PDBParser, ShrakeRupley
import statsmodels.api as sm
MAX={'ALA':129,'ARG':274,'ASN':195,'ASP':193,'CYS':167,'GLN':225,'GLU':223,'GLY':104,'HIS':224,'ILE':197,'LEU':201,'LYS':236,'MET':224,'PHE':240,'PRO':159,'SER':155,'THR':172,'TRP':285,'TYR':263,'VAL':174}
def bits(f):
    h=open(f).read().strip();b=''.join(format(int(c,16),'04b') for c in h);return {i for i,c in enumerate(b) if c=='1'}
cfg={'SCN1A':('../data/pockets/SCN1A_7DTD_fpocket_pockets.pdb',{3,50,101}),'SCN2A':('../data/pockets/SCN2A_6J8E_fpocket_pockets.pdb',{1,48,79,95,99}),
     'KCNQ2':('../data/pockets/KCNQ2_7CR0_fpocket_pockets.pdb',{23,33,42,67,36,66,72,79,47,60,63,73,38,69,74,75})}
import sys
rows=[]
for g,(f,ids) in cfg.items():
    atoms=[];sph=collections.defaultdict(list);lines=[]
    for l in open(f):
        if l.startswith('ATOM') and l[12:16].strip()[0]!='H' and (g!='KCNQ2' or l[21]=='A'):
            lines.append(l)
        elif l.startswith('HETATM'): sph[int(l[22:26])].append([float(l[30:38]),float(l[38:46]),float(l[46:54])])
    open(f'{g}_prot.pdb','w').writelines(lines+['END\n'])
    s=PDBParser(QUIET=True).get_structure('x',f'{g}_prot.pdb');ShrakeRupley().compute(s,level='R')
    res={r.id[1]:(r.resname,r.sasa) for r in s[0]['A']}
    allatoms=[(int(l[22:26]),[float(l[30:38]),float(l[38:46]),float(l[46:54])]) for l in open(f) if l.startswith('ATOM') and l[12:16].strip()[0]!='H' and (g!='KCNQ2' or l[21] in 'ABCD')]
    X=np.array([a[1] for a in allatoms]);lin=set()
    for p in ids:
        S=np.array(sph[p]);d=np.linalg.norm(X[:,None]-S[None],axis=2).min(1);lin|={allatoms[i][0] for i in np.where(d<=4)[0]}
    path=bits(f'../data/variants/{g}.path.hex');ben=bits(f'../data/variants/{g}.ben.hex') if g!='KCNQ2' else set()
    for n,(rn,sa) in res.items():
        rows.append((g,n,rn,min(sa/MAX.get(rn,200),1.0),n in lin,n in path,n in ben))
    print(g,'modeled',len(res),'pocket-lining',len(lin&set(res)),'path pos in structure',len(path&set(res)),'ben pos in structure',len(ben&set(res)))
import pickle;pickle.dump(rows,open('rows.pkl','wb'))
