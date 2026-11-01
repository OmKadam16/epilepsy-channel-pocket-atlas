import pickle,numpy as np,statsmodels.api as sm
from scipy.stats import fisher_exact
rows=pickle.load(open('rows.pkl','rb'));rng=np.random.default_rng(3)
def bits(f):
    h=open(f).read().strip();b=''.join(format(int(c,16),'04b') for c in h);return {i for i,c in enumerate(b) if c=='1'}
gn={g:bits('../data/variants/'+g+'.gn.hex') for g in ['SCN1A','SCN2A','KCNQ2']}
def run(name,R,W=None):
    pos=np.array([r[1] for r in R]);rsa=np.array([r[3] for r in R]);pk=np.array([r[4] for r in R],int);y=np.array([r[5] for r in R],int)
    a=((pk==1)&(y==1)).sum();b=((pk==1)&(y==0)).sum();c=((pk==0)&(y==1)).sum();d=((pk==0)&(y==0)).sum()
    orr,p=fisher_exact([[a,b],[c,d]])
    X=sm.add_constant(np.column_stack([pk,rsa]));m=sm.Logit(y,X).fit(disp=0)
    rb=np.digitize(rsa,np.quantile(rsa,[.33,.66]));key=(pos//(W or 10**6))*3+rb
    groups=[np.where(key==k)[0] for k in np.unique(key)];obs=y[pk==1].mean();cnt=0;N=5000
    for _ in range(N):
        yy=y.copy()
        for i in groups: yy[i]=rng.permutation(y[i])
        cnt+=yy[pk==1].mean()>=obs
    print(f'{name}: n={len(R)} pocket={pk.sum()} path-frac in pocket {obs:.2f} vs {y[pk==0].mean():.2f} | Fisher OR={orr:.2f} p={p:.3g} | RSA-adj OR={np.exp(m.params[1]):.2f} p={m.pvalues[1]:.3g} | perm(RSA{",win%d"%W if W else ""}) p={(cnt+1)/(N+1):.4f}')
for g in ['SCN1A','SCN2A','KCNQ2']:
    R=[(r[0],r[1],r[2],r[3],r[4],1 if r[5] else 0) for r in rows if r[0]==g and (r[5] or r[1] in gn[g])]
    print(g,'labelled',len(R),'path',sum(r[5] for r in R),'ctrl',sum(1-r[5] for r in R),'overlap path&gnomAD',sum(1 for r in rows if r[0]==g and r[5] and r[1] in gn[g]))
    run(g,R);run(g,R,150 if g!='KCNQ2' else 60)
