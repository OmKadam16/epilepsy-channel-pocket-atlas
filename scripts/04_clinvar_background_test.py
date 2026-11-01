import pickle,numpy as np,statsmodels.api as sm
from scipy.stats import fisher_exact
rows=pickle.load(open('rows.pkl','rb'));rng=np.random.default_rng(1)
def run(name,R):
    rsa=np.array([r[3] for r in R]);pk=np.array([r[4] for r in R],int);y=np.array([r[5] for r in R],int)
    a=((pk==1)&(y==1)).sum();b=((pk==1)&(y==0)).sum();c=((pk==0)&(y==1)).sum();d=((pk==0)&(y==0)).sum()
    orr,p=fisher_exact([[a,b],[c,d]])
    X=sm.add_constant(np.column_stack([pk,rsa]));m=sm.Logit(y,X).fit(disp=0)
    # permutation within RSA quintile bins
    bins=np.digitize(rsa,np.quantile(rsa,[.2,.4,.6,.8]));obs=y[pk==1].mean();cnt=0;N=5000
    for _ in range(N):
        yy=y.copy()
        for b_ in range(5):
            i=np.where(bins==b_)[0];yy[i]=rng.permutation(y[i])
        cnt+=yy[pk==1].mean()>=obs
    print(f'{name}: n={len(R)} pocket={pk.sum()} path_in_pocket={obs:.2f} vs outside={y[pk==0].mean():.2f} | Fisher OR={orr:.2f} p={p:.3g} | adj OR(pocket)={np.exp(m.params[1]):.2f} p={m.pvalues[1]:.3g} | RSA-perm p={(cnt+1)/(N+1):.4f} | mean RSA pocket {rsa[pk==1].mean():.2f} vs {rsa[pk==0].mean():.2f}')
for g in ['SCN1A','SCN2A','KCNQ2']: run(g,[r for r in rows if r[0]==g])
run('Nav pooled',[r for r in rows if r[0]!='KCNQ2'])
run('All pooled',rows)
