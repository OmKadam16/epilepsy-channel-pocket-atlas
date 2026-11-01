import pickle,numpy as np
rows=pickle.load(open('rows.pkl','rb'));rng=np.random.default_rng(2)
def run(name,R,W):
    pos=np.array([r[1] for r in R]);rsa=np.array([r[3] for r in R]);pk=np.array([r[4] for r in R],int);y=np.array([r[5] for r in R],int)
    rb=np.digitize(rsa,np.quantile(rsa,[.33,.66]));key=(pos//W)*3+rb
    obs=y[pk==1].mean();N=5000;cnt=0
    groups=[np.where(key==k)[0] for k in np.unique(key)]
    for _ in range(N):
        yy=y.copy()
        for i in groups: yy[i]=rng.permutation(y[i])
        cnt+=yy[pk==1].mean()>=obs
    print(f'{name} (window {W}, RSA terciles): obs {obs:.2f}  perm p={(cnt+1)/(N+1):.4f}')
for g in ['SCN1A','SCN2A','KCNQ2']:
    for W in ([150,75] if g!='KCNQ2' else [60,30]): run(g,[r for r in rows if r[0]==g],W)
