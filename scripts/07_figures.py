"""Figures for the README. Run after 03_build_table.py (needs rows.pkl)."""
import pickle, numpy as np, os
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
rows = pickle.load(open('rows.pkl','rb'))
def bits(f):
    h=open(f).read().strip(); b=''.join(format(int(c,16),'04b') for c in h); return {i for i,c in enumerate(b) if c=='1'}
gn = {g: bits(f'../data/variants/{g}.gn.hex') for g in ['SCN1A','SCN2A','KCNQ2']}
BLUE,GRAY,INK,MUTED,ORANGE='#2a78d6','#b8b7b0','#0b0b0b','#52514e','#eb6834'
NAMES={'SCN1A':'Nav1.1 (SCN1A)','SCN2A':'Nav1.2 (SCN2A)','KCNQ2':'KCNQ2'}
plt.rcParams.update({'font.family':'DejaVu Sans','axes.spines.top':False,'axes.spines.right':False,'axes.edgecolor':MUTED,'text.color':INK,'axes.labelcolor':MUTED,'xtick.color':MUTED,'ytick.color':MUTED})
os.makedirs('../figures',exist_ok=True)

# ---- Fig 1: pathogenic fraction, pocket vs elsewhere (labelled positions only)
res={}
for g in NAMES:
    R=[r for r in rows if r[0]==g and (r[5] or r[1] in gn[g])]
    pk=[r for r in R if r[4]]; ot=[r for r in R if not r[4]]
    res[g]=(np.mean([bool(r[5]) for r in pk]),len(pk),np.mean([bool(r[5]) for r in ot]),len(ot))
p_adj={'SCN1A':'0.018','SCN2A':'0.009','KCNQ2':'0.0009'}; p_reg={'SCN1A':'0.007','SCN2A':'0.0004','KCNQ2':'0.31'}
fig,ax=plt.subplots(figsize=(8,4.6),dpi=160); x=np.arange(3); w=0.36
for i,g in enumerate(NAMES):
    a,na,b,nb=res[g]
    ax.bar(i-w/2-0.01,a*100,w,color=BLUE,label='In druggable pocket' if i==0 else None)
    ax.bar(i+w/2+0.01,b*100,w,color=GRAY,label='Elsewhere' if i==0 else None)
    ax.text(i-w/2-0.01,a*100+1.5,f'{a*100:.0f}%\n(n={na})',ha='center',va='bottom',fontsize=8,color=INK)
    ax.text(i+w/2+0.01,b*100+1.5,f'{b*100:.0f}%\n(n={nb})',ha='center',va='bottom',fontsize=8,color=INK)
    ax.text(i,-17,f'p (burial-adj.) = {p_adj[g]}\np (burial + region) = {p_reg[g]}',ha='center',va='top',fontsize=8,color=MUTED)
ax.set_xticks(x); ax.set_xticklabels(NAMES.values(),fontsize=10); ax.tick_params(axis='x',length=0,pad=2)
ax.set_ylim(0,125); ax.set_yticks([0,25,50,75,100]); ax.set_ylabel('Positions that are pathogenic (%)')
ax.set_title('Pathogenic positions are over-represented in druggable pockets',loc='left',fontsize=12,fontweight='bold',color=INK)
ax.legend(frameon=False,loc='upper right',ncol=2,fontsize=9); ax.yaxis.grid(True,color='#e6e5e0',lw=0.8); ax.set_axisbelow(True)
fig.text(0.01,0.01,'Pathogenic = ClinVar P/LP missense positions; comparison = gnomAD missense positions (AC>=2). Exploratory pilot.',fontsize=7,color=MUTED)
fig.subplots_adjust(bottom=0.27,left=0.09,right=0.98,top=0.9); fig.savefig('../figures/fig1_enrichment.png'); plt.close(fig)

# ---- Fig 2: sequence tracks
fig,axs=plt.subplots(3,1,figsize=(10,5.4),dpi=160)
for ax,g in zip(axs,NAMES):
    R=[r for r in rows if r[0]==g]; pos=np.array([r[1] for r in R])
    path=[r[1] for r in R if r[5]]; ctrl=[r[1] for r in R if (not r[5]) and r[1] in gn[g]]; pock=[r[1] for r in R if r[4]]
    ax.vlines(path,0.55,1.0,color=ORANGE,lw=0.7); ax.vlines(ctrl,0.05,0.5,color=GRAY,lw=0.7)
    ax.vlines(pock,-0.35,0,color=BLUE,lw=1.4)
    ax.set_xlim(pos.min()-5,pos.max()+5); ax.set_ylim(-0.4,1.05); ax.set_yticks([])
    ax.set_title(NAMES[g],loc='left',fontsize=10,fontweight='bold'); ax.spines['left'].set_visible(False)
    ax.set_xlabel('Residue number' if g=='KCNQ2' else '',fontsize=8)
for xx,(lab,col) in zip([0.55,0.72,0.89],[('pathogenic (ClinVar)',ORANGE),('population (gnomAD)',MUTED),('pocket lining',BLUE)]):
    fig.text(xx,0.955,'\u25A0 '+lab,color=col,fontsize=8.5,ha='left',va='center')
fig.text(0.01,0.955,'Variants and druggable pockets along each channel',ha='left',va='center',fontsize=12,fontweight='bold')
fig.tight_layout(rect=(0,0,1,0.94)); fig.savefig('../figures/fig2_sequence_tracks.png'); plt.close(fig)

# ---- Fig 3: KCNQ2 retigabine control
ret={104,236,239,240,243,272,299,300,303,304,305,308,312}
lin={93,96,97,100,211,212,217,233,236,237,240,243,272,299,300,305}
both=ret&lin
fig,ax=plt.subplots(figsize=(9,2.9),dpi=160)
for r in range(70,331):
    pass
ax.vlines(sorted(ret),0.55,1.0,color=ORANGE,lw=2.2); ax.vlines(sorted(lin),0.0,0.45,color=BLUE,lw=2.2)
for r in sorted(both): ax.plot([r,r],[0.45,0.55],color=INK,lw=1)
ax.set_xlim(70,330); ax.set_ylim(-0.22,1.15); ax.set_yticks([]); ax.set_xlabel('KCNQ2 residue number (7CR0 / 7CR2)')
ax.spines['left'].set_visible(False)
ax.text(71,1.06,'Retigabine contacts, 7CR2 (13 residues)',color=ORANGE,fontsize=9,fontweight='bold')
ax.text(71,-0.19,"fpocket's best druggable pocket, apo 7CR0 (16 residues)",color=BLUE,fontsize=9,fontweight='bold')
ax.text(158,0.5,f'{len(both)} of 13 shared (54%): 236, 240, 243, 272, 299, 300, 305',ha='center',va='center',fontsize=8.5,color=INK,bbox=dict(fc='white',ec='none',pad=2))
ax.set_title('Positive control: fpocket finds the retigabine site without being told',loc='left',fontsize=11,fontweight='bold')
fig.tight_layout(); fig.savefig('../figures/fig3_retigabine_control.png'); plt.close(fig)
print('ok')
