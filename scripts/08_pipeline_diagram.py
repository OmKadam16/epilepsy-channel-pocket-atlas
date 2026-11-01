"""Pipeline diagram -> figures/fig0_pipeline.png"""
import matplotlib; matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch
INK,MUTED,BLUE,ORANGE,GRAY='#0b0b0b','#52514e','#2a78d6','#eb6834','#8a8982'
fig,ax=plt.subplots(figsize=(12,5.4),dpi=160); ax.set_xlim(0,12.6); ax.set_ylim(0,5.4); ax.axis('off')
def box(x,y,w,h,title,sub,col):
    ax.add_patch(FancyBboxPatch((x,y),w,h,boxstyle='round,pad=0.02,rounding_size=0.12',fc='white',ec=col,lw=1.8))
    ax.text(x+w/2,y+h-0.2,title,ha='center',va='top',fontsize=9.5,fontweight='bold',color=INK)
    ax.text(x+w/2,y+h-0.55,sub,ha='center',va='top',fontsize=7.8,color=MUTED,linespacing=1.35)
def arrow(x1,y1,x2,y2,col=GRAY):
    ax.add_patch(FancyArrowPatch((x1,y1),(x2,y2),arrowstyle='-|>',mutation_scale=12,color=col,lw=1.5,shrinkA=0,shrinkB=0))
ax.set_ylim(0,5.4); 
ax.text(0.05,5.25,'STRUCTURE SIDE',fontsize=8,color=BLUE,fontweight='bold'); ax.text(0.05,3.55,'VARIANT SIDE',fontsize=8,color=ORANGE,fontweight='bold'); ax.text(0.05,1.85,'BURIAL',fontsize=8,color=GRAY,fontweight='bold')
box(0.05,3.95,2.1,1.2,'Channel structures','PDB 7DTD (Nav1.1)\n6J8E (Nav1.2)\n7CR0 (KCNQ2)',BLUE)
box(2.6,3.95,2.1,1.2,'fpocket on Vecura','Voronoi pocket detection\nscore, druggability, volume',BLUE)
box(5.15,3.95,2.1,1.2,'Druggable pockets','druggability >= 0.5\nsymmetry copies collapsed',BLUE)
box(7.7,3.95,2.1,1.2,'Lining residues','protein atoms within 4 \u00c5\nof pocket alpha-spheres',BLUE)
for x in (2.15,4.7,7.25): arrow(x+0.03,4.55,x+0.42,4.55,BLUE)
box(0.05,2.2,2.1,1.2,'ClinVar','missense, pathogenic /\nlikely pathogenic',ORANGE)
box(2.6,2.2,2.1,1.2,'gnomAD v4','missense seen in >= 2\nalleles (comparison)',ORANGE)
box(5.15,2.2,4.65,1.2,'Unique positions','per gene, one count per residue (symmetry copies not counted twice)',ORANGE)
arrow(2.18,2.8,2.57,2.8,ORANGE); arrow(4.73,2.8,5.12,2.8,ORANGE)
box(5.15,0.45,4.65,1.2,'Burial (relative SASA)','Shrake-Rupley on the same structures; used as a covariate',GRAY)
ax.add_patch(FancyBboxPatch((10.5,0.45),2.0,4.7,boxstyle='round,pad=0.02,rounding_size=0.12',fc='#f4f3ef',ec=INK,lw=1.8))
ax.text(11.5,4.95,'Enrichment tests',ha='center',va='top',fontsize=9.5,fontweight='bold',color=INK)
ax.text(11.5,4.45,'Fisher exact\n\nlogistic regression\n(adjusted for burial)\n\npermutation within\nburial + region bins\n\nchecked against the\nretigabine site (7CR2)',ha='center',va='top',fontsize=7.8,color=MUTED,linespacing=1.3)
arrow(9.83,4.55,10.47,4.55,BLUE); arrow(9.83,2.8,10.47,2.8,ORANGE); arrow(9.83,1.05,10.47,1.05,GRAY)
ax.text(0.05,0.2,'Exploratory pilot',fontsize=8,color=MUTED,style='italic')
fig.tight_layout(pad=0.4); fig.savefig('../figures/fig0_pipeline.png'); print('ok')
