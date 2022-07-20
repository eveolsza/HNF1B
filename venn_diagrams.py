
#%%

import os
import matplotlib.pyplot as plt
from matplotlib_venn import venn2, venn2_circles
import pandas as pd
import numpy as np

# change working directory
os.chdir('/home/eve/Dropbox/doctorado/hnf1b/bulk_rna-seq/plots/venn diagrams/')

# función para ver qué genes están en común (de la lista total) y de los markers que estamos evaluando y luego plotea el diagrama de venn
# seleccionar toda la función y correrla al principio
def venn(het, hom, title):
    
    global genes_het
    global genes_hom
    global genes_both
    
    # elegir la columna con los nombres de los genes, del data frame total
    het = het['GENEID'].tolist()
    hom = hom['GENEID'].tolist()
    
    # ver qué genes están en común (de la lista total) y los imprime en la consola de python para poder verlos   
    genes_het = set(het) - set(hom)
#    print('genes het: ', genes_het)
    genes_hom = set(hom) - set(het)
#    print('genes hom: ', genes_hom)
    genes_both = set(het) & set(hom)
#    print('genes both: ', genes_both)
    
    # ver qué markers (de la lista que le pasamos) están en común   
    markers_in_het = set(genes_het) & set(markers)
    markers_in_hom = set(genes_hom) & set(markers)
    markers_in_both = set(genes_both) & set(markers)
    
    # diagrama de venn 
    fig = plt.figure()
    v = venn2([set(het), set(hom)], set_labels = ('Het vs WT', 'Hom vs WT'), set_colors = ('white', 'white'), alpha = 0.5)
    c = venn2_circles(subsets = [set(het), set(hom)], linewidth=1, color='b')
    
    # título del plot
    plt.title(title, y=1.2)
    
# =============================================================================
#     # anotar los markers en cada sector 
#     plt.annotate('\n'.join(markers_in_het), fontweight='bold', xy=v.get_label_by_id('10').get_position() + np.array([-0.05, 0]), xytext=(-50,0),
#                  ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.1),
#                  arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.3',color='gray'))
#     
#     if v.get_label_by_id('11') != None:    
#         plt.annotate('\n'.join(markers_in_both), fontweight='bold', xy=v.get_label_by_id('11').get_position() + np.array([0, 0.05]), xytext=(-80,120),
#                      ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.1),
#                      arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.3',color='gray'))
#        
#     plt.annotate('\n'.join(markers_in_hom), fontweight='bold', xy=v.get_label_by_id('01').get_position() + np.array([0.05, 0.05]), xytext=(120,0),
#                  ha='center', textcoords='offset points', bbox=dict(boxstyle='round,pad=0.5', fc='white', alpha=0.1),
#                  arrowprops=dict(arrowstyle='->', connectionstyle='arc3,rad=0.3',color='gray'))
#     
#     plt.show()
# =============================================================================
    
    # guardar el plot con el título como nombre del archivo
    fig.savefig(title+'.pdf', bbox_inches='tight')
#acá termina la función


   
# markers a evaluar, modificarlos a piaccere
markers = ['SOX4','SOX11', 'HMGA2', 'PDX1', 'HNF1B', 'SOX9', 'PTF1A', 'NKX6-1', 
           'ONECUT1', 'SOX2', 'DLK1', 'ONECUT2', 'RFX6', 'FGFR2', 'TTR', 'CDX2', 'SOX21', 
           'FRZB', 'TOP2A', 'AURKB', 'CPA2', 'NEUROG3', 'DLL1', 'FOXA2', 'BICC2', 'INS', 'NEUROD1',
           'NOTCH2', 'RBPJ', 'FEV', 'GCG', 'GATA4','SOX17', 'HNF4A', 'BICC1']


#%%

# =============================================================================
# SI QUERÉS HACER TODOS LOS GRÁFICOS JUNTOS, CORRER LO SIGUIENTE:

# nombres de todas las tablas
tablas = [('results_D6het_D6wt_Down_norm.csv', 'results_D6hom_D6wt_Down_norm.csv', 'Downregulated genes - Day 6'),
          ('results_D6het_D6wt_Up_norm.csv', 'results_D6hom_D6wt_Up_norm.csv', 'Upregulated genes - Day 6'),
          ('results_D8het_D8wt_Down_norm.csv', 'results_D8hom_D8wt_Down_norm.csv', 'Downregulated genes - Day 8'),
          ('results_D8het_D8wt_Up_norm.csv', 'results_D8hom_D8wt_Up_norm.csv', 'Upregulated genes - Day 8'),
          ('results_D13het_D13wt_Down_norm.csv', 'results_D13hom_D13wt_Down_norm.csv', 'Downregulated genes - Day 13'), 
          ('results_D13het_D13wt_Up_norm.csv', 'results_D13hom_D13wt_Up_norm.csv', 'Upregulated genes - Day 13'),
          ('results_D16het_D16wt_Down_norm.csv', 'results_D16hom_D16wt_Down_norm.csv', 'Downregulated genes - Day 16'),
          ('results_D16het_D16wt_Up_norm.csv', 'results_D16hom_D16wt_Up_norm.csv', 'Upregulated genes - Day 16'),
          ('results_D27het_D27wt_Down_norm.csv', 'results_D27hom_D27wt_Down_norm.csv', 'Downregulated genes - Day 27'), 
          ('results_D27het_D27wt_Up_norm.csv', 'results_D27hom_D27wt_Up_norm.csv', 'Upregulated genes - Day 27'),]


for i in range(len(tablas)):
    
    het = pd.read_csv('/home/eve/Dropbox/doctorado/hnf1b/Table for figure 4 HNF1B paper/'+tablas[i][0], sep=',')
    hom = pd.read_csv('/home/eve/Dropbox/doctorado/hnf1b/Table for figure 4 HNF1B paper/'+tablas[i][1], sep=',')
       
    venn(het, hom, tablas[i][2])
# =============================================================================

#%%


# =============================================================================
# SI QUERÉS HACER SÓLO ALGUNO, CORRER LO SIGUIENTE:
# por ejemplo para day6, Downregulated:
het = pd.read_csv('/home/eve/Dropbox/doctorado/hnf1b/bulk_rna-seq/Table for figure 4 HNF1B paper/results_D13het_D13wt_Down_norm.csv', sep=',')
hom = pd.read_csv('/home/eve/Dropbox/doctorado/hnf1b/bulk_rna-seq/Table for figure 4 HNF1B paper/results_D13hom_D13wt_Down_norm.csv', sep=',')
venn(het, hom, 'Downregulated genes - Day 13')

het = pd.read_csv('/home/eve/Dropbox/doctorado/hnf1b/bulk_rna-seq/Table for figure 4 HNF1B paper/results_D13het_D13wt_Up_norm.csv', sep=',')
hom = pd.read_csv('/home/eve/Dropbox/doctorado/hnf1b/bulk_rna-seq/Table for figure 4 HNF1B paper/results_D13hom_D13wt_Up_norm.csv', sep=',')
venn(het, hom, 'Upregulated genes - Day 13')

# =============================================================================

#%%

trust = pd.read_csv('/home/eve/Dropbox/doctorado/hnf1b/bulk_rna-seq/plots/trrust_rawdata.human.csv', sep='\t') 

TFs = trust.TF

TFs_uniques = [] 
for i in TFs: 
    if i not in TFs_uniques: 
        TFs_uniques.append(i) 

TFs_uniques
len(TFs_uniques)


for i in markers:
    print(i, i in TFs_uniques)



TFs_het = []
for i in genes_het:
    if i in TFs_uniques:
        TFs_het.append(i)
TFs_het
for i in TFs_het:
    print(i)    


TFs_hom = []
for i in genes_hom:
    if i in TFs_uniques:
        TFs_hom.append(i)
TFs_hom    
for i in TFs_hom:
    print(i)    



TFs_both = []
for i in genes_both:
    if i in TFs_uniques:
        TFs_both.append(i)
TFs_both    
for i in TFs_both:
    print(i)  















