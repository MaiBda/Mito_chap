import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

number_pathways = {'kinases': 56, 'DNA replication': 28, 'RTK signaling': 46,
                   'PI3K/MTOR signaling': 47, 'Metabolism': 10,
                   'ERK MAPK signaling': 21, 'Genome integrity': 20, 'Mitosis': 19, 'Apoptosis regulation': 22,
                   'Chromatin other': 11, 'Chromatin histone acetylation': 17, 'Cell cycle': 27,
                   'Protein stability and degradation': 10, 'Cytoskeleton': 10,
                   'ErbBs signaling ': 16, 'WNT signaling': 12}

new_indx_dict = {'DNA replication': 'DNA rep.', 'RTK signaling': 'RTK sign.', 'PI3K/MTOR signaling': 'PI3K/mTOR sign.',
                 'ERK MAPK signaling': 'ERK/MAPK sign.', 'Apoptosis regulation': 'Apoptosis reg.',
                 'Chromatin histone acetylation': 'Chromatin acet.',
                 'Protein stability and degradation': 'Protein stab.&deg.',
                 'ErbBs signaling ': 'ErbBs sign.', 'WNT signaling': 'WNT sign.', 'kinases': 'Kinases'}

gene_list = ['ERBB2', 'HSPD1', 'CLPP', 'YME1L1', 'SPG7', 'HTRA2', 'LONP1', 'CLPX', 'HSPA9', 'HSPE1',
             'DNAJA3', 'TRAP1', 'GRPEL2', 'HSCB', 'DNAJC19', 'AFG3L2']

pos_list = []
neg_list = []
pathway_score_list = []

for gene in gene_list:

    df = pd.read_excel(gene + '_PanCanMatrix_sub_Tissue_10_Zscore' + '.xlsx')

    pathway_list = df.Target_pathway.unique()
    pos_neg_df = pd.DataFrame(columns=['pos', 'neg'], index=number_pathways.keys())
    df.rename(columns={'Unnamed: 0': 'Drug'}, inplace=True)
    df = df.set_index('Drug')
    for pathway in pathway_list:
        if pathway in number_pathways.keys():
            pos = []
            neg = []
            filter_df = df[df['Target_pathway'] == pathway]
            filter_df = filter_df.loc[:, (filter_df != 0).any(axis=0)]  # remove columns with zeros only
            filter_df = filter_df[(filter_df.T != 0).any()]  # remove raws with zeros only
            filter_df = filter_df.iloc[:, :-2]

            for col in filter_df.columns:
                for i in filter_df[col]:
                    if i >= 1.7:
                        pos.append(i)
                    if i <= -1.7:
                        neg.append(i)
            if len(pos) > 0 and len(filter_df.columns) > 0:
                pos_neg_df.loc[pathway, 'pos'] = len(pos) / (number_pathways[pathway] * len(filter_df.columns)) * 100
            if len(neg) > 0 and len(filter_df.columns) > 0:
                pos_neg_df.loc[pathway, 'neg'] = len(neg) / (number_pathways[pathway] * len(filter_df.columns)) * 100

    pos_neg_df = pos_neg_df.fillna(0)
    pos_neg_df = pos_neg_df[(pos_neg_df.T != 0).any()]  # remove raws with zeros only
    pos_neg_df = pos_neg_df.rename(index=lambda s: new_indx_dict[s] if s in new_indx_dict.keys() else s)
    for j in pos_neg_df.pos:
        if j > 0:
            pos_list.append(j)

    for m in pos_neg_df.neg:
        if m > 0:
            neg_list.append(m)

    y = pos_neg_df.index
    x1 = pos_neg_df.pos
    x2 = pos_neg_df.neg

    plt.title(gene, size=10)
    bar1 = plt.barh(y, x1, color='indianred')
    bar2 = plt.barh(y, -x2, color='skyblue')
    plt.legend(handles=[mpatches.Patch(color='indianred', label='% pos correlations'),
                        mpatches.Patch(color='skyblue', label='% neg correlations')], loc='upper left', fontsize=6)
    plt.axvline(0, color="black")
    plt.axvline(7.9, color="black", linestyle="dashed")
    plt.axvline(-8.4, color="black", linestyle="dashed")
    plt.savefig(gene + '_biplot_zscore.png')
    plt.tight_layout()
    plt.show()

print('done')
