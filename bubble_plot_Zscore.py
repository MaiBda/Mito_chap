import pandas as pd
import matplotlib.pyplot as plt

number_pathways = {'kinases': 56, 'DNA replication': 28, 'RTK signaling': 46,
                   'PI3K/MTOR signaling': 47, 'Metabolism': 10,
                   'ERK MAPK signaling': 21, 'Genome integrity': 20, 'Mitosis': 19, 'Apoptosis regulation': 22,
                   'Chromatin other': 11, 'Chromatin histone acetylation': 17, 'Cell cycle': 27,
                   'Protein stability and degradation': 10, 'Cytoskeleton': 10,
                   'ErbBs signaling ': 16, 'WNT signaling': 12}

new_col_dict = {'lung_small_cell_carcinoma': 'Lung_SCLC', 'lung_NSCLC_large_cell': 'Lung_Nsclc_LC',
                'lung_NSCLC_squamous_cell_carcinoma': 'Lung_Nsclc_SCC', 'lung_NSCLC_adenocarcinoma': 'Lung_Nsclc_AC',
                'head_and_neck': 'H&N', 'lung_NSCLC_not_specified': 'Lung_Nsclc_NS'}

new_indx_list = ['Kinases', 'DNA rep.', 'RTK sign.', 'PI3K/mTOR sign.', 'Metabolism',
                 'ERK MAPK sign.', 'Genome integrity', 'Mitosis', 'Apoptosis reg.', 'Chromatin other',
                 'Chromatin acet.', 'Cell cycle', 'Protein stab.&deg.', 'Cytoskeleton',
                 'ErbBs sign.', 'WNT sign.']

gene_list = ['CLPP']

indx_list = []

col_list = []
score_list = []
for gene in gene_list:
    df = pd.read_excel(gene + '_PanCanMatrix_sub_Tissue_10_Zscore.xlsx')
    df.rename(columns={'Unnamed: 0': 'Drug'}, inplace=True)
    pos_df = pd.DataFrame(index=number_pathways.keys(), columns=df.columns[1:-2])
    neg_df = pd.DataFrame(index=number_pathways.keys(), columns=df.columns[1:-2])
    for pathway in number_pathways.keys():
        ERBB2_df = df[df.Target_pathway == pathway]
        for col in ERBB2_df.columns[1:-2]:
            pos = []
            neg = []
            for i in ERBB2_df[col]:
                if i >= 1.7:
                    pos.append(i)
                if i <= -1.7:
                    neg.append(i)

            pos_df.loc[pathway, col] = round(len(pos) / number_pathways[pathway], 2) * 100
            neg_df.loc[pathway, col] = round(len(neg) / number_pathways[pathway], 2) * -100
            if col not in col_list:
                col_list.append(col)

    for datafrm in [pos_df, neg_df]:
        datafrm['new_idx'] = new_indx_list
        datafrm = datafrm.set_index('new_idx')

        for col in datafrm.columns:
            if col in new_col_dict.keys():
                datafrm.rename({col: new_col_dict[col]}, axis=1, inplace=True)
            else:
                datafrm.rename({col: col.capitalize()}, axis=1, inplace=True)
        x_lst = []
        y_lst = []
        S_list = []
        C_list = []
        bubble_df = pd.DataFrame(columns=['X', 'Y', 'S'])
        for col in datafrm.columns:

            for i in range(len(datafrm.index)):
                x_lst.append(col)
                y_lst.append(datafrm.index[i])
                value = datafrm.at[datafrm.index[i], col]
                if value > 0:
                    S_list.append(value * 5)
                    C_list.append('indianred')
                if value < 0:
                    S_list.append(abs(value * 5))
                    C_list.append('skyblue')
                if value == 0:
                    S_list.append(value)
                    C_list.append('white')
        bubble_df['X'] = x_lst
        bubble_df['Y'] = y_lst
        bubble_df['S'] = S_list
        bubble_df['C'] = C_list
        for i in bubble_df.S:
            score_list.append(i / 5)

        dfu = bubble_df
        fig, ax = plt.subplots()
        plt.rcParams['ytick.right'] = plt.rcParams['ytick.labelright'] = True
        plt.rcParams['ytick.left'] = plt.rcParams['ytick.labelleft'] = False
        sc = plt.scatter(x="X", y="Y", s="S", c="C", data=dfu)

        plt.xticks(rotation='vertical')
        for line in range(0, bubble_df.X.shape[0]):
            if bubble_df.S[line] != 0:
                plt.text(bubble_df.X[line], bubble_df.Y[line], int(round(bubble_df.S[line] / 5, 2)),
                         horizontalalignment='center', size=5)

        plt.tight_layout()
        plt.show()
