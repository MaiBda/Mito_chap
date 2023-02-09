import numpy as np
import pandas as pd
from matplotlib import pyplot as plt
from sklearn.datasets import _samples_generator as sg
from sklearn.cluster import SpectralBiclustering
from pandas import ExcelWriter
from sklearn.metrics import consensus_score
import numpy as np

drug_df = pd.read_excel('Drug_listWed Jul 28 09_48_04 2021.xlsx', sheet_name='Drug_listWed Jul 28 09_48_04 20')
gene_list = ['CLPP', 'YME1L1', 'SPG7', 'HTRA2', 'LONP1', 'CLPX', 'HSPA9', 'HSPE1', 'HSPD1',
             'DNAJA3', 'TRAP1', 'GRPEL2', 'HSCB', 'DNAJC19', 'AFG3L2']

for gene in gene_list:
    cells_df = pd.read_excel(gene + '_PanCan_CCLE.xlsx')
    AUC_df = pd.read_excel('Alldrugs_PANCAN_' + gene + '_Sub-Tissue_20.2.22.xlsx', sheet_name=None)
    # df = pd.read_excel(gene + '_Score_matrix.xlsx')
    df = pd.read_excel('chap_correlation_matrix.xlsx')
    # df.rename(columns={'Unnamed: 0': 'Cell_line'}, inplace=True)
    df.set_index('Drugs', inplace=True)
    df = df.fillna(0)  # fill empty cells with 0
    df = df.loc[:, (df != 0).any(axis=0)]  # remove columns with zeros only
    df = df[(df.T != 0).any()]  # remove raws with zeros only
    df = df.T
    data = pd.DataFrame(df).to_numpy()

    n_clusters = (3)
    # plt.matshow(data, cmap=plt.cm.seismic)
    # plt.title("Original dataset" + "_" + gene)

    rng = np.random.RandomState(0)
    row_idx = rng.permutation(data.shape[0])
    col_idx = rng.permutation(data.shape[1])
    data = data[row_idx][:, col_idx]
    data, row_idx, col_idx = sg._shuffle(data, random_state=0, )
    # plt.matshow(data, cmap=plt.cm.seismic)
    # plt.title("Shuffled dataset")
    model = SpectralBiclustering(n_clusters=n_clusters, random_state=0)
    model.fit(data)

    row_labels = model.row_labels_
    col_labels = model.column_labels_

    fit_data = data[np.argsort(model.row_labels_)]
    fit_data = fit_data[:, np.argsort(model.column_labels_)]
    # plt.matshow(fit_data, cmap=plt.cm.seismic)
    # plt.title("After biclustering; rearranged to show biclusters")

    plt.matshow(np.outer(np.sort(model.row_labels_) + 1, np.sort(model.column_labels_) + 1), cmap=plt.cm.seismic)
    # plt.title("Checkerboard structure of rearranged data" + "_" + gene)
    plt.tight_layout()
    plt.show()

    # test = model.get_indices(15)
    # print(test)

    col_labels = list(model.column_labels_)
    row_labels = list(model.row_labels_)
    col_clusters = []
    row_clusters = []
    for i in range(0, len(col_labels)):
        col_clusters.append((col_labels[i], df.columns[i]))
    for j in range(0, len(row_labels)):
        row_clusters.append((row_labels[j], df.index[j]))

    col_dict = {}
    row_dict = {}
    for a, b in (col_clusters):
        col_dict[b] = a
    for a, b in (row_clusters):
        row_dict[b] = a
    # print(col_dict)
    col_label_df = pd.DataFrame()
    row_label_df = pd.DataFrame()
    keys_list = [col_dict.keys()]
    value_list = [col_dict.values()]
    keys_list_row = [row_dict.keys()]
    value_list_row = [row_dict.values()]
    col_label_df['Drugs'] = list(col_dict.keys())
    col_label_df['Cluster'] = list(col_dict.values())
    row_label_df['Cell_lines'] = list(row_dict.keys())
    row_label_df['Cell_Cluster'] = list(row_dict.values())
    col_label_df.set_index('Drugs', inplace=True)
    row_label_df.set_index('Cell_lines', inplace=True)
    for i in range(len(drug_df.Name)):
        for drug in col_label_df.index:
            if drug_df.Name[i] == drug:
                col_label_df.at[drug, 'Pathway'] = drug_df['Targets'][i]
                col_label_df.at[drug, 'Target_pathway'] = drug_df['Target pathway'][i]

    for j in range(len(cells_df['Cell Line Name'])):
        for cell in row_label_df.index:
            if cells_df['Cell Line Name'][j] == cell:
                row_label_df.at[cell, 'cell_type'] = cells_df['Primary Disease'][j]
                row_label_df.at[cell, 'mRNA'] = cells_df['Expression 21Q2 Public'][j]

    with ExcelWriter(('chap' + '_' + str(n_clusters) + '_clusters.xlsx')) as writer:
        col_label_df.to_excel(writer, sheet_name='Drugs')
        row_label_df.to_excel(writer, sheet_name='Cancers')



