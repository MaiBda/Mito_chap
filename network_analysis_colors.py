import pandas as pd
import numpy as np
import networkx as nx
import matplotlib.pyplot as plt

gene_df = pd.read_excel('HSPD1_client_10.10.xlsx',sheet_name='PD1_PE1_mutual')

HSPD1 = list(gene_df.only_HSPD1)
HSPE1 = list(gene_df.only_HSPE1)
mutual = list(gene_df.Mutual)
all_HSPD1 = list(gene_df.HSPD1_all)
all_HSPE1 = list(gene_df.HSPE1_all)

all_HSPD1_with_blanks = []
all_HSPE1_with_blanks = []

HSPD1_with_blank = []
HSPE1_with_blank = []
mutual_with_blank = []

for gene in HSPD1:
    gene_with_blank = str(gene) + ' '
    HSPD1_with_blank.append(gene_with_blank)
for gene in HSPE1:
    gene_with_blank = str(gene) + ' '
    HSPE1_with_blank.append(gene_with_blank)
for gene in mutual:
    gene_with_blank = str(gene) + ' '
    mutual_with_blank.append(gene_with_blank)
for gene in all_HSPD1:
    gene_with_blank = str(gene) + ' '
    all_HSPD1_with_blanks.append(gene_with_blank)
for gene in all_HSPE1:
    gene_with_blank = str(gene) + ' '
    all_HSPE1_with_blanks.append(gene_with_blank)
df = pd.read_excel('Matrix_for_network_analysis_sorted_allgenes.xlsx', sheet_name=None)
cancer_type = ['Ovary', 'Lung_Nsclc_Adenocarcinoma']

for cancer in cancer_type:
    new_dict = {}
    for sheet in df.keys():
        # new_dict[sheet] = list(df[sheet].Lung_Nsclc_Adenocarcinoma)
        if sheet in all_HSPD1_with_blanks or sheet in all_HSPE1_with_blanks:
            new_dict[sheet] = list(df[sheet][cancer])

            new_dict['Drug'] = list(df[sheet].Drug)
    new_df = pd.DataFrame(new_dict)
    new_df['Drug'] = list(df['HSPD1 '].Drug)
    new_df = new_df.set_index('Drug')

    # pos_df.to_excel(cancer + "_matrix_sorted_allgenes.xlsx")
    # Calculate the correlation between individuals. We have to transpose first, because the corr function calculate the pairwise correlations between columns.
    corr = new_df.corr()

    # Transform it in a links data frame (3 columns only):
    links = corr.stack().reset_index()
    links.columns = ['var1', 'var2', 'value']
    # Keep only correlation over a threshold and remove self correlation (cor(A,A)=1)
    links_filtered = links.loc[(links['value'] > 0.5) & (links['var1'] != links['var2'])]

    #  Build graph
    G = nx.from_pandas_edgelist(links_filtered, 'var1', 'var2')

    # set node color
    color_map = []
    for node in list(G):
        if node in HSPD1_with_blank:
            color_map.append('blue')
        if node in HSPE1_with_blank:
            color_map.append('green')
        if node in mutual_with_blank:
            color_map.append('red')
    # Plot the network:
    nx.draw_spring(G, node_color=color_map, with_labels=True, node_size=40, edge_color='black', linewidths=0.09,
                   font_size=0.000000000000000000001, width=0.05)
    ax = plt.gca()
    plt.axis("off")

    plt.savefig('Network_analysis_' + cancer + '_PD+PE_colors.pdf')
    plt.show()
