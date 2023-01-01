from math import sqrt
import numpy as np
import pandas as pd
from scipy.stats import pearsonr

deleted_tissue = ['acute_myeloid_leukaemia', 'anaplastic_large_cell_lymphoma', 'B_cell_leukemia', 'B_cell_lymphoma',
                  'bone_other', 'Burkitt_lymphoma', 'chondrosarcoma', 'chronic_myeloid_leukaemia', 'ewings_sarcoma',
                  'fibrosarcoma', 'haematopoietic_neoplasm_other', 'hairy_cell_leukaemia', 'Hodgkin_lymphoma',
                  'leukemia', 'lymphoblastic_leukemia', 'lymphoblastic_T_cell_leukaemia', 'lymphoid_neoplasm_other',
                  'myeloma', 'osteosarcoma', 'rhabdomyosarcoma', 'soft_tissue_other', 'T_cell_leukemia', 'mesothelioma',
                  'lung_NSCLC_not_specified']

zscore_list = []

### open files ###
df1 = pd.read_excel('Merged_datasets_final_30.1.22.xlsx', sheet_name=None)
gene_df = pd.read_excel('CCLE_expression_with_names.xlsx', sheet_name='Sheet1')
drug_df = pd.read_excel('Drug_listWed Jul 28 09_48_04 2021_norep.xlsx', sheet_name='Drug_listWed Jul 28 09_48_04 20',
                        index_col='Name')

for gene in list(gene_df.columns[1:]):
    final_df = pd.DataFrame(
        columns=['glioma', 'lung_small_cell_carcinoma', 'neuroblastoma', 'large_intestine', 'breast', 'melanoma',
                 'bladder', 'cervix', 'lung_NSCLC_large_cell', 'lung_NSCLC_squamous_cell_carcinoma',
                 'lung_NSCLC_adenocarcinoma', 'pancreas', 'oesophagus', 'head_and_neck', 'kidney', 'ovary', 'thyroid',
                 'stomach', 'endometrium', 'liver'])

    final_df['Drug'] = list(df1.keys())
    final_df.set_index('Drug', inplace=True)

    ### open dicts for mRNA ###
    try:
        df2 = gene_df[[gene, 'Cell_line_name']]

        ### open dicts for mRNA ###
        dict1 = df2.set_index('Cell_line_name')[gene].to_dict()

        ### loop over master excel file with sheets and extract data ###
        for sheet in df1.keys():  # loop over drugs (sheets
            new_dict = {"Cell_line": [], "Tissue": [], "Sub_Tissue": [], "AUC": [],
                        gene + "_mRNA": []}  # empty dict with multipule columns

            cells_only_letters = []  # remove '-' from cellines
            for i in df1[sheet]['Cell line name']:
                i = str(i)
                new_name = i.replace('-', '')
                cells_only_letters.append(new_name)
            df1[sheet]['Cell_line_name'] = cells_only_letters
            len_col = len(df1[sheet]['Cell_line_name'])
            for i in range(len_col):  # loop over columns (cell_line)
                if df1[sheet]['Cell_line_name'][i] in dict1.keys():  # and df1[sheet]['Cell_line'][i] in dict2.keys()
                    new_dict['Cell_line'].append(df1[sheet]['Cell_line_name'][i])
                    new_dict['Tissue'].append(df1[sheet]['Tissue'][i])
                    new_dict['Sub_Tissue'].append(df1[sheet]['Tissue sub-type'][i])
                    new_dict['AUC'].append(df1[sheet]['AUC'][i])
                    new_dict[gene + "_mRNA"].append(dict1[df1[sheet]['Cell_line_name'][i]])
            pos_df = pd.DataFrame.from_dict(new_dict)  # convert dict into DF

            ## building colleration curve ###
            tissue_list = pos_df['Sub_Tissue'].unique()
            for tissue in tissue_list:
                if tissue not in deleted_tissue:
                    filter_df = pos_df[pos_df.Sub_Tissue == tissue]
                    x = filter_df[gene + '_mRNA']
                    y = filter_df['AUC']
                    if len(x) and len(y) >= 10:
                        pearson_coef, p_value = pearsonr(x, y)
                        z_score = np.arctanh(pearson_coef)
                        final_z_score = z_score * sqrt(len(x) - 3)
                        final_df.loc[sheet, tissue] = final_z_score
                        zscore_list.append(final_z_score)

        final_df = final_df.fillna(0)

        ### add target payhway for each drug ###
        final_df['Target'] = ""  # add empty columns to dataframe
        final_df['Target_pathway'] = ""
        for drug in final_df.index:
            target = str(drug_df._get_value(drug, 'Pooled Targets'))
            pathway = str(drug_df._get_value(drug, 'Target pathway'))

            final_df.loc[drug, 'Target'] = target
            final_df.loc[drug, 'Target_pathway'] = pathway
        final_df.to_excel(gene + '_PanCanMatrix_sub_Tissue_10_Zscore_1.1.23' + '.xlsx')

    except:
        print(gene + ' has no values')

### for Zscore distripution ###
# zscore_arr = np.array(zscore_list)
# ax1 = sns.displot(zscore_arr, kde=True)
# plt.suptitle(str(np.percentile(zscore_arr, 95)) + ' ' + str(np.percentile(zscore_arr, 5)))
# plt.savefig('zscore_dist.pdf')
# plt.show()

print('Done')
