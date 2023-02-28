import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

mito_gene_df = pd.read_excel('HSPD1_client_10.10.xlsx', sheet_name='all_mito_genes')
zscore_list = []
gene_df = pd.read_excel('CCLE_expression_with_names.xlsx', sheet_name='Sheet1')

for gene in list(gene_df.columns[1:]):
    try:
        df = pd.read_excel(gene + '_PanCanMatrix_sub_Tissue_10_Zscore' + '.xlsx')
        if not df.empty:
            df.rename(columns={'Unnamed: 0': 'Drug'}, inplace=True)
            df = df.set_index('Drug')
            df = df[df.columns[0:-2]]
            for col in df.columns:
                for i in df[col]:
                    # if i != 0 and 5>=i>=-5:
                    if i != 0:
                        zscore_list.append(i)
    except:
        print(gene + ' not in CCLE')
zscore_arr = np.array(zscore_list)
new_zscore_arr = zscore_arr[np.isfinite(zscore_arr)]

cutoff_range = np.asarray([round(x, 2) for x in np.arange(1, 10 + 0.1, 0.1)])
corr_vc = new_zscore_arr

prop_gn_pairs_up = []
prop_gn_pairs_dw = []
for cutoff in cutoff_range:
    up = len(np.where(corr_vc > cutoff)[0]) / len(corr_vc)
    dw = len(np.where(corr_vc < -cutoff)[0]) / len(corr_vc)
    prop_gn_pairs_up.append(up)
    prop_gn_pairs_dw.append(dw)

prop_gn_pairs_up = np.asarray(prop_gn_pairs_up)
prop_gn_pairs_dw = np.asarray(prop_gn_pairs_dw)

ix_005_r = min(np.where(prop_gn_pairs_up <= 0.05)[0])
ix_005_l = min(np.where(prop_gn_pairs_dw <= 0.05)[0])
zscore_r = round(cutoff_range[ix_005_r], 2)
zscore_l = round(cutoff_range[ix_005_l], 2)
print(zscore_l, zscore_r)

# plot dist#
ax1 = sns.displot(new_zscore_arr, kde=True)
plt.suptitle(str(zscore_r) + '          ' + str(-zscore_l))
plt.axvline(zscore_r, color="black", linestyle="dashed")
plt.axvline(-zscore_l, color="black", linestyle="dashed")
plt.xlabel('Z-score')
plt.savefig('Zscore_dist.pdf')
plt.show()
print('Done')
