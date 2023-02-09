import pandas as pd
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test
from matplotlib.offsetbox import AnchoredText
gene_list = ['CLPP', 'YME1L1', 'SPG7', 'HTRA2', 'LONP1', 'CLPX', 'HSPA9', 'HSPE1', 'HSPD1',
              'DNAJA3', 'TRAP1', 'GRPEL2', 'HSCB', 'DNAJC19', 'AFG3L2']

df1 = pd.read_excel('PanCan_TCGA_30.8.xlsx', sheet_name=None)

for gene in gene_list:
    df2 = pd.read_excel('splited_RNAseq_' + gene + '.xlsx', sheet_name=None)
    for sheet in df2.keys():
        try:
            quantil = df2[sheet][gene].quantile(0.75)
        except:
            quantil = df2[sheet][gene + '_RNAseqV2RSEM'].quantile(0.75)

        try:
            df2[sheet].loc[df2[sheet][gene] >= quantil, 'expression'] = 3  # higher than quantil
            df2[sheet].loc[df2[sheet][gene] < quantil, 'expression'] = 2  # lower than quantil
        except:
            df2[sheet].loc[df2[sheet][gene + '_RNAseqV2RSEM'] >= quantil, 'expression'] = 3  # higher than quantil
            df2[sheet].loc[df2[sheet][gene + '_RNAseqV2RSEM'] < quantil, 'expression'] = 2  # lower than quantil

    for sheet_1 in df1.keys():
        for sheet_2 in df2.keys():
            if sheet_1 == sheet_2:
                for i in range(len(df1[sheet_1]['bcr_patient_barcode'])):
                    if df1[sheet_1]['pharmaceutical_therapy_type'][i] == 'Chemotherapy':
                    # if any(map((lambda value: value == df1[sheet_1]['pharmaceutical_therapy_drug_name'][i]),
                    #            ('5-Fluorouracil','Etoposide','Oxaliplatin','Paclitaxel','Gemcitabine'))):
                        for m in range(len(df1[sheet_1]['Patient ID'])):
                            if df1[sheet_1]['bcr_patient_barcode'][i] == df1[sheet_1]['Patient ID'][m]:
                                for j in df2[sheet_2].Mod_SAMPLE_ID:
                                    if df1[sheet_1]['Patient ID'][m] == j:
                                        df2[sheet_2].loc[df2[sheet_2]['Mod_SAMPLE_ID'] == j, 'survival'] = \
                                            df1[sheet_1]['OS_MONTHS'][m]
                                        df2[sheet_2].loc[df2[sheet_2]['Mod_SAMPLE_ID'] == j, 'status'] = \
                                            str(df1[sheet_1]['OS_STATUS'][m])[0]

        df3 = df2[sheet_1]
        df3 = df3.dropna(subset=["survival"])

        if len(df3['survival']) > 10:
            kmf_high = KaplanMeierFitter()
            kmf_low = KaplanMeierFitter()
            high = df3.query('expression == 3')
            low = df3.query('expression == 2')

            kmf_high.fit(durations=high['survival'], event_observed=high['status'], label=gene + ' high expression')
            kmf_low.fit(durations=low['survival'], event_observed=low['status'], label=gene + ' low expression')
            p = round(logrank_test(durations_A=high['survival'], durations_B=low['survival'],
                                   event_observed_hi=high['status'],
                                   event_observed_lo=low['status'], weightings='wilcoxon').p_value, 3)
            if p <= 0.050:
                kmf_high.plot_survival_function(ci_show=True)
                kmf_low.plot_survival_function(ci_show=True)
                plt.xlabel('Months')
                plt.ylabel('Survival')
                plt.title('p-value=' + str(p) + " " + sheet_1)
                # plt.savefig(gene +'_'+ sheet_1 +'_'+'chemo_survival.pdf')
                plt.show()
