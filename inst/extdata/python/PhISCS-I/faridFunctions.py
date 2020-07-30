import pandas as pd
from collections import defaultdict

def give_me_muts_to_filter(f_i, f_o, kmax):
    report = pd.DataFrame(columns=['mut','a','b','c'])
    def get_decendent(df_o, gene1, gene2):
        df1 = df_o[[gene1]]
        df3 = df1.rename(columns={gene1: "X1"})
        df2 = df_o[[gene2]]
        df4 = df2.rename(columns={gene2: "X2"})
        df = pd.concat([df3, df4], axis=1, join='outer')
        b = df.loc[df['X1'] >= df['X2']].shape[0]
        return b == df.shape[0]
    
    df_i = pd.read_csv(f_i, index_col=0, sep='\t')
    df_o = pd.read_csv(f_o, index_col=0, sep='\t')
    
    subset = defaultdict(list)
    for mut1 in df_o.columns:
        for mut2 in df_o.columns:
            if mut1 != mut2:
                if get_decendent(df_o, mut1, mut2):
                    subset[mut1].append(mut2)
    for k, v in subset.items():
        original_not_present_but_of_it_dcendent_present = 0
        original_present = 0
        for index, row in df_i.iterrows():
            if row[k] == 1:
                original_present += 1
            if row[k] == 0:
                founded = False
                for mut in v:
                    if row[mut] == 1:
                        founded = True
                if founded:
                    original_not_present_but_of_it_dcendent_present += 1
        p = '{:.3f}'.format(1.0*original_not_present_but_of_it_dcendent_present/original_present)
        report.loc[len(report)] = [k,original_not_present_but_of_it_dcendent_present,original_present,p]
    
    report = report.sort_values(by=['c','a'], ascending=[False,False])
    report = report[report.a != 0]
    return [int(item.replace('mut','')) for item in list(report['mut'])[:kmax*3]]