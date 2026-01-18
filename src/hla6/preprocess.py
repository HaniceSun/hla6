import numpy as np
import pandas as pd

class DataPreprocessor:
    def __init__(self, ref_bim='Pan-Asian_REF.bim', sample_bim='1958BC.bim', ref_phased='Pan-Asian_REF.bgl.phased'):
        self.ref_bim = pd.read_table(ref_bim, header=None, sep='\t')
        self.ref_phased = pd.read_table(ref_phased, header=None, sep=' ')
        self.sample_bim = pd.read_table(sample_bim, header=None, sep='\t')

        bim_cols = ['chrom', 'id', 'cm', 'pos', 'A1', 'A2']
        self.ref_bim.columns = [f'{x}_ref' for x in bim_cols]
        self.sample_bim.columns = [f'{x}_sample' for x in bim_cols]

        cols = list(self.ref_phased.columns)
        cols[0] = 'I'
        cols[1] = 'id_ref'
        for i in range(2, self.ref_phased.shape[1], 2):
            cols[i] = f'A1_s{i//2}'
            cols[i + 1] = f'A2_s{i//2}'
        self.ref_phased.columns = cols

        self.hla_filter = ['HLA']
        wh = []
        for n in range(self.ref_phased.shape[0]):
            id_ref = self.ref_phased['id_ref'].iloc[n]
            flag = False
            for hf in self.hla_filter:
                if str(id_ref).find(hf) != -1:
                    flag = True
                    break
            wh.append(flag)
        wh = np.array(wh)
        self.ref_phased_hla = self.ref_phased[wh]
        self.ref_phased_non_hla = self.ref_phased[~wh]

    def make_features(self, by_id=True, by_pos=False, check_alleles=True, subset_sample_bim=None, out_file='features.txt'):
        # Subset sample bim to the HLA region
        if subset_sample_bim is not None:
            wh = []
            for n in range(self.sample_bim.shape[0]):
                flag = False
                ch = self.sample_bim['chrom_sample'].iloc[n]
                pos = self.sample_bim['pos_sample'].iloc[n]
                if str(ch) == str(subset_sample_bim[0]):
                    if pos < subset_sample_bim[2] and pos > subset_sample_bim[1]:
                        flag = True
                        break
                wh.append(flag)
            self.sample_bim = self.sample_bim[wh]
            print(f"Subsetted sample bim to the region: {subset_sample_bim}")

        # Get shared variants between reference and sample
        if by_id:
            df = pd.merge(self.ref_bim, self.sample_bim, left_on='id_ref', right_on='id_sample')
        elif by_pos:
            df = pd.merge(self.ref_bim, self.sample_bim, left_on=['chrom_ref', 'pos_ref'], right_on=['chrom_sample', 'pos_sample'])
        else:
            raise ValueError("Either by_id or by_pos must be True.")

        if check_alleles:
            wh = ( ( (df['A1_ref'] == df['A1_sample']) & (df['A2_ref'] == df['A2_sample']) ) | 
                   ( (df['A1_ref'] == df['A2_sample']) & (df['A2_ref'] == df['A1_sample']) ) )
            df = df[wh]

        # Merge with phased data
        df = pd.merge(df, self.ref_phased_non_hla, on='id_ref') 
        L = []
        L.append(df[['id_ref', 'pos_ref']])
        for col in df.columns:
            if str(col).startswith('A') and str(col).find('_ref') == -1 and str(col).find('_sample') == -1:
                d = {col:(df[col]==df['A1_ref']).astype(int)}
                L.append( pd.DataFrame(d))
        df = pd.concat(L, axis=1)
        df = df.sort_values(by='pos_ref').reset_index(drop=True)
        print(df)

        M = np.zeros((int(df.shape[1]/2 - 1), df.shape[0] * 2), dtype=int)
        H = []
        S = []
        for n in range(2, df.shape[1], 2):
            S.append(df.columns[n].split('_')[-1])
            for m in range(df.shape[0]):
                for j in range(2):
                    M[n // 2 - 1, m * 2 + j] = df.iloc[m, n + j]
                if n == 2:
                    H.append(f'A1_v{m+1}')
                    H.append(f'A2_v{m+1}')
        df = pd.DataFrame(M)
        df.index = S
        df.index.name = 'sample'
        df.columns = H
        df.to_csv(out_file, sep='\t', index=True, header=True)

    def make_labels(self, out_file='labels.txt', encoding_file='label_encoding.txt'):
        self.get_label_encoding(encoding_file)
        print(self.ref_phased_hla)
        n_heads = (max(self.encoding['head_idx'].unique()) + 1) * 2
        heads = self.encoding.sort_values(by='head_idx')['head'].unique()
        S = []
        H = []
        for h in heads:
            H.append(f"A1_{h}")
            H.append(f"A2_{h}")
        M = np.zeros((int(self.ref_phased_hla.shape[1]/2 - 1), n_heads), dtype=int)
        for n in range(2, self.ref_phased_hla.shape[1], 2):
            S.append(self.ref_phased_hla.columns[n].split('_')[-1])
            for j in range(2):
                variants = self.ref_phased_hla.loc[self.ref_phased_hla.iloc[:, n + j] == 'P', 'id_ref'].values
                df_sub = self.encoding[self.encoding['allele'].isin(variants)]
                for m in range(df_sub.shape[0]):
                    head_idx = df_sub['head_idx'].iloc[m]
                    label = df_sub['label'].iloc[m]
                    M[n // 2 - 1, head_idx * 2 + j] = label
        df = pd.DataFrame(M)
        df.index = S
        df.index.name = 'sample'
        df.columns = H
        df.to_csv(out_file, sep='\t', index=True, header=True)

    def get_label_encoding(self, encoding_file='label_encoding.txt'):
        D = {}
        for n in range(self.ref_phased_hla.shape[0]):
            id_ref = self.ref_phased_hla['id_ref'].iloc[n]
            head = '_'.join(id_ref.split('_')[0:2])
            dl = len(id_ref.split('_')[-1])
            if dl >= 4:
                head = id_ref[0:-2]
            D.setdefault(head, [])
            if id_ref not in D[head]:
                D[head].append(id_ref)

        H = {}
        for k in sorted(D):
            H[k] = sorted(D[k])

        encoding = []
        for head in H:
            for allele in H[head]:
                digit = len(allele.split('_')[-1])//2*2
                encoding.append([digit, allele, H[head].index(allele) + 1, head])
        encoding = pd.DataFrame(encoding, columns=['digit', 'allele', 'label', 'head'])
        encoding.sort_values(by=['digit', 'head'], inplace=True)

        heads = []
        for digit in sorted(encoding['digit'].unique()):
            df_sub = encoding[encoding['digit'] == digit]
            for n in range(df_sub.shape[0]):
                head = df_sub['head'].iloc[n]
                if head not in heads:
                    heads.append(head)
        encoding['head_idx'] = [heads.index(x) for x in encoding['head']]

        parent = []
        parent_idx = []
        parent_value = []
        for n in range(encoding.shape[0]):
            p = '.'
            p_idx = -1
            p_val = -1
            digit = encoding['digit'].iloc[n]
            head = encoding['head'].iloc[n]
            if digit == 4:
                p = '_'.join(head.split('_')[0:2])
            elif digit >= 6:
                p = head[0:-2]
            if p != '.':
                p_idx = heads.index(p)
                p_val = H[p].index(head) + 1
            parent.append(p)
            parent_idx.append(p_idx)
            parent_value.append(p_val)
        encoding['parent'] = parent
        encoding['parent_idx'] = parent_idx
        encoding['parent_val'] = parent_value
        encoding.to_csv(encoding_file, sep='\t', index=False)
        self.encoding = encoding

dp = DataPreprocessor()
dp.make_features()
dp.make_labels()
