import os
import gzip
import pandas as pd
import subprocess
from importlib import resources

class Array():
    def __init__(self):
        self.HLA = ['HLA-A', 'HLA-B', 'HLA-C', 'HLA-DPA1', 'HLA-DPB1', 'HLA-DQA1', 'HLA-DQB1', 'HLA-DRB1']

    def run_snp2hla(self, in_file='1958BC', ref_file='HM_CEU_REF', snp2hla_dir=None, heap_size=2000, window_size=1000):
        if not snp2hla_dir:
            snp2hla_dir = f'{resources.files("hla6").parent.parent}/vendor/SNP2HLA/home'

        out_file = os.path.abspath(f'{in_file}_{ref_file}')

        if os.path.exists(f'{in_file}.bed'):
            in_file = os.path.abspath(in_file)
        else:
            raise FileNotFoundError(f'Input file {in_file}.bed not found in current directory')

        if os.path.exists(f'{ref_file}.bed'):
            ref_file = os.path.abspath(ref_file)
        else:
            raise FileNotFoundError(f'Reference file {ref_file}.bed not found in current directory')

        os.chdir(snp2hla_dir)
        cmd = f'tcsh SNP2HLA.csh {in_file} {ref_file} {out_file} plink {heap_size} {window_size}'
        subprocess.run(cmd, shell=True, check=True)

    def run_deephla(self, mode='train', in_file='1958BC_Pan-Asian_REF', ref_file='Pan-Asian_REF', subset=[], model_json=None, model_dir='model', deephla_dir=None):
        if not deephla_dir:
            deephla_dir = f'{resources.files("hla6").parent.parent}/vendor/DEEP-HLA'

        if mode == 'train':
            if not os.path.exists(f'{in_file}.bgl.phased') or not os.path.exists(f'{ref_file}.bgl.phased'):
                print('Please use hla6 run-snp2hla first to get the phased bgl files')
                return

            if subset:
                region = subset.replace('-', ':').split(':')
                df = pd.read_table(f'{in_file}.bgl.phased', sep=' ', header=None)
                wh = []
                for n in range(df.shape[0]):
                    chrom = str(df.iloc[n, 0])
                    pos = int(df.iloc[n, 1])
                    if chrom == region[0] and pos >= int(region[1]) and pos <= int(region[2]):
                        wh.append(True)
                    else:
                        wh.append(False)
                df = df.loc[wh, ] 
                df.to_csv(f'{in_file}.bgl.phased', header=False, index=False, sep=' ')

            hla_json = f'{ref_file}.hla.json'
            if not os.path.exists(hla_json):
                print('Generating HLA info JSON file...')
                cmd = f'conda run -n DEEP-HLA python {deephla_dir}/make_hlainfo.py --ref {ref_file} --out {ref_file}.hla.json'
                subprocess.run(cmd, shell=True, check=True)

            if not model_json:
                model_json = f'{ref_file}.model.json'
            if os.path.exists(model_json):
                cmd = f'conda run -n DEEP-HLA python {deephla_dir}/train.py --ref {ref_file} --sample {in_file} --model {model_json.replace(".model.json", "")} --hla {hla_json.replace(".hla.json", "")} --model-dir {model_dir}'
                print(cmd)
                subprocess.run(cmd, shell=True, check=True)
            else:
                raise FileNotFoundError(f'{model_json} not found')

        elif mode == 'impute':
            pass

    def format_output(self, in_file='data/1958BC_Euro.bgl.phased', out_file='data/1958BC_Euro_digit4.txt', digit=4, in_type='snp2hla'):
        if in_type == 'snp2hla':
            if in_file.endswith('.phased'):
                header = ['SampleID', 'HLA', 'Allele1', 'Allele2']
                df = pd.read_table(in_file, sep=' ', skiprows=1)
                df = df.loc[df.iloc[:, 1].str.startswith('HLA'), ]
    
                wh = [len(x.split('_')[2]) == digit for x in df.iloc[:, 1]]
                df = df.loc[wh, ]
    
                L = []
                for n in range(2, df.shape[1], 2):
                    sample_id = df.columns[n]
                    allele1 = {}
                    allele2 = {}
                    for m in range(df.shape[0]):
                        allele = df.iloc[m, 1]
                        if digit == 4:
                            allele = f'{allele[0:-2]}:{allele[-2:]}'
                        fields = allele.split('_')
                        k = '-'.join(fields[0:2])
                        if df.iloc[m, n] == 'P':
                            allele1.setdefault(k, [])
                            allele1[k].append(':'.join([k] + fields[2:]))
                        if df.iloc[m, n + 1] == 'P':
                            allele2.setdefault(k, [])
                            allele2[k].append(':'.join([k] + fields[2:]))
    
                    for hla in self.HLA:
                        L.append([sample_id, hla, ','.join(allele1.get(hla, 'X')), ','.join(allele2.get(hla, 'X'))])
    
                df = pd.DataFrame(L)
                df.columns = header
                df.to_csv(out_file, header=True, index=False, sep='\t')
                print('Formatted SNP2HLA output saved to', out_file)

            elif in_file.endswith('.dosage'):
                df = pd.read_table(in_file.replace('.dosage', '.bgl.phased'), sep=' ', skiprows=1)
                sample_name_list = [df.columns[n] for n in range(2, df.shape[1], 2)]
                header = ['SampleID', 'HLA', 'Allele1', 'Allele2', 'Dosage1', 'Dosage2']
                df = pd.read_table(in_file, sep='\t', header=None)
                df = df.loc[df.iloc[:, 0].str.startswith('HLA'), ]
        
                wh = [len(x.split('_')[2]) == digit for x in df.iloc[:, 0]]
                df = df.loc[wh, ]
        
                L = []
                S = []
                for n in range(3, df.shape[1]):
                    D = {}
                    sample_id = sample_name_list[n - 3]
                    if sample_id not in S:
                        S.append(sample_id)
                    for m in range(df.shape[0]):
                        allele = df.iloc[m, 0]
                        if digit == 4:
                            allele = f'{allele[0:-2]}:{allele[-2:]}'
                        fields = allele.split('_')
                        k = '-'.join(fields[0:2])
                        D.setdefault(k, [])
                        D[k].append([':'.join([k] + fields[2:]), df.iloc[m, n]])
        
                    for hla in self.HLA:
                        allele1 = 'X'
                        allele2 = 'X'
                        dosage1 = 'X'
                        dosage2 = 'X'
                        if hla in D:
                            alleles = sorted(D[hla], key=lambda x:float(x[1]), reverse=True)
                            if len(alleles) >= 2:
                                if (float(alleles[0][1]) <= 1.5 and float(alleles[0][1]) > 0.5) \
                                        and (float(alleles[1][1]) <= 1.5 and float(alleles[1][1]) > 0.5):
                                    allele1 = alleles[0][0]
                                    dosage1 = alleles[0][1]
                                    allele2 = alleles[1][0]
                                    dosage2 = alleles[1][1]
                                elif float(alleles[0][1]) > 1.5:
                                    allele1 = alleles[0][0]
                                    dosage1 = alleles[0][1]
                                    allele2 = alleles[0][0]
                                    dosage2 = alleles[0][1]
                            elif len(alleles) == 1:
                                if float(alleles[0][1]) > 1.5:
                                    allele1 = alleles[0][0]
                                    dosage1 = alleles[0][1]
                                    allele2 = alleles[0][0]
                                    dosage2 = alleles[0][1]
                                elif float(alleles[0][1]) <= 1.5 and float(alleles[0][1]) > 0.5:
                                    allele1 = alleles[0][0]
                                    dosage1 = alleles[0][1]
                        L.append([sample_id, hla, allele1, allele2, dosage1, dosage2])
                df = pd.DataFrame(L)
                df.columns = header
                df.to_csv(out_file, header=True, index=False, sep='\t')
                print('Formatted SNP2HLA output saved to', out_file)

if __name__ == "__main__":
    ar = Array()
    ar.run_snp2hla()
