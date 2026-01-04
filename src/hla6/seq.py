import os
import subprocess
from importlib import resources

class Seq():
    def __init__(self):
        pass

    def run_xHLA(self, xHLA_sif=None, input_bam='test.bam', output_dir='test', sample_id='test'):
        if xHLA_sif is None:
            xHLA_sif = f'{resources.files("hla6").parent.parent}/vendor/xHLA.sif'

            if not os.path.isfile(xHLA_sif):
                raise FileNotFoundError(f'xHLA Singularity image not found at {xHLA_sif}')

        cmd = f'singularity exec -B `pwd`:`pwd` --pwd `pwd` {xHLA_sif} run.py --sample_id {sample_id} --input_bam_path {input_bam} --output_path {output_dir}'
        subprocess.run(cmd, shell=True, check=True)

