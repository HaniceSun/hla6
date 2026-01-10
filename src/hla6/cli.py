import argparse
from importlib import resources
from .array import Array
from .seq import Seq

def get_parser():
    formatter_class = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=formatter_class)
    subparsers = parser.add_subparsers(dest='command', required=True)

    p1 = subparsers.add_parser("run-snp2hla", help="run SNP2HLA on array data")
    p1.add_argument('--input', type=str, default='1958BC', help='input file prefix')
    p1.add_argument('--ref', type=str, default='HM_CEU_REF', help='reference panel prefix, can be HM_CEU_REF or Pan-Asian_REF currently')

    p2 = subparsers.add_parser("format-output", help="format the output to allele table")
    p2.add_argument('--input', type=str, default='data/1958BC_Euro.bgl.phased', help='input file')
    p2.add_argument('--input_type', type=str, default='snp2hla', help='where does the input come from')
    p2.add_argument('--output', type=str, default='data/1958BC_Euro_digit4.txt', help='output file')
    p2.add_argument('--digit', type=int, default=4, help='digit level for HLA alleles')

    p3 = subparsers.add_parser("run-deephla", help="run CNN-based DEEP*HLA, to be implemented")
    p3.add_argument('--mode', type=str, default='train', help='mode: train or impute')
    p3.add_argument('--input', type=str, default='1958BC_Pan-Asian_REF', help='input file prefix')
    p3.add_argument('--output', type=str, default='1958BC_Pan-Asian_REF_DEEPHLA', help='output file prefix')
    p3.add_argument('--ref', type=str, default='Pan-Asian_REF', help='reference panel prefix, can be HM_CEU_REF or Pan-Asian_REF currently')
    p3.add_argument('--subset', type=str, default=None, help='subset the input to the HLA regions according to the reference genome, e.g., chr6:28510120-33480577 on GRCh37')
    p3.add_argument('--model_json', type=str, default='Pan-Asian_REF.model.json', help='the config file of the model')
    p3.add_argument('--model_dir', type=str, default='model', help='the output directory of the trained model')
 
    p4 = subparsers.add_parser("run-hlarimnt", help="run Transformer-based HLARIMNT, to be implemented")

    p5 = subparsers.add_parser("run-xHLA", help="run xHLA on sequencing data")
    p5.add_argument('--input', type=str, default='test.bam', help='input bam file')
    p5.add_argument('--output_dir', type=str, default='test', help='output directory')
    p5.add_argument('--sample_id', type=str, default='test', help='sample id used in output files')
    p5.add_argument('--singularity_sif', type=str, default=None, help='path to xHLA singularity image file, using the default one if not provided')
    return parser

def main():
    parser = get_parser()
    args = parser.parse_args()
    ar = Array()
    if args.command == 'run-snp2hla':
        ar.run_snp2hla(in_file=args.input, ref_file=args.ref)
    if args.command == 'run-deephla':
        ar.run_deephla(mode=args.mode, in_file=args.input, ref_file=args.ref, subset=args.subset,
                       model_json=args.model_json, model_dir=args.model_dir, out_file=args.output)
    elif args.command == 'format-output':
        ar.format_output(in_file=args.input, in_type=args.input_type, out_file=args.output, digit=args.digit)
    elif args.command == 'run-xHLA':
        seq = Seq()
        seq.run_xHLA(input_bam=args.input, out_dir=args.output_dir, sample_id=args.sample_id, singularity_sif=args.singularity_sif)

if __name__ == '__main__':
    main()
