#!/usr/bin/env python3

"""Entry point"""


def main():

    from argparse import ArgumentParser
    from . import SequenceBouncer


    # Load files, receive parameters, and provide assistance

    ap = ArgumentParser()
    ap.add_argument('-i','--input_file',required=True,type=str,help='Input file in FASTA format.\n')
    ap.add_argument('-o','--output_file',required=False,type=str,default=None,help="Output filestem [do not include extensions] (default will be '<input_filestem>.ext').\n")
    ap.add_argument('-g','--gap_percent_cut',required=False,type=float,default=2.0,help='For columns with a greater fraction of gaps than the selected value, expressed in percent, data will be ignored in calculations (default is 2).\n')
    ap.add_argument('-k','--IQR_coefficient',required=False,type=float,default=1.0,help='Coefficient multiplied by the interquartile range that helps to define an outlier sequence (default is 1.0).\n')
    ap.add_argument('-n','--subsample_size',required=False,type=int,default=0,help='|> Available for large alignments | The size of a single sample taken from the full dataset (default is entire alignment, but try a subsample size of 50 or 100 for large alignments).\n')
    ap.add_argument('-t','--trials',required=False,type=int,default=1,help='|> Available for large alignments | Number of times each sequence is sampled and tested (default is to examine all sequences in one single trial, but 5 or 10 trials may work well when subsamples are taken from large alignments).\n')
    ap.add_argument('-s','--stringency',required=False,type=int,default=2,help='|> Available for large alignments | 1: Minimal stringency 2: Moderate stringency 3: Maximum stringency (default is moderate stringency).\n')
    ap.add_argument('-r','--random_seed',required=False,type=int,default=None,help='Random seed (integer) to be used during a sampling-based approach (default is that the seed is randomly selected). The user can use this seed to obtain reproducible output and should note it in their publications. \n')
    args = vars(ap.parse_args())

    args['input'] = args['input_file']
    del args['input_file']

    SequenceBouncer(**args)()


if __name__=='__main__':
    main()
