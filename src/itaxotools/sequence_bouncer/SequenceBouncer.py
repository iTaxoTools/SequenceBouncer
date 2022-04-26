#!/usr/bin/env python
# coding: utf-8

# Funding received from the Sigrid Jusélius Foundation contributed to the development of this software.
# Author: Cory Dunn
# Institution: University of Helsinki
# Author Email: cory.dunn@helsinki.fi
# Version: 1.23
version = '1.23'
# License: GPLv3

# Modified for iTaxoTools by Stefanos Patmanidis


import time
import sys
import math
import gc
import random
import logging

from typing import Dict, Iterator, Optional, Union
from io import TextIOWrapper
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import rcParams
from matplotlib.backends.backend_pdf import PdfPages
from Bio import AlignIO, SeqIO
from Bio.Align import AlignInfo
from Bio.SeqRecord import SeqRecord


# Logging to console

stream = logging.StreamHandler()
stream.setLevel(logging.INFO)
streamformat = logging.Formatter("%(message)s")
stream.setFormatter(streamformat)

logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
logger.addHandler(stream)


class AttrDict(dict):
    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.__dict__ = self


class InvalidInput(Exception):
    def __init__(self):
        super().__init__(
            'Input must be of type: str, Path, TextIOWrapper or InputSequence.')


class AllColumnsRemovedAsGaps(Exception):
    def __init__(self):
        super().__init__((
            'All columns were removed as gaps. '
            'Choose a larger value for --gap_percent_cut to continue.'
            ))


class InputSequence:
    def __init__(self, name: str, iterator: Iterator[SeqRecord]):
        self.name = name
        self.iterator = iterator

    def __iter__(self):
        return (record for record in self.iterator)

    def __next__(self):
        return next(self.iterator)


class SequenceBouncer():
    """
    Detect and remove outlier sequences from a multiple sequence alignment.
    Initialize with parameters, then call to execute.
    """

    def __init__(
        self,
        input: Union[str, Path, TextIOWrapper],
        output_file: Optional[Union[str, Path]] = None,
        gap_percent_cut: float = 2.0,
        IQR_coefficient: float = 1.0,
        subsample_size: int = 0,
        trials: int = 1,
        stringency: int = 2,
        random_seed: int = None,
        write_log_file: bool = True,
        write_sequence_files: bool = True,
        write_table_files: bool = True,
        write_plot_files: bool = True,
        write_none: bool = False,
        input_format: str = 'fasta',
    ):

        self.vars = AttrDict()
        self.params = AttrDict(
            input=input,
            output_file=output_file,
            gap_percent_cut=gap_percent_cut,
            IQR_coefficient=IQR_coefficient,
            subsample_size=subsample_size,
            trials=trials,
            stringency=stringency,
            random_seed=random_seed,
            write_log_file=write_log_file,
            write_sequence_files=write_sequence_files,
            write_table_files=write_table_files,
            write_plot_files=write_plot_files,
            write_none=write_none,
            input_format=input_format,
        )

    def __call__(self) -> Dict[str, bool]:
        """
        Returns a dictionary where the keys are the sequence names,
        and the values are True if accepted and False for outliers.
        """
        v = self.vars
        p = self.params

        v.start_time = time.time()

        # Initialize

        self.parse_params()

        self.initialize_logger()

        self.parse_input()

        # Prepare dataframe for storage of trial results (these columns are stripped away later if only one trial is performed)

        v.record_sequence_trial_results = pd.DataFrame(v.alignment_record_name_list, columns=['Accession'])
        v.record_sequence_trial_results['Total_trials'] = 0
        v.record_sequence_trial_results['Outlier_instances'] = 0

        # Print input summary

        self.logger.info("Analyzing '" + v.input_name + "'.")
        self.logger.info('Flags are --IQR_coefficient: ' + str(v.multiplier_on_interquartile_range) + ', -subsample_size: ' + str(v.number_in_small_test) + ', --gap_percent_cut: ' + str(v.gap_value_cutoff))
        if v.min_trials_for_each_sequence != 1:
            self.logger.info('          --stringency: ' + str(v.stringency_flag) + ', --trials: ' + str(v.min_trials_for_each_sequence) + ', --random_seed: ' + str(v.seed))

        self.logger.info('Input alignment length is: ' + str(v.length_of_alignment) + ' characters.')
        self.logger.info("Input alignment depth is: " + str(v.depth_of_alignment) + " sequences.")

        # Calculate Shannon entropy and fraction of column gapped

        self.logger.info('Calculating Shannon entropy values and gap metrics across all input sequences.')
        v.entropy_record = []
        v.gap_record = []
        sequence_columns = len(v.sequence_dataframe.axes[1])
        for i in range(sequence_columns):
            column_fractions_S = v.sequence_dataframe[i].value_counts(normalize='True')
            shannon_entropy_column = 0
            gap_fraction_column = 0
            for character, val in  column_fractions_S.iteritems():
                shannon_entropy_column +=  val * math.log(val,2)
                if character == 45:
                    gap_fraction_column = val
            shannon_entropy_column *= -1
            v.entropy_record.append(shannon_entropy_column)
            v.gap_record.append(gap_fraction_column)
        v.entropylist_S = pd.Series(v.entropy_record)
        v.gap_fraction_S = pd.Series(v.gap_record)

        if p.write_plot_files:
            self.write_gap_plot()

        # Generate boolean based upon gap values

        v.gap_value_cutoff_float = float(v.gap_value_cutoff/100)
        gap_percent_bool_series_remove = v.gap_fraction_S > v.gap_value_cutoff_float
        v.gap_percent_bool_index_remove = gap_percent_bool_series_remove[gap_percent_bool_series_remove].index

        # Remove gapped positions

        v.entropylist_S_gap_considered = v.entropylist_S.drop(v.gap_percent_bool_index_remove)

        if v.entropylist_S_gap_considered.size == 0:
            self.logger.info('All columns were removed as gaps.')
            self.logger.info('Choose a larger value for --gap_percent_cut to continue.')
            raise AllColumnsRemovedAsGaps()

        max_entropy_before_gaps = pd.Series.max(v.entropylist_S)
        self.logger.info('Maximum Shannon entropy alignment score before gap % considered: ' + str(round(max_entropy_before_gaps,2)))
        v.max_entropy_after_gaps = pd.Series.max(v.entropylist_S_gap_considered)
        self.logger.info('Maximum Shannon entropy alignment score after gap % considered: ' + str(round(v.max_entropy_after_gaps,2)))

        v.entropy_record_numpy = v.entropylist_S_gap_considered.to_numpy()
        v.entropy_record_numpy.shape = (-1,len(v.entropylist_S_gap_considered))
        self.logger.info('Removing gapped positions from analysis set.')
        v.sequence_dataframe_gap_considered = v.sequence_dataframe.drop(v.gap_percent_bool_index_remove,axis=1)
        self.logger.info("Elapsed time: ~ " + str(int(time.time() - v.start_time)) + " seconds.")
        self.logger.info('Alignment positions analyzed after ' + str(v.gap_value_cutoff) + '% gap cutoff: ' + str(v.length_of_alignment-len(v.gap_percent_bool_index_remove)))

        self.logger.info('Preparing sequences for comparison.')

        # Comparison time warning

        comparison_time_full_table_seconds = v.depth_of_alignment * v.depth_of_alignment * (v.length_of_alignment-len(v.gap_percent_bool_index_remove)) * 3.14E-8
        if comparison_time_full_table_seconds > 1800 and v.number_in_small_test == v.depth_of_alignment:
            self.logger.warning('\n***WARNING: An input alignment of this size may take a considerable amount of time')
            self.logger.warning('   if all pairwise sequence comparisons are performed.')
            self.logger.warning('   A sampling-based approach may be considered.')
            self.logger.warning('   For a sampling-based approach, take advantage of the -n, -t, and -s flags.\n')

        # Clear out used items from memory

        del v.sequence_dataframe
        gc.collect()

        # Avoid empty source dataframe

        if v.depth_of_alignment//v.number_in_small_test == v.depth_of_alignment/v.number_in_small_test:
            v.times_to_sample_max_keep = (v.depth_of_alignment//v.number_in_small_test)
        else:
            v.times_to_sample_max_keep = (v.depth_of_alignment//v.number_in_small_test) + 1

        # Enter trial loop

        self.logger.info('Beginning sequence trials.')

        for trial in range(v.min_trials_for_each_sequence):

            if v.min_trials_for_each_sequence > 1:
                self.logger.info("Trial: " + str(trial+1) + " of " + str(v.min_trials_for_each_sequence))

            v.sequence_dataframe_gap_considered = v.sequence_dataframe_gap_considered.sample(frac=1,random_state = v.seed) # shuffle master sequence dataframe, use user-defined or standard random seed
            v.sequence_dataframe_gap_considered_max_keep = v.sequence_dataframe_gap_considered # copy shuffled version for work below

            for j in range(v.times_to_sample_max_keep):
                if v.number_in_small_test != v.depth_of_alignment and (j+1)//50 == (j+1)/50:
                    self.logger.info('\rSample: '+str((j+1)) + ' of ' +str(v.times_to_sample_max_keep) + ' | Trial: ' + str(trial+1))
                max_times_tested = v.record_sequence_trial_results['Total_trials'].max()

                if max_times_tested > trial:

                    v.sequence_dataframe_gap_considered_max_keep = v.sequence_dataframe_gap_considered.loc[v.record_sequence_trial_results['Total_trials'] != max_times_tested]

                length_sequence_dataframe_gap_considered_max_keep = len(v.sequence_dataframe_gap_considered_max_keep)

                if length_sequence_dataframe_gap_considered_max_keep >= v.number_in_small_test:
                    number_to_choose = v.number_in_small_test

                elif length_sequence_dataframe_gap_considered_max_keep < v.number_in_small_test:
                    number_to_choose = length_sequence_dataframe_gap_considered_max_keep

                table_sample = v.sequence_dataframe_gap_considered_max_keep.iloc[0:number_to_choose,:]
                v.table_sample_numpy = table_sample.to_numpy() # convert pandas dataframe to numpy array
                v.table_sample_numpy = v.table_sample_numpy.astype(np.int8) # change datatype in an attempt to reduce memory immylogs.info
                v.table_sample_numpy_rows, table_sample_numpy_columns = v.table_sample_numpy.shape

            # Initiate numpy array for entropy calculation values

                v.entropy_array = np.empty((number_to_choose,number_to_choose),dtype=float)
                v.entropy_array[:] = np.nan

            # Calculations of match or not, and sum entropy values

                self.engine()

                v.entropy_DF = pd.DataFrame(v.entropy_array,index=np.arange(number_to_choose), columns=np.arange(number_to_choose))
                maximum_entropy_sum = v.entropy_DF.max(skipna=True)
                v.entropy_DF -= maximum_entropy_sum
                v.entropy_DF *= -1.0
                v.entropy_DF.columns = table_sample.index
                v.entropy_DF.index = table_sample.index
                entropy_DF_analysis_empty = np.empty((number_to_choose,1))
                entropy_DF_analysis_empty[:] = np.nan
                v.entropy_DF_analysis = pd.DataFrame(data = entropy_DF_analysis_empty, index=v.entropy_DF.index, columns=['Median'])

                for z in v.entropy_DF:
                    v.entropy_DF_analysis.loc[z,'Median'] = v.entropy_DF.loc[z,:].median(skipna=True)

                v.record_sequence_trial_results.loc[v.entropy_DF_analysis.index,'Total_trials'] += 1

        # Calculate interquartile range and outlier cutoff

                v.entropy_DF_analysis_values_list = v.entropy_DF_analysis.values.tolist()
                q25, q75 = np.nanpercentile(v.entropy_DF_analysis_values_list, 25), np.nanpercentile(v.entropy_DF_analysis_values_list, 75)
                iqr = q75 - q25

                CIQR = iqr * v.multiplier_on_interquartile_range
                v.lower_cutoff, v.upper_cutoff = q25 - CIQR, q75 + CIQR

                if p.write_plot_files:
                    self.write_median_plot()

        # Identify the outlier sequences using the interquartile range cutoff

                entropy_DF_analysis_above_cutoff = v.entropy_DF_analysis > v.upper_cutoff
                entropy_median_too_high = entropy_DF_analysis_above_cutoff.loc[entropy_DF_analysis_above_cutoff['Median'] == True]
                v.record_sequence_trial_results.loc[entropy_median_too_high.index,'Outlier_instances'] += 1

            self.logger.info("Elapsed time: ~ " + str(int(time.time() - v.start_time)) + " seconds.")
            self.logger.info("Estimated total time for analysis: ~ " + str(int(((time.time() - v.start_time))/(1+trial)*v.min_trials_for_each_sequence)) + " seconds.")

        if p.write_table_files:
            self.write_full_table()

        v.record_sequence_trial_results['Fraction_positive'] = v.record_sequence_trial_results['Outlier_instances'] / v.record_sequence_trial_results['Total_trials']

        if p.write_sequence_files:
            self.write_sequence_files()

        if p.write_plot_files:
            self.write_sampling_trials_plot()

        v.record_sequence_trial_results['Selected_for_retention'] = ''
        if v.stringency_flag == 1: # minimal stringency
            v.record_sequence_trial_results['Selected_for_retention'] = np.where(v.record_sequence_trial_results['Fraction_positive'] != 1.0, True, False)

        if v.stringency_flag == 2: # moderate stringency
            v.record_sequence_trial_results['Selected_for_retention'] = np.where(v.record_sequence_trial_results['Fraction_positive'] <= 0.5, True, False)

        if v.stringency_flag == 3: # maximal stringency
            v.record_sequence_trial_results['Selected_for_retention'] = np.where(v.record_sequence_trial_results['Fraction_positive'] == 0.0, True, False)

        if p.write_table_files:
            self.write_output_table()

        # Provide total runtime

        self.logger.info("Analysis complete.")
        self.logger.info("Total time for analysis: ~ " + str(int(((time.time() - v.start_time))/(1+trial)*v.min_trials_for_each_sequence)) + " seconds.")

        # Return final verdict for each sample

        return {x['Accession']: x['Selected_for_retention'] for x in v.record_sequence_trial_results.to_dict('records')}

    def initialize_logger(self):
        "Also adds file handler if required"
        v = self.vars
        p = self.params

        self.logger = logger.getChild(str(id(self)))
        if p.write_log_file:
            log_file = logging.FileHandler(v.output_entry + '_output.log')
            self.logger.addHandler(log_file)

        # Print header

        self.logger.info('\nSequenceBouncer: A method to remove outlier entries from a multiple sequence alignment\n')
        self.logger.info('Cory Dunn')
        self.logger.info('University of Helsinki')
        self.logger.info('cory.dunn@helsinki.fi')
        self.logger.info('Version: ' + version)
        self.logger.info('Please cite DOI: 10.1101/2020.11.24.395459')
        self.logger.info('___\n')

    def parse_params(self):
        "Parse parameters passed to class constructor"
        v = self.vars
        p = self.params

        if isinstance(p.input, InputSequence):
            v.input_name = p.input.name
            v.input_sequence = iter(p.input)
        elif isinstance(p.input, str) or isinstance(p.input, Path):
            v.input_name = str(Path(p.input))
            v.input_sequence = SeqIO.parse(p.input, p.input_format)
        elif isinstance(p.input, TextIOWrapper):
            v.input_name = p.input.name
            v.input_sequence = SeqIO.parse(p.input, p.input_format)
        else:
            raise InvalidInput()

        v.stringency_flag = p.stringency
        v.min_trials_for_each_sequence = p.trials
        v.multiplier_on_interquartile_range = p.IQR_coefficient
        v.number_in_small_test = p.subsample_size
        v.gap_value_cutoff = p.gap_percent_cut
        v.output_entry = p.output_file
        v.seed = p.random_seed or random.randint(0,1000)

        if v.output_entry is None:
            sep = '.'
            input_strip = v.input_name.split(sep, 1)[0]
            v.output_entry = input_strip
        v.output_sequence = input_strip + '_output_clean.fasta'
        v.output_rejected = input_strip + '_output_rejected.fasta'
        v.output_tabular = input_strip + '_output_analysis.csv'
        v.output_full_table = input_strip + '_full_comparison_table.csv'

        if p.write_none:
            p.write_log_file = False
            p.write_sequence_files = False
            p.write_table_files = False
            p.write_plot_files = False

    def parse_input(self):
        "Parse all input sequences"
        v = self.vars
        p = self.params

        self.logger.info('Generating sequence dataframe.')

        v.alignment_record_name_list = []
        v.record_seq_convert_to_string = []
        sequence_records = []

        for record in v.input_sequence:

            v.alignment_record_name_list.append(record.name)

            # Prepare dataframe to generate FASTA files

            if p.write_sequence_files:
                v.record_seq_convert_to_string.append(str(record.seq))

            # Load sequences from alignment into list and control case

            record_x_toward_seq_dataframe = list(record.seq)
            record_x_toward_seq_dataframe_lower = [x.lower() for x in record_x_toward_seq_dataframe]
            record_x_toward_seq_dataframe_ASCII = [ord(x) for x in record_x_toward_seq_dataframe_lower]
            sequence_records.append(record_x_toward_seq_dataframe_ASCII)

        # Generate dataframe of alignment from list

        v.sequence_dataframe = pd.DataFrame(sequence_records, dtype='int8')

        # Update input information

        v.depth_of_alignment = len(v.alignment_record_name_list)
        v.length_of_alignment = len(record.seq)

        if v.number_in_small_test == 0:
            v.number_in_small_test = v.depth_of_alignment

        if v.number_in_small_test == v.depth_of_alignment:
            v.min_trials_for_each_sequence = 1

    def engine(self):
        "Define the calculation engine"
        v = self.vars

        for counter_x in range(v.table_sample_numpy_rows):
            counter_x_numpy_row = v.table_sample_numpy[counter_x:(counter_x+1),:]
            if v.depth_of_alignment < 1000 and ((counter_x+1)/25) == ((counter_x+1)//25):
                self.logger.info('\rSequences analyzed: '+str(counter_x+1))
            elif v.depth_of_alignment < 10000 and ((counter_x+1)/250) == ((counter_x+1)//250):
                self.logger.info('\rSequences analyzed: '+str(counter_x+1))
            elif v.depth_of_alignment < 100000 and ((counter_x+1)/2500) == ((counter_x+1)//2500):
                self.logger.info('\rSequences analyzed: '+str(counter_x+1))
            for counter_y in range((counter_x+1)):
                counter_y_numpy_row = v.table_sample_numpy[counter_y:(counter_y+1),:]
                comparison_bool_series_match = counter_x_numpy_row == counter_y_numpy_row
                comparison_bool_series_NOT_match = counter_x_numpy_row != counter_y_numpy_row
                entropy_record_match = v.entropy_record_numpy[(comparison_bool_series_match)]
                entropy_record_NOT_match = v.entropy_record_numpy[(comparison_bool_series_NOT_match)]
                match_entropy_total = entropy_record_match.sum(axis=0)
                NOT_match_entropy_minus_max_entropy = entropy_record_NOT_match - v.max_entropy_after_gaps
                NOT_match_entropy_total =  NOT_match_entropy_minus_max_entropy.sum(axis=0)
                total_entropy_recorded = match_entropy_total + NOT_match_entropy_total
                v.entropy_array[counter_x, counter_y] = total_entropy_recorded
                v.entropy_array[counter_y, counter_x] = total_entropy_recorded
        return v.entropy_array

    def write_sequence_files(self):
        "Save clean and rejected FASTA files"
        v = self.vars

        acc_records_S = pd.Series(v.alignment_record_name_list)
        sequence_records_S = pd.Series(v.record_seq_convert_to_string)
        frame = { 'Accession': acc_records_S, 'Sequence': sequence_records_S }
        FASTA_output_unclean_DF = pd.DataFrame(frame)

        # Generating clean dataframes

        if v.stringency_flag == 1: # minimal stringency
            FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[v.record_sequence_trial_results['Fraction_positive'] != 1]
        if v.stringency_flag == 2: # moderate stringency
            FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[v.record_sequence_trial_results['Fraction_positive'] <= 0.5]
        if v.stringency_flag == 3: # maximal stringency
            FASTA_output_clean_DF = FASTA_output_unclean_DF.loc[v.record_sequence_trial_results['Fraction_positive'] == 0]

        # Generating rejection dataframes

        if v.stringency_flag == 1: # minimal stringency
            FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[v.record_sequence_trial_results['Fraction_positive'] == 1]
        if v.stringency_flag == 2: # moderate stringency
            FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[v.record_sequence_trial_results['Fraction_positive'] > 0.5]
        if v.stringency_flag == 3: # maximal stringency
            FASTA_output_reject_DF = FASTA_output_unclean_DF.loc[v.record_sequence_trial_results['Fraction_positive'] != 0]

        # Save clean FASTA file

        self.logger.info('Writing cleaned alignment as FASTA.')
        self.logger.info(FASTA_output_clean_DF)
        FASTA_output_final_acc_list = FASTA_output_clean_DF.loc[:,'Accession'].values.tolist()
        FASTA_output_final_seq_list = FASTA_output_clean_DF.loc[:,'Sequence'].values.tolist()

        with open(v.output_sequence, "w") as ofile:
            for seqi in range(len(FASTA_output_final_acc_list)):
                ofile.write(">" + FASTA_output_final_acc_list[seqi] + "\n" + FASTA_output_final_seq_list[seqi] + "\n")

        # Save rejected FASTA file

        self.logger.info('Writing rejected sequences to FASTA.')
        self.logger.info(FASTA_output_reject_DF)
        FASTA_output_rejected_acc_list = FASTA_output_reject_DF.loc[:,'Accession'].values.tolist()
        FASTA_output_rejected_seq_list = FASTA_output_reject_DF.loc[:,'Sequence'].values.tolist()

        with open(v.output_rejected, "w")as ofile:
            for seqi in range(len(FASTA_output_rejected_acc_list)):
                ofile.write(">" + FASTA_output_rejected_acc_list[seqi] + "\n" + FASTA_output_rejected_seq_list[seqi] + "\n")

    def write_gap_plot(self):
        "Plot gap fractions for alignment positions"
        v = self.vars

        v.gap_record.sort()
        plotgaplistrange = np.arange(len(v.gap_record))
        plt.plot(plotgaplistrange, v.gap_record, 'o', ms=1,c='slategrey')
        plot_cutoff_label = 'Selected gap fraction cut-off: ' + str(v.gap_value_cutoff/100) + ' (' + str(v.gap_value_cutoff) + ' %)'
        plt.axhline(y=v.gap_value_cutoff/100, color='teal', linestyle='--', label=plot_cutoff_label)
        rcParams['font.family'] = 'sans-serif'
        rcParams['font.sans-serif'] = ['Arial']
        plt.xlabel('Input alignment column (Ordered by gap fraction)', fontsize=12)
        plt.ylabel('Gap fraction', fontsize=12)
        plt.legend()
        plt.savefig(v.output_entry + '_gap_plot.pdf', format="pdf", bbox_inches="tight")
        self.logger.info('Printed gap distribution of input alignment to file: ' + v.output_entry + '_gap_plot.pdf')
        plt.close()

    def write_sampling_trials_plot(self):
        "Plot trial results from sampling-based approach"
        v = self.vars

        if v.number_in_small_test != v.depth_of_alignment:
            record_sequence_trial_results_list = v.record_sequence_trial_results['Fraction_positive'].tolist()
            record_sequence_trial_results_list.sort()
            plottrialrange = np.arange(len(record_sequence_trial_results_list))
            plt.plot(plottrialrange, record_sequence_trial_results_list, 'o', ms=1,c='deepskyblue')
            rcParams['font.family'] = 'sans-serif'
            rcParams['font.sans-serif'] = ['Arial']
            plt.xlabel('Sequence', fontsize=12)
            plt.ylabel('Fraction of times sequence called aberrant (In order of increasing positive calls)', fontsize=12)
            plt.savefig(v.output_entry + '_sampling_trials_plot.pdf', format="pdf", bbox_inches="tight")
            self.logger.info('Printed plot of sampling trial results to file: ' + v.output_entry + '_sampling_trials_plot.pdf')
            plt.close()

    def write_median_plot(self):
        "Plot comparison values, along with selected cut-off IQR cut-off value for full analysis"
        v = self.vars

        if v.number_in_small_test == v.depth_of_alignment:
            v.entropy_DF_analysis_values_list.sort()
            plotlistrange = np.arange(len(v.entropy_DF_analysis_values_list))
            plt.plot(plotlistrange, v.entropy_DF_analysis_values_list, 'o', ms=1,c='darkgreen')
            plot_cutoff_label = 'Selected IQR cut-off:  ' + str(v.multiplier_on_interquartile_range)
            plt.axhline(y=v.upper_cutoff, color='red', linestyle='--', label=plot_cutoff_label)
            rcParams['font.family'] = 'sans-serif'
            rcParams['font.sans-serif'] = ['Arial']
            plt.xlabel('Sequence', fontsize=12)
            plt.ylabel('Median across pairwise comparisons (Ordered by median value)', fontsize=12)
            plt.legend()
            plt.savefig(v.output_entry + '_median_plot.pdf', format="pdf", bbox_inches="tight")
            self.logger.info('Printed median values of sequence comparisons from full analysis to file ' + v.output_entry + '_median_plot.pdf')
            plt.close()

    def write_full_table(self):
        "Print full distance matrix for analysis and generate a plot only if a single test of all sequences versus all sequences was performed"
        v = self.vars

        if v.number_in_small_test == v.depth_of_alignment:

            self.logger.info('Cut-off value for median taken across comparisons (full-alignment pairwise analysis): ' + str(round(v.upper_cutoff,1)))
            v.entropy_DF.sort_index(axis=0,inplace=True,ascending=True)
            v.entropy_DF.sort_index(axis=1,inplace=True,ascending=True)
            v.entropy_DF.index = v.alignment_record_name_list
            v.entropy_DF.columns = v.alignment_record_name_list
            v.entropy_DF['Median_across_pairwise_comparisons'] = v.entropy_DF.median(axis=1) # add column calculating median across pairwise comparisons
            first_column = v.entropy_DF.pop('Median_across_pairwise_comparisons')
            v.entropy_DF.insert(0,'Median_across_pairwise_comparisons',first_column)
            v.entropy_DF.to_csv(v.output_full_table)

    def write_output_table(self):
        "Save analysis file as CSV"
        v = self.vars

        if v.min_trials_for_each_sequence == 1:
            del v.record_sequence_trial_results['Total_trials']
            del v.record_sequence_trial_results['Outlier_instances']
            del v.record_sequence_trial_results['Fraction_positive']

        v.record_sequence_trial_results.to_csv(v.output_tabular,index=False)
