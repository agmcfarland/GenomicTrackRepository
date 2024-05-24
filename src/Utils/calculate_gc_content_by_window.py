from Bio import SeqIO
import os
from os.path import join as pjoin
import pandas as pd
import itertools
import sys

def chunk_sequence(sequence, window_size):
	"""
	Split the string into equal-sized chunks of size 'size'.
	"""
	return [sequence[i:i + window_size] for i in range(0, len(sequence), window_size)]

def ignore_record(record_id):
	"""
	Chromosome names (records) that are non-standard are set as ignored via True/False.
	"""
	for i in ['chrUn', '_fix', '_alt', '_random']:
		if record_id.find(i) > -1:
			return True
	return False

def main(genomic_fasta_file_path: str, output_file_path: str, window_size: int):

	df_gc_content = []

	for record in SeqIO.parse(genomic_fasta_file_path, 'fasta'):

		if ignore_record(record.id):
			continue

		print(record.id)

		record_bp_length = len(str(record.seq))

		start_window = 0

		for genomic_window_ in chunk_sequence(sequence = str(record.seq), window_size = window_size):
			a_count = genomic_window_.count('A')
			t_count = genomic_window_.count('T')
			c_count = genomic_window_.count('C')
			g_count = genomic_window_.count('G')

			try:
				gc_percentage = (c_count + g_count) / (a_count + t_count + c_count + g_count) * 100 # Don't divide by length but by counts of A,C,T,G because of Ns.
			except ZeroDivisionError:
				gc_percentage = 0  # Case denominator is all Ns

			end_window = start_window + window_size - 1

			if end_window > record_bp_length:
				end_window = record_bp_length

			df_gc_content.append([record.id, start_window, end_window, gc_percentage])

			start_window = end_window + 1


	df_gc_content_df = pd.DataFrame(df_gc_content, columns=['chrom', 'start_window', 'end_window', 'gc_percentage'])

	df_gc_content_df.to_csv(output_file_path, index = None)

if __name__ == '__main__':
	genomic_fasta_file_path = sys.argv[1] #pjoin('/data/GenomicTrackRepository/data/raw/hg38.fa')
	output_file_path = sys.argv[2] #pjoin('/data/GenomicTrackRepository/data/raw/output_file.csv')
	window_size = int(sys.argv[3]) # 100000000

	# genomic_fasta_file_path = pjoin('/data/GenomicTrackRepository/data/raw/hg38.fa')
	# output_file_path = pjoin('/data/GenomicTrackRepository/data/raw/output_file.csv')
	# window_size = 100000000

	print('Reading from:', genomic_fasta_file_path)
	print('Writing to:', output_file_path)
	print('Window size:', window_size)

	main(genomic_fasta_file_path, output_file_path, window_size)	





