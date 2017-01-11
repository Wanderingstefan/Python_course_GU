#!/usr/local/env python2.7

try:
	import sys, os, multiprocessing, subprocess, argparse, re
	import pandas as pd
	from argparse import RawTextHelpFormatter

except Exception, e:
	print '%s\nMissing dependencies!'
	sys.exit()
	

"""Converts fastq to fasta files and then maps nucleotide reads to CARD protein db. The best hits are chosen, reads per arg counted and written to csv. Takes directory containing the reads and number of processes to run as input."""


#Converts trimgalore output files (which are recognized by their ending) to fasta files for further analysis. 
def fastq_to_fasta(fq_file):
	print 'converting fastq to fasta...'
	file_end=fq_file.split('.')[-1]
	conversion='fastq_to_fasta -i %s -o %s' % (fq_file, fq_file.replace(file_end, '_reads.fasta'))
	subprocess.call(conversion, shell=True)

#Calls a subprocess using diamond for blasting fasta nucleotide reads against card protein database, saves mapping info to csv.
def map_reads(read_file):
	print 'mapping reads against CARD Resistance Gene Database...'
	CARD_diamond='/storage/stefan/PhD_project/context_db/bacterial_genomes_ncbi/test_genomes/blast_resdbs/card/card_protein.dmnd'
	diamond='/storage/stefan/PhD_project/context_db/tools/diamond/diamond blastx -q %s -p 10 -d %s -o %s -f 6 qseqid sseqid evalue score stitle slen' % (read_file, CARD_diamond, read_file.replace('_reads.fasta', '_reads_mapped.csv'))
	subprocess.call(diamond, shell=True, env=dict(ENV='~/.bashrc'))
	print 'reads mapped for %s' % read_file

#Uses the functions below to pick best read alignment for each read and count reads per ARG
def length_normalize_hits(mapped_reads):
	raw_mapped=pd.read_csv(mapped_reads, sep='\t', names = ['qseqid', 'sseqid', 'evalue', 'score', 'stitle', 'slen'])
	filtered_reads=pick_best_alignment(raw_mapped)
	counts=count_reads(filtered_reads, mapped_reads)
	
#Takes information from read mapping and selects the best hit for each read, returning a dataframe in which each read maps only once	
def pick_best_alignment(raw_mapped):
	filtered_reads=pd.DataFrame()
	grouped_reads=raw_mapped.groupby(['qseqid'])
	for read, group in grouped_reads:
		group_df=pd.DataFrame()
		for index, row in group.iterrows():
			if row['score']==max(group['score']):
				group_df=group_df.append(row)
		if len(group_df)>1:
			best_hit=group_df.sample(n=1)
		else:
			best_hit=group_df
		filtered_reads=filtered_reads.append(best_hit)
	print 'Best hits determined...'
	return filtered_reads

#Counts how many reads in total and how many reads per base map to each resistance gene
def count_reads(filtered_reads, mapped_reads):
	ARGs=filtered_reads.groupby(['sseqid'])
	count_frame=pd.DataFrame()
	for name, group in ARGs:
		data=[name, len(group), max(group['slen']), round(float(len(group)/max(group['slen'])), 2)]
		count_frame=count_frame.append(pd.Series(data, index=['ARG', 'total_read_count', 'ARGlen', 'per_base_read_count']), ignore_index=True)
	count_frame.to_csv(mapped_reads.replace('_reads_mapped.csv', '_arg_count.csv'), sep=',')
	print 'count data written to .csv ...'

def parse_arguments():
	man_description='%r\n\nConverts fastq to fasta files and then maps nucleotide reads to CARD protein db. The best hits are chosen, reads per arg counted and written to csv. Takes directory containing the reads and number of processes to run as input\n%r' % ('_'*80, '_'*80)
	parser=argparse.ArgumentParser(description=man_description.replace("'", ""), formatter_class=RawTextHelpFormatter)
	parser.add_argument('target_dir', help='directory containing the files to be analyzed.')
	parser.add_argument('-db', '--database', help='path to diamond sequence database.')
	parser.add_argument('-dmnd', '--path_to_diamond', help='path to diamond execultable.')
	parser.add_argument('-p', '--processes', help='number of processes excuting the script in parallel.', type=int)
	args=parser.parse_args()
	return args			

###MAIN###

args=parse_arguments()
fq_root=args.target_dir.rstrip('/')

fq_files=[]
for file in os.listdir(fq_root):
	galore_output=re.search(r'val_[1-2]\.(fq|fastq)$', file)
	if galore_output:
		fq_files.append(fq_root+'/'+file)

pool=multiprocessing.Pool(args.processes)
pool.map(fastq_to_fasta, fq_files)
pool.close()
pool.join()

read_files=[]
for file in os.listdir(fq_root):
	if file.endswith('_reads.fasta'):
		read_files.append(fq_root+'/'+file)

for read_file in read_files:
	map_reads(read_file)


mapped_readfiles=[]
for file in os.listdir(fq_root):
	if file.endswith('_reads_mapped.csv'):
		mapped_readfiles.append(fq_root+'/'+file)

pool=multiprocessing.Pool(args.processes)
pool.map(length_normalize_hits, mapped_readfiles)
pool.close()
pool.join()

print 'All samples processed!'

