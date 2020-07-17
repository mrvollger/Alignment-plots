#!/usr/bin/env python
import argparse
import os 
import sys
import pysam
import pandas as pd
import numpy as np
# global var for inputs
args=None 

if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("--gff", help="string option")
	parser.add_argument("--ref", help="string option")
	parser.add_argument("--gene", help="string option")
	parser.add_argument("--ogff", help="string option", )
	parser.add_argument("--ofasta", help="string option", )
	parser.add_argument("-n", "--number", help="numeric option", type=int, default=5)
	parser.add_argument("-l", "--list", nargs="*", help="list with zero or more entries")
	parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()
	

	fasta = pysam.FastaFile(args.ref)
	names = ['contig', 'source', 'feature', 'start', 'end', 'score', 'strand', 'frame', 'attributes']
	df = pd.read_csv(args.gff, sep='\t', comment='#', names = names)
	#df.drop_duplicates(inplace=True)
	ofasta = open(args.ofasta, 'w+')

	for contig, gf in df.groupby('contig'):
		mn = gf.start[gf.attributes.str.contains(args.gene)].min()
		mx = gf.end[gf.attributes.str.contains(args.gene)].max()
		if(np.isnan(mn) or np.isnan(mx)):
			continue
		gf = gf.copy()
		gf = gf[ (gf.start >= mn) & (gf.end <= mx)  ]

		seq = fasta.fetch(reference = contig, start = mn -1, end = mx)
		o = f'>{contig}:{mn}-{mx}\n{seq}\n'
		gf['start'] = gf.start - mn + 1
		gf['end'] = gf.end - mn + 1

		gf.to_csv(args.ogff, mode='a', sep='\t', index=False, header=None)
		ofasta.write(o)

		if(False):
			#print(gf.attributes.str.split('([^=;]+)=([^=;]+);',expand=True))
			df['gene_id'] = gf.attributes.str.extract('gene_id=([^=;]+);',expand=True)
			df['gene_name'] = gf.attributes.str.extract('gene_name=([^=;]+);',expand=True)

			for gene_id, gene in df.groupby('gene_id'):
				print(gene[gene.feature=='exon'])

			
