#!/usr/bin/env python
import argparse
import os 
import sys
import numpy as np
import re
import pandas as pd
import matplotlib
matplotlib.use('agg')
from matplotlib import cm
import matplotlib.pyplot as plt
import seaborn as sns
import scipy.sparse as sparse

def read_tbl(FILE, NCOLORS):
	df = pd.read_csv(FILE, sep = "\t")
	MINID = int(df.perID_by_all.quantile(.005))
	df = df.loc[df.perID_by_all >= MINID]
	df["q_pos"] = df["query_name"].str.extract(".*:(\d+)-\d+").astype(int)
	df["r_pos"] = df["reference_name"].str.extract(".*:(\d+)-\d+").astype(int)
	df["first_pos"] = df[["q_pos", "r_pos"]].min(axis=1)
	df["second_pos"] = df[["q_pos", "r_pos"]].max(axis=1)
	df["cut"] = pd.qcut( df["perID_by_all"], NCOLORS, duplicates="drop")
	df["cutid"] = pd.qcut( df["perID_by_all"], NCOLORS, duplicates="drop", labels=False)
	NCOLORS = len(df.cutid.unique())
	cmap = cm.get_cmap('plasma', NCOLORS+1)
	df["color"] = df["cutid"].map(cmap)
	# modal alignment length 
	return(df, cmap)

def plot(df, OUT, cmap):
	fig, (axh, ax) = plt.subplots(2, 1, gridspec_kw={'height_ratios': [1, 4]}, figsize=(8,11), dpi=600)
	
	# get a marker size that is relative to the x axis size
	ALNSIZE = (df.reference_end-df.reference_start).mode()[0]
	dim = max(df.first_pos.max(), df.second_pos.max()) + 1
	s = ALNSIZE*((ax.get_window_extent().width  / (dim) * 72./fig.dpi) ** 2)
		

	# make the histogram
	step=.25
	axh.hist( [g.perID_by_all for n, g in df.groupby("color")], bins=np.arange(0,100+step, step), 
			 histtype='bar', stacked=True, color=cmap(range(df.cutid.max() + 1)), label=df.color, log=True)
	axh.set_xlim(int(df.perID_by_all.min()), 100)

	# make the dotplot 
	# draw diagonal
	ax.set_xlim(0,dim); ax.set_ylim(0,dim)
	ax.plot(ax.get_xlim(), ax.get_ylim(), color = cmap(df.cutid.max()), markersize=s)
	# draw off diagonal
	for color, group in df.groupby("color"):
		z, x, y = group.perID_by_all.to_numpy(), group.first_pos.to_numpy(), group.second_pos.to_numpy()   
		tmat  = sparse.coo_matrix((z, (y,x)), shape=(dim, dim))
		tmat2 = sparse.coo_matrix((z, (x,y)), shape=(dim, dim))
		ax.spy( tmat,  color=color, markersize= s , origin="lower")
		ax.spy( tmat2, color=color, markersize= s , origin="lower")

	ax.get_xaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))
	ax.get_yaxis().set_major_formatter(matplotlib.ticker.FuncFormatter(lambda x, p: format(int(x), ',')))

	plt.tight_layout()
	plt.savefig(OUT)




if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("infile", help="positional input")
	parser.add_argument("outfile", help="positional input")
	parser.add_argument("-n", "--ncolors", help="numeric option", type=int, default=10)
	parser.add_argument("-l", "--list", nargs="*", help="list with zero or more entries")
	parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()
	

	df, cmap = read_tbl(args.infile, args.ncolors)
	plot(df, args.outfile, cmap)



