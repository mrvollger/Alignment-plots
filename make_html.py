#!/usr/bin/env python
import argparse
import os 
import sys
import pandas as pd
import re
import altair as alt
from vega_datasets import data
import json
alt.themes.enable("default")

# global var for inputs
args=None 


if __name__ == "__main__":
	parser = argparse.ArgumentParser(description="", formatter_class=argparse.ArgumentDefaultsHelpFormatter)
	parser.add_argument("infile", help="positional input")
	parser.add_argument("outfile", help="output html")
	parser.add_argument("-s", "--string", help="string option")
	parser.add_argument("-t", "--height", help="", type=int, default=800)
	parser.add_argument("-w", "--width", help="", type=int, default=1000)
	parser.add_argument("-l", "--list", nargs="*", help="list with zero or more entries")
	parser.add_argument("-l2", "--list2", nargs="+", help="list one or more entries")
	parser.add_argument('-d', help="store args.d as true if -d",  action="store_true", default=False)
	args = parser.parse_args()
	

	H=args.height
	W=args.width
	tbl = pd.read_csv(args.infile, sep="\t")
	tbl.cutid = tbl.cutid.astype('category')
	aln=tbl.loc[:,[ "cutid","identity","aln_len","x1","x2","y1","y2","q_name","t_name"] ]
	aln["label"] = ( "(" + aln.x1.astype(str) + "," + aln.y1.astype(str) + ") (" 
					+ aln.x2.astype(str) + ", " + aln.y2.astype(str) + "); " 
					+ aln.identity.astype(str) + "%" )


	contigs = sorted(list(set(list(aln.q_name) + list(aln.t_name))))
	

	#
	# altair code 
	#
	chart = alt.Chart(aln)

	brush = alt.selection(type='interval', encodings=["x"])
	color = alt.Color('cutid:Q', sort="descending", scale=alt.Scale(scheme='spectral'), legend=None)
	mycolor = alt.condition(brush, color, alt.value('lightgray'))
	single_nearest = alt.selection_single(on='mouseover', nearest=True, empty='none')
	scales = alt.selection_interval(bind='scales')

	t_name = alt.binding_select(options=contigs)
	q_name = alt.binding_select(options=contigs)
	select_t = alt.selection_single(fields=['t_name'], bind=t_name, name="Target", init={"t_name":contigs[0]})
	select_q = alt.selection_single(fields=['q_name'], bind=q_name, name="Query", init={"q_name":contigs[0]})




	segs = chart.mark_line().encode(
		x=alt.X('x1:Q', title='Target position (bp)'),
		x2="x2:Q",
		y=alt.Y("y1:Q", title='Query position (bp)'),
		y2="y2:Q",
		color=mycolor,
		size=alt.condition(single_nearest, alt.value(3), alt.value(1.5))
	)
	labels = segs.mark_text(
		align='left',
		baseline='middle',
		dx=3,
	).encode(
		text=alt.condition(single_nearest, "label:N", alt.value(' ')),
		size=alt.value(10),
		color=alt.value("black")
	).add_selection(
		single_nearest
	)
	lines = (segs + labels).properties(
		width=W,
		height=H
	).add_selection(
		scales
	)




	points = chart.mark_circle().encode(
		x=alt.X('identity:Q', title='% identity',  scale=alt.Scale(domain= [min(aln.identity),100] ) ),
		y=alt.Y('aln_len:Q', title='Alignment length (bp)',  scale=alt.Scale(type='log')),
		color=mycolor,
		size=alt.condition(single_nearest, alt.value(128), alt.value(16))
	).properties(
		width=W/2,
		height=H/3,
	).add_selection(
		brush,
	).transform_filter(
		scales
	)


	bars = chart.mark_bar().encode(
		x=alt.X("sum(aln_len):Q", title="Total aligned bases"),
		y=alt.Y("cutid:N", title="% identity"),
		color=color,
	)
	text = bars.mark_text(
		align='left',
		baseline='middle',
		dx=3
	).encode(text="sum(aln_len):Q")

	hist = (bars + text).properties(
		width=W/2,
		height=H/3
	).transform_filter(
		scales
	).transform_filter(
		brush
	)

	plot = ( lines &  (points | hist) ).add_selection(select_q, select_t).transform_filter(select_q).transform_filter(select_t)
	plot.save(args.outfile)








