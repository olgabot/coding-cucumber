#!/usr/bin/env python

# Parse command line arguments
import argparse

# R like data structures
import pandas as pd

# Scientific matrices
import numpy as np

# Plotting
import matplotlib.pyplot as plt

# for setting parameters
import matplotlib as mpl

# For highlighting ranges of plot greater than 0.5
from itertools import *
from operator import *

'''
Author: Olga Botvinnik
Date created: 11/15/12 09:12:52

The purpose of this program is to take in pondrfit data and make a pretty
plot using matplotlib

Example run:
python pondrfit_plots.py.py ...
'''

#######################################################################
# Class: CommandLine
#######################################################################
class CommandLine(object):
	def __init__(self, inOpts=None):
		self.parser = argparse.ArgumentParser(description = 
			'''
			Import a pondrfit data file and a title
			''',
			add_help=True, prefix_chars='-', 
			usage='%(prog)s --file FILE [options]')
		self.parser.add_argument('--file', '-f', action='store',
			help='the pondrfit datafile to import and query', required=True)
		self.parser.add_argument('--title', '-t', action='store',
			type=str, default=10,
			help='title of the plot.')
		self.parser.add_argument('--color', '-c', action='store',
			help='color to use for the plot')
		self.parser.add_argument('--show-range-numbers', '-s', action='store_true',
			help='For the predicted disordered regions, show the range numbers')
		# self.parser.add_argument('--sampleA', '-sA', action='store',
		# 	help='For finding differentially expressed genes between sample A\
		# 	and sample B at time point X')
		# self.parser.add_argument('--sampleB', '-sB', action='store',
		# 	help='For finding differentially expressed genes between sample A\
		# 	and sample B at time point X')
		# self.parser.add_argument('--timeX', '-tX', action='store',
		# 	help='For finding differentially expressed genes between sample A\
		# 	and sample B at time point X')
		if inOpts is None:
			self.args = vars(self.parser.parse_args())
		else:
			self.args = vars(self.parser.parse_args(inOpts))

	def do_usage_and_die (self, str):
		'''
		If a critical error is encountered, where it is suspected that the 
		program is not being called with consistent parameters or data, this
		method will write out an error string (str), then terminate execution 
		of the program.
		'''
		import sys
		print >>sys.stderr, str
		self.parser.print_usage()
		return 2


#######################################################################
# Class: Usage
#######################################################################
class Usage(Exception):
	'''
	Used to signal a Usage error, evoking a usage statement and eventual
	exit when raised
	'''
	def __init__(self, msg):
		self.msg = msg

#######################################################################
# Function: plot_lines
#######################################################################
def plot_lines(df, title, color):
	'''
	Plots the PONDR-FIT data using a Tufte-style plot
	'''
	n = len(df.X0)

	ax = plt.axes((0.15, 0.15, 0.8, 0.8))
	plt.plot(df.X0, df.X2, color=color)
	plt.plot(df.X0, 0.5*np.ones(n), color='black')

	ax.set_xlim((df.X0[0], n))
	ax.set_ylim(( 0,1))
	#ax.set_visible(
	ax.spines['top'].set_visible(False)
	ax.spines['right'].set_visible(False)
	ax.spines['bottom'].set_position(('outward', 20))
	ax.spines['left'].set_position(('outward', 30))
	step = int(np.round(n/4.0, -1))
	xticks = range(0,n+1,step)
	xticks[0] = 1
	if xticks[-1] != n:
		xticks.append(n)
	ax.xaxis.set_ticks(xticks)
	ax.yaxis.set_ticks_position('left')
	ax.xaxis.set_ticks_position('bottom')
	plt.title('%s disorder prediction via PONDR-FIT' % title)
	plt.xlabel('position along protein')
	plt.ylabel('Probability of disorder according to PONDR-FIT')
	return ax

#######################################################################
# Function: highlight_disordered
#######################################################################
def highlight_disordered(df, ax, show_range_numbers):
	'''
	Highlight the disordered regions of the protein, aka the ones that
	have a disorder prediction score greater than 0.5
	'''
	# highlight values greater than 0.5
	df_disordered = df.ix[df.X2>0.5,:]

	ranges = []
	for k, g in groupby(enumerate(df_disordered.index), lambda (i,x):i-x):
		inds = map(itemgetter(1), g)
	#    print inds
		ranges.append([inds[0], inds[-1]])
	prev_y = 0.125
	for i in range(len(ranges)):
		ind1 = ranges[i][0]
		ind2 = ranges[i][1]
		#print 'ind1:', ind1, '  ind2:', ind2
		start = df_disordered.X0[ind1]
		stop = df_disordered.X0[ind2]
		# print 'start:', start, '\tstop:', stop
		plt.axvspan(start, stop, color='grey', alpha=0.25)
		
		if show_range_numbers:
			diff = start - ranges[i-1][0]
			max_diff = df.shape[0]/4
			if i == 0:
				ha = 'left'
				x = start
			elif i == len(ranges)-1:
				ha = 'right'
				x = stop
			elif start == stop:
				ha = 'center'
				x = start
			elif diff < max_diff and prev_y == 0.05:
				ha = 'right'
				x = stop
			elif start == 269:
				ha = 'center'
				x = (start+stop)/2.0
			else:
				ha = 'left'
				x = start
				
			# if the previous range was close to this range

			# print 'diff:', diff, '\tmax_diff:', max_diff
			if start == stop:
				y = 0.2
			elif i > 0 and diff < max_diff and prev_y == 0.05:
				y = 0.125
			elif diff < max_diff and prev_y == 0.125:
				y = 0.2
			else:
				y = 0.05
				
			prev_y = y
			ax.text(x=x, y=y, s='%d-%d' % (start, stop),
				color='black', size=10, ha=ha, va='top', 
				backgroundcolor='white')


#######################################################################
# Function: main
#######################################################################
def main():
	'''
	This function is invoked when 
	'''
	cl = CommandLine()
	try:
		# Debugging statement: print out the command line arguments
		# to make sure they were parsed correctly
		print cl.args

		# mpl.rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})

		# mpl.rcParams['font.sans-serif'] = 'Helvetica'

		df = pd.read_table(cl.args['file'], header=None, sep=' +')

		ax = plot_lines(df, cl.args['title'], cl.args['color'])

		highlight_disordered(df, ax, cl.args['show_range_numbers'])

		plt.savefig('pondrfit.pdf')

	# If not all the correct arguments are given, break the program and
	# show the usage information
	except Usage, err:
		cl.do_usage_and_die(err.msg)


if __name__ == '__main__':
	main()