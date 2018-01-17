#!/usr/bin/python
"""
Combines FASC files from PyRosetta job distributor and converts them into a 
tabular form
"""
import argparse
from glob import glob
from os.path import basename, join

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("directory", type=str, 
		help="Read fasc files from what directory?")
	args = parser.parse_args()
	return args


def main():
	# Getting user inputs
	args = parse_args()

	# Getting fasc files
	folder = args.directory
	folder_search = join(folder,"*.fasc")
	fasc_files = glob(folder_search)
	fasc_files.sort()

	# Making combined report
	report_name = folder.rstrip('/') + '_combined_reports.fasc'
	with open(report_name, 'w') as r:
		for f in fasc_files:
			# Reading in individual FASC file
			with open(f, 'r') as read:
				f_lines = read.readlines()
				f_lines.pop(0) # First line is not useful
				lines_data = [i.split()[1::2] for i in f_lines]
				lines_data.sort()

			# Making template and header
			if f == fasc_files[0]:
				headline = f_lines[0].split()[::2]
				template = '{:50s}' + '{:25s}' * (len(headline) - 3) 
				template += '{:110s}' * 2 + '\n'
				r.write(template.format(*headline))

			# Adding in lines
			for line in lines_data:
				line[0] = basename(line[0])
				r.write(template.format(*line))


if __name__ == '__main__':
	main()