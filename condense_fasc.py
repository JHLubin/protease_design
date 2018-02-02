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
	# Prevent self-reference if re-running
	for n, i in enumerate(fasc_files):
		if 'combined_reports' in i:
			fasc_files.pop(n)

	# Making combined report
	report_name = folder.rstrip('/') + '_combined_reports.fasc'
	report_name = join(folder, report_name)
	with open(report_name, 'w') as r:
		for f in fasc_files:
			# Reading in individual FASC file
			with open(f, 'r') as read:
				f_lines = read.readlines()
				f_lines.pop(0) # First line is not useful
				lines_data = []
				for i in f_lines:
					start_mutations = i.find('mutations')
					scores_section = i[:start_mutations].split()
					prot_mutations_section = i[start_mutations:].split()
					line_data = scores_section[1::2] + prot_mutations_section
					if i == f_lines[0]:
						headline = scores_section[::2]
					lines_data.append(line_data)
				lines_data.sort()

			if f == fasc_files[0]:
				# Making template and header
				head_length = len(headline)
				template = '{:50s}' + '{:25s}' * (head_length - 1)
				r.write(template.format(*headline) + '\n')

			# Adding in lines
			for line in lines_data:
				line[0] = basename(line[0])
				line_out = template.format(*line[:head_length])
				line_out += '   '.join(line[head_length + 1:])
				r.write(line_out + '\n')

	print report_name

if __name__ == '__main__':
	main()