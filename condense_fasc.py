#!/usr/bin/python
"""
Combines FASC files from PyRosetta job distributor and converts them into a 
tabular form. 
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


def cleanup_mutations_section(report_lines, start_point):
	""" 
	Some sets are missing mutable residues present in others. This function 
	maintiains column alignment across all sets by inserting what will appear
	as blank cells where a mutable residue is missing.
	"""
	max_len = max([len(line) for line in report_lines])
	res_columns = range(start_point, max_len, 4)

	for c in res_columns:
		first_des_res = min([int(line[c]) for line in report_lines])
		for line in report_lines:
			if int(line[c]) != first_des_res:
				for i in range(4):
					line.insert(c, '=""')


def main():
	# Getting user inputs
	args = parse_args()

	# Getting fasc files
	folder = args.directory
	folder_search = join(folder, "*.fasc")
	fasc_files = glob(folder_search)
	fasc_files.sort()
	# Prevent self-reference if re-running
	for n, i in enumerate(fasc_files):
		if 'combined_reports' in i:
			fasc_files.pop(n)

	# Collecting fasc lines
	headline = []
	report_lines = []
	mut_section_ind = 0
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
				if i == f_lines[0] and f == fasc_files[0]:
					headline = scores_section[::2]
					mut_section_ind = len(scores_section) / 2 + 1
				lines_data.append(line_data)
			lines_data.sort(key=lambda x: x[1]) # Sorting by total score
			report_lines += lines_data

	# Making combined report
	report_name = folder.rstrip('/') + '_combined_reports.fasc'
	report_name = join(folder, report_name)

	with open(report_name, 'w') as r:
		# Making template and header
		head_length = len(headline)
		template = '{:50s}' + '{:25s}' * (head_length - 1)
		r.write(template.format(*headline) + '\n')

		# Adding in lines
		cleanup_mutations_section(report_lines, mut_section_ind)
		for line in report_lines:
			line[0] = basename(line[0])
			line_out = template.format(*line[:head_length])
			line_out += '   '.join(line[head_length + 1:])
			r.write(line_out + '\n')

		print report_name

if __name__ == '__main__':
	main()