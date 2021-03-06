#!/usr/bin/python
""" 
Generate PyMOL scripts for viewing design_protease models
Intended to be run after condense_fasc.py
Joseph Lubin, 2017
"""
import argparse
from glob import glob
from os.path import basename, join

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("folder", help="Pick folder to check")
	args = parser.parse_args()
	return args


def cluster_designs(score_lines, ind):
	""" 
	From a list of score lines from design_protease.py, groups decoys by a 
	profile of mutations. Output is still a list of scores lines.
	"""
	mutation_profiles = {}
	for l in score_lines:
		m_cols = l[ind:]
		m_prof = ''

		# Looping through mutations columns to determine profile
		for i in range(2, len(m_cols), 4):
			if m_cols[i] not in ['NO', '=""']:
				locus = m_cols[i-2]
				mutation = m_cols[i]
				m_prof += locus + mutation + '_'

		# Adding line to appropriate group by profile
		m_prof = m_prof.strip('_')
		if m_prof != '': # Ignore unchanged models
			if m_prof in mutation_profiles:
				mutation_profiles[m_prof].append(l)
			else:
				mutation_profiles[m_prof] = [l]		

	return mutation_profiles


class design_group():
	def __init__(self, profile, score_lines):
		self.profile = profile
		self.scorelines = score_lines
		self.score_head = score_lines.pop(0)
		self.mutated_loci = [int(i[:-1]) for i in self.profile.split('_')]

		# Don't need to view too many PDBs
		self.included_pdbs = []
		self.excluded_pdbs = []
		self.prune_pdbs()
		self.excluded_pdbs.sort()

		self.pdbs = [i[0] for i in self.scorelines]
		relax_pdbs = [i.replace('design', 'relax') for i in self.pdbs]
		self.pdbs = relax_pdbs + self.pdbs

	def prune_pdbs(self):
		"""
		Checks list of score lines. If the list is longer than 5, it will cut
		back down based on scores. The best two decoys 
		each for total score, residue score, and ddG will be selected. In some
		cases, this will be redundant, so less than 6 may be represented.
		"""
		s = self.scorelines
		if len(self.scorelines) <= 5:
			for line in s:
				decoy = line[0]
				decoy_number = decoy.split('_')[-1].replace('.pdb', '')
				self.included_pdbs.append(int(decoy_number))			
			return

		else:
			# Identifying two best scores for each category
			total_column = self.score_head.index('total_score:')
			total_scores = [float(line[total_column]) for line in s]
			total_select = sorted(total_scores)[:2]

			res_column = self.score_head.index('residue_scores:')
			res_scores = [float(line[res_column]) for line in s]
			res_select = sorted(res_scores)[:2]

			ddg_column = self.score_head.index('ddG:')
			ddg_scores = [float(line[ddg_column]) for line in s]
			ddg_select = sorted(ddg_scores)[:2]

			# Keeping lines with best scores
			lines_to_keep = []
			for line in s:
				if any([float(line[total_column]) in total_select, 
					float(line[res_column]) in res_select, 
					float(line[ddg_column]) in ddg_select]):

					lines_to_keep.append(line)

			# Listing rejected lines
			for line in s:
				decoy = line[0]
				decoy_number = decoy.split('_')[-1].replace('.pdb', '')
				if line in lines_to_keep:
					self.included_pdbs.append(int(decoy_number))
				else:
					self.excluded_pdbs.append(int(decoy_number))

			self.scorelines = lines_to_keep
			self.included_pdbs.sort()
			self.excluded_pdbs.sort()
			return


def write_pymol(odir, seq, clusters):
	""" 
	Writes a PyMOL script to load given relaxed and designed PDB files and 
	format the view nicely for review of mutations. 
	"""
	comments = {}
	cmds = ['set seq_view, 1\n']
	for c in clusters:
		#cmds.append('hide everything')

		# Load in PDB files
		loadfiles = [join(odir, p) for p in c.pdbs]
		for p in loadfiles:
			cmds.append('load ' + p)

		# Making selections list for full PDB objects
		p_objs = [p.replace('.pdb', '') for p in c.pdbs]
		cmds.append('select current_pdbs, ' + ' + '.join(p_objs))
		first_pdb = p_objs[0]
		cmds.append('select first_pdb, ' + first_pdb)

		# Selecting main features
		cmds.append('select relaxed, *relaxed* and current_pdbs')
		cmds.append('select designed, *designed* and current_pdbs')
		cmds.append('select peptide, chain B and current_pdbs')
		cmds.append('select pep_recognition, chain B and res 2-6 and current_pdbs')
		cmds.append('select cat_res, res 72+96+154')

		# Selecting mutated and designable residues
		drs = [58, 70, 138, 147, 150, 151, 152, 153, 169, 170, 
				171, 172, 173, 174, 175, 176, 177, 183]
		for n, i in enumerate(drs):
			if i in c.mutated_loci:
				drs.pop(n)
		cmds.append('select des_res, current_pdbs and res ' + '+'.join([str(i) for i in c.mutated_loci]))
		cmds.append('select contact_res, current_pdbs and res ' + '+'.join([str(i) for i in drs]))

		# Coloring
		cmds.append('color green, relaxed')
		cmds.append('color cyan, designed')
		cmds.append('color magenta, cat_res')
		cmds.append('color yellow, peptide and designed')
		cmds.append('color yellow, peptide and first_pdb and n. CA+C+N')
		cmds.append('color purpleblue, des_res and designed')
		cmds.append('util.cnc')

		# Showing desired view
		cmds.append('hide')
		cmds.append('show cartoon, first_pdb')
		cmds.append('show lines, cat_res and first_pdb')
		cmds.append('show lines, pep_recognition des_res contact_res')
		cmds.append('hide lines, name C+N+O')
		cmds.append('hide everything, elem H')
		cmds.append('label (first_pdb and des_res and n. CA), "(%s, %s)" % (resi, resn)')

		# Saving scene
		cmds.append('deselect')
		cmds.append('scene ' + c.profile + ', store')

		# Adding report comments for which PDBs are in each set, and which were excluded
		comments[c.profile] = [c.included_pdbs, c.excluded_pdbs]

		# Hiding and disabling to be clean for next set
		for p in p_objs:
			cmds.append('hide everything, ' + p)
			cmds.append('disable ' + p)

		cmds.append('')

	# Outputting comments
	for p, [i, e] in comments.items():
		comment = '"' + p + ' included decoys: '
		comment += ', '.join([str(n) for n in i])
		if len(e) == 0:
			comment += '; no exclusions'
		else:
			comment += '; excluded decoys: '
			comment += ', '.join([str(n) for n in e])
		comment += '"'
		cmds.append('print ' + comment)

	# Adding line breaks
	cmds = [i + '\n' for i in cmds]

	# Writing pml script
	script_name = join(odir, 'view_' + seq + '.pml')
	with open(script_name, 'w') as w:
		w.writelines(cmds)
	print script_name


def main(folder):
	# Collecting decoys
	search_space = join(folder, '*.pdb')
	pdbs = glob(search_space)

	# Getting scores file
	base_name = basename(folder.rstrip('/'))
	report_name = join(folder, base_name + '_combined_reports.fasc')
	with open(report_name, 'r') as r:
		all_scores = [line.split() for line in r.readlines()]
	score_head = all_scores.pop(0)
	head_len = len(score_head)

	# Get list of peptide sequences
	sequences = []
	for p in pdbs:
		for i in p.split('_'):
			if 'SMHL' in i:
				if i not in sequences:
					sequences.append(i)

	# Looking at decoys by sequence
	for s in sequences:
		# Getting scores lines for this sequence
		scores_lines = [l for l in all_scores if s in l[0]]

		# Clustering scores lines by mutation profile
		m_profiles = cluster_designs(scores_lines, head_len)
		clusters = []
		for pro, sco in m_profiles.items():
			clusters.append(design_group(pro, [score_head] + sco))
		
		# Writing PyMOL script
		write_pymol(folder, s, clusters)


if __name__ == '__main__':
	args = parse_args()
	main(args.folder)
