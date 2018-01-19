#!/usr/bin/python
""" 
PyRosetta4, Python 2.7
Joseph Lubin, 2017
"""
import argparse
from glob import glob
from numpy import mean
from os.path import basename, join
from pyrosetta import *
from pyrosetta.teaching import EMapVector
from pyrosetta.rosetta import *

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("folder", help="Pick folder to check")
	args = parser.parse_args()
	return args


class mutation_collection:
	"""
	This class object will include information on the mutations in decoys 
	produced by design_protease.py
	"""
	def __init__(self, folder, sequence):
		# Identifying decoy location
		self.folder = folder

		# Identifying peptide sequence, and whether it is cleaved
		if 'un' in sequence:
			self.sequence = sequence.replace('un','')
			self.cleaved = 'uncleaved'
		else:
			self.sequence = sequence
			self.cleaved = 'cleaved'
		print '\n', self.sequence, self.cleaved

		# Getting lists of designed and relaxed decoys
		self.des_search = join(folder, '*' + self.sequence + '_designed*.pdb')
		self.designed_pdbs = sorted(glob(self.des_search))
		self.rel_search = join(folder, '*' + self.sequence + '_relaxed*')
		self.relaxed_pdbs = sorted(glob(self.rel_search))
		assert len(self.designed_pdbs) == len(self.relaxed_pdbs)
		self.decoy_count = len(self.designed_pdbs)

		# Defining relevant residues
		self.des_res = [52, 56, 57, 58, 70, 73, 124, 138, 150, 151, 
						152, 153, 169, 171, 172, 173, 174, 175, 183]
		self.pep_res = range(197, 208)

		# Setting score function
		self.sf = create_score_function('ref2015')

		# Identifying mutations and their frequency and energies
		self.mutations = []
		self.mut_freq = []
		self.mut_location = []
		self.mut_res_es = []
		self.mut_int_es = []
		for i in range(self.decoy_count):
			print 'analyzing decoy', i
			self.ident_mutations(self.relaxed_pdbs[i], self.designed_pdbs[i])

		# Getting mutation frequencies, energy averages
		self.mut_rate = [round(float(i)/self.decoy_count,3) for i in self.mut_freq]
		self.mut_res_average_e = [mean(i) for i in self.mut_res_es]
		self.mut_int_average_e = [mean(i) for i in self.mut_int_es]

	def ident_mutations(self, rel_pdb, des_pdb):
		"""
		Compares the sequences of a decoy before and after design at specified 
		residues, and identifies the differences. Mutations are identified in 
		the format of N_A_B, where A is the starting residue, N is the residue 
		number, and B is the ending residue. The residue energy difference for 
		the mutation is also reported, as is the sum of two-residue 
		interactions between the mutated residue and the peptide residues (if 
		the mutated residue is on the protease) or between the mutated 
		residue and the mutable protease residues (if the residue is on the 
		peptide).
		"""
		# Getting starting and end poses, excerpts from residue energy tables
		start_pose = pose_from_pdb(rel_pdb)
		end_pose = pose_from_pdb(des_pdb)
		self.sf(start_pose)
		start_e_tab = str(start_pose.energies()).split('\n')
		self.sf(end_pose)
		end_e_tab = str(end_pose.energies()).split('\n')

		# Evaluating protease mutations
		for i in self.des_res + self.pep_res:
			# Selecting residues to compare
			start_res = start_pose.residue(i)
			start_res_name = start_res.name1()
			end_res = end_pose.residue(i)
			end_res_name = end_res.name1()
			
			# If they are different
			if start_res_name != end_res_name:
				# Identifying mutant
				if i < 100:
					l_num = '0' + str(i)
				else:
					l_num = str(i)
				mut_string = '_'.join([l_num, start_res_name, end_res_name])

				# Getting residue energy change
				res_delta_e = self.res_energy_dif(i, start_e_tab, end_e_tab)

				# Getting location, sum residue interaction energy with peptide residues
				if i in self.des_res:
					mut_loc = 'Protease'
					res_start_int_e = self.interact_energy(start_pose, i, self.pep_res)
					res_end_int_e = self.interact_energy(end_pose, i, self.pep_res)
					res_delta_int_e = res_end_int_e - res_start_int_e

				if i in self.pep_res:
					mut_loc = 'Peptide'
					# Getting sum residue interaction energy with peptide residues
					res_start_int_e = self.interact_energy(start_pose, i, self.des_res)
					res_end_int_e = self.interact_energy(end_pose, i, self.des_res)
					res_delta_int_e = res_end_int_e - res_start_int_e

				# Aggregating results
				if mut_string in self.mutations:
					# Mutation has been identified on previous decoy
					mut_ind = self.mutations.index(mut_string)
					assert self.mut_location[mut_ind] == mut_loc
					self.mut_freq[mut_ind] += 1
					self.mut_res_es[mut_ind].append(res_delta_e)
					self.mut_int_es[mut_ind].append(res_delta_int_e)
				else:
					# New mutation or first decoy
					self.mutations.append(mut_string)
					self.mut_location.append(mut_loc)
					self.mut_freq.append(1)
					self.mut_res_es.append([res_delta_e])
					self.mut_int_es.append([res_delta_int_e])

	def res_energy_dif(self, res, start_energies, end_energies):
		""" 
		From a residue number and two pose_energies tables, determines the 
		residue energy difference 
		"""
		res_start_energies = start_energies[res].split()[1:]
		res_end_energies = end_energies[res].split()[1:]
		res_start_tot_energy = sum([float(j) for j in res_start_energies])
		res_end_tot_energy = sum([float(j) for j in res_end_energies])
		res_energy_dif = round(res_end_tot_energy - res_start_tot_energy, 3)

		return res_energy_dif

	def interact_energy(self, pose, res, int_set):
		""" 
		Takes a given pose, residue, and set of interaction residues to check 
		against, and sums the weighted residue-pair interaction energies 
		between the input residue and each member of the interaction set.
		"""
		# Nonzero component weights for REF_2015
		ref_wts = {
			'fa_atr': 1, 'fa_rep': 0.55, 'fa_sol': 1, 'fa_elec': 1,
			'fa_intra_rep': 0.005, 'fa_intra_sol_xover4': 1, 'lk_ball_wtd': 1, 
			'pro_close': 1.25, 'hbond_sr_bb': 1, 'hbond_lr_bb': 1, 
			'hbond_bb_sc': 1, 'hbond_sc': 1, 'dslf_fa13': 1.25, 'omega': 0.4, 
			'fa_dun': 0.7, 'p_aa_pp': 0.6, 'yhh_planarity': 0.625, 'ref': 1, 
			'rama_prepro': 0.45}

		# Collecting interaction energies
		sum_int_energy = 0
		r1 = pose.residue(res)
		for i in int_set:
			r2 = pose.residue(i)

			# Getting emap
			emap = EMapVector()
			self.sf.eval_ci_2b(r1, r2, pose, emap) # Includes many zero terms

			# Cleaning up emap
			org_emp = {}
			emap_nz = str(emap.show_nonzero()).split()
			for j in range(len(emap_nz)/2): # Converting to a dict
				org_emp[emap_nz[2 * j].rstrip(':')] = float(emap_nz[2 * j + 1])

			# Adding weighted scores
			for j in org_emp:
				if j in ref_wts: # Emap contains terms with 0 weight in REF
					sum_int_energy += org_emp[j] * ref_wts[j]

		return sum_int_energy


def isolate_sequence_from_fasc(fasc_name):
	""" 
	Reads the name of a FASC file output by design_protease.py, isolates and 
	returns the peptide sequence 
	"""
	short_name = basename(fasc_name)
	minus_tail = short_name.replace('_designed.fasc', '')
	sequence = minus_tail.replace('cleaved_ly104_wt_','')
	return sequence


def main():
	# Getting user inputs
	args = parse_args()
	folder = args.folder

	# Intitialize Rosetta
	ros_opts = '-mute core -mute protocols -mute basic'
	init(ros_opts)	

	# Formatting report content, writing header
	report_name = join(folder, folder.rstrip('/') + '_mutation_summary.txt')
	template = '{:10s}' + '{:15s}' + '{:10s}' * 3 + '{:18s}' * 2 + '\n'
	head = ['Cleaved', 'Sequence', 'Location', 'Mutation', 'Frequency', 
			'delta_E_residue', 'delta_E_interaction']
	with open(report_name, 'w') as r:
		r.write(template.format(*head))

	# Making mutation_collection objects for each sequence
	report_files = sorted(glob(join(folder, '*.fasc')))
	peptide_sequences = [isolate_sequence_from_fasc(i) for i in report_files]
	seq_mutants = []
	for i in peptide_sequences:
		seq_muts = mutation_collection(folder, i)
		seq_mutants.append(seq_muts)

		# Writing to report
		report_lines = []
		for j in range(len(seq_muts.mutations)):
			line = [str(k) for k in [seq_muts.cleaved, seq_muts.sequence, 
					seq_muts.mut_location[j], seq_muts.mutations[j], seq_muts.mut_rate[j], 
					seq_muts.mut_res_average_e[j], seq_muts.mut_int_average_e[j]]]
			report_lines.append(line)
		report_lines.sort()

		with open(report_name, 'a') as r:
			for l in report_lines:
				r.write(template.format(*l))
	

if __name__ == '__main__':
	main()