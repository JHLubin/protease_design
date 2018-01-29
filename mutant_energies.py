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
from design_protease import *

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
		self.sequence = sequence
		self.seq_and_cleavage()

		# Getting lists of designed and relaxed decoys
		self.collect_models()

		# Defining relevant residues
		self.select_residues()

		# Setting score function
		self.sf = get_fa_scorefxn()

		# Identifying mutations and their frequency and energies
		other_class_values = ['mutations', 'mut_freq', 'mut_location', 
			'mut_res_es', 'mut_int_es']
		for i in other_class_values:
			setattr(self, i, [])
		self.mutation_objects = []

		for i in range(self.decoy_count):
			print '\tanalyzing decoy', i
			self.ident_mutations(self.relaxed_pdbs[i], self.designed_pdbs[i])

		# Getting mutation frequencies, energy averages
		self.mut_rate = \
			[round(float(i)/self.decoy_count,3) for i in self.mut_freq]
		self.mut_res_average_e = [mean(i) for i in self.mut_res_es]
		self.mut_int_average_e = [mean(i) for i in self.mut_int_es]

	##########################################################################
	def classify_residues(self):
		""" 
		Looks at the residues in a sequence and determines whether they are 
		positively charged (+), negatively charged (-), neutral-polar (N), or
		hydrophobic (O).
		"""
		res_types = {'A':'O', 'C':'N', 'D':'-', 'E':'-', 'F':'O', 'G':'O', 
			'H':'+', 'I':'O', 'K':'+', 'L':'O', 'M':'O', 'N':'N', 'P':'O', 
			'Q':'N', 'R':'+', 'S':'N', 'T':'N', 'V':'O', 'W':'O', 'Y':'O'}

		self.peptide_res_types = ''
		self.peptide_charge = 0
		self.peptide_polarity = 0

		for res in self.short_sequence:
			res_type = res_types[res]
			self.peptide_res_types += res_type

			if res_type == '+':
				self.peptide_charge += 1
			if res_type == '-':
				self.peptide_charge -= 1

			if res_type in ['+', '-', 'N']:
				self.peptide_polarity += 1

	def seq_and_cleavage(self):
		""" 
		Determines the input peptide sequence and whether it was from a 
		cleaved or uncleaved substrate. The sequences fed into this object
		by isolate_sequence_from_fasc will have either SEQUENCE or unSEQUENCE,
		depending on whether the sequence is cleaved or not.
		"""
		if 'un' in self.sequence:
			self.sequence = self.sequence.replace('un','')
			self.cleaved = 'uncleaved'
		else:
			self.cleaved = 'cleaved'
		print '\n', self.sequence, self.cleaved

		self.short_sequence = self.sequence[1:7]
		self.classify_residues()
		print self.short_sequence
		print self.peptide_res_types
		print "Peptide charge:", self.peptide_charge
		print "Peptide polarity:", self.peptide_polarity

	def collect_models(self):
		""" Collects threaded, designed, and relaxed decoys """
		# Finding the threaded pre-relax model
		thread_search = join(self.folder, '*' + self.sequence + '.pdb.gz')
		threaded_structures = glob(thread_search)
		assert len(threaded_structures) == 1 
			# There should only be one threaded model for each sequence
		self.threaded_pdb = threaded_structures[0]

		# Finding relaxed decoys
		relaxed_search = join(self.folder, '*' + self.sequence + '_relaxed*.pdb')
		self.relaxed_pdbs = sorted(glob(relaxed_search))

		# Finding designed decoys
		designed_search = join(self.folder, '*' + self.sequence + '_designed*.pdb')
		self.designed_pdbs = sorted(glob(designed_search))

		# Counting decoys
		assert len(self.designed_pdbs) == len(self.relaxed_pdbs)
		self.decoy_count = len(self.designed_pdbs)

	def select_residues(self):
		""" 
		Uses residue selectors to determine the set of designable residues.
		mutable_residues_selector and selector_to_list are from 
		design_protease. The option to select peptide residues is passed as 
		true so that any mutations will be collected, regardless of whether 
		the peptide was designable.
		"""
		peptide_selector = ResidueIndexSelector('198-203')
		designable_selector = mutable_residues_selector(True)

		threaded_pose = pose_from_pdb(self.threaded_pdb)
		self.pep_res = selector_to_list(threaded_pose, peptide_selector)
		self.des_res = selector_to_list(threaded_pose, designable_selector)

	##########################################################################
	def track_mutation(self, locus, start, end):
		""" 
		Makes a mutation_location object to keep track of a mutation and its 
		frequency of occurrence.
		"""
		ml = mutation_location(self.short_sequence, self.decoy_count, 
								locus, start, end)
		return ml

	def residue_before_and_after(self, res, before, after):
		""" Determines the AA in the same position in two poses """
		start_res = before.residue(res)
		start_res_name = start_res.name1()

		end_res = after.residue(res)
		end_res_name = end_res.name1()

		mut_string = '_'.join([str(res), start_res_name, end_res_name])

		return start_res_name, end_res_name, mut_string
		
	def interact_energy(self, pose, res, interaction_set):
		""" 
		Takes a given pose, residue, and set of interaction residues to check 
		against, and sums the weighted residue-pair interaction energies 
		between the input residue and each member of the interaction set.
		ref_wts from design_protease.
		"""
		# Collecting interaction energies
		sum_int_energy = 0
		r1 = pose.residue(res)
		for i in interaction_set:
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

	def residue_pairwise_energy(self, res, before, after):
		""" 
		Theis function determines whether a mutation is on the protease or  
		peptide, and conversely the relevant set of residues contributing to 
		the mutated residue's score. Then it runs interact_energy accordingly
		on the residue in the relaxed and designed poses, and calculates the
		difference.
		"""
		if res in self.pep_res:
			mut_loc = 'Peptide'
			partner = self.des_res
		else:
			mut_loc = 'Protease'
			partner = self.pep_res

		# Getting sum residue interaction energy with contact residues
		res_start_int_e = self.interact_energy(before, res, partner)
		res_end_int_e = self.interact_energy(after, res, partner)
		res_delta_int_e = res_end_int_e - res_start_int_e

		return mut_loc, res_delta_int_e

	def ident_mutations(self, rel_pdb, des_pdb):
		"""
		Compares the sequences of a decoy before and after design at specified 
		residues, and identifies the differences. Mutations are identified in 
		the format of N_A_B, where N is the residue number, A is the starting 
		residue, and B is the ending residue. The residue energy difference for 
		the mutation is also collected, as is the sum of two-residue 
		interactions between the mutated residue and the peptide residues (if 
		the mutated residue is on the protease) or between the mutated 
		residue and the mutable protease residues (if the residue is on the 
		peptide). relax_res_energy ia s function from design_protease.
		"""
		# Getting starting and end poses, excerpts from residue energy tables
		relax_pose = pose_from_pdb(rel_pdb)
		design_pose = pose_from_pdb(des_pdb)
		self.sf(relax_pose)
		self.sf(design_pose)

		# Evaluating protease mutations
		for i in self.des_res:
			# Selecting residues to compare
			start_res_name, end_res_name, mut_string = \
				self.residue_before_and_after(i, relax_pose, design_pose)
			
			# If they are different, i.e. there was a mutation
			if start_res_name != end_res_name:
				# Getting overall residue energy change
				relax_res_energy = res_scores(relax_pose, [i], self.sf)[0]
				design_res_energy = res_scores(design_pose, [i], self.sf)[0]
				res_delta_e = design_res_energy - relax_res_energy

				# Getting pairwise residue energy change
				mut_loc, res_delta_int_e = \
					self.residue_pairwise_energy(i, relax_pose, design_pose)

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


class mutation_location:
	""" 
	Stores a mutation's location, before-and-after, and rate of occurrence 
	"""
	def __init__(self, pep_sequence, num_models, locus, initial_aa, final_aa):
		self.peptide_sequence = pep_sequence
		self.mutated_position = locus
		self.native_residue = initial_aa
		self.mutant_residue = final_aa
		self.number_decoys = num_models
		self.observed_occurrences = 0.0

	def mutation_frequency(self):
		""" Returns the rate in the sample set of a mutation's occurrence """
		return self.observed_occurrences / self.number_decoys


def isolate_sequence_from_fasc(fasc_name):
	""" 
	Reads the name of a FASC file output by design_protease.py, isolates and 
	returns the peptide sequence 
	"""
	short_name = basename(fasc_name)
	minus_tail = short_name.replace('_designed.fasc', '')
	sequence = minus_tail.replace('cleaved_ly104_wt_','')
	return sequence


def mutations_by_seq_report(name, mc, head=False):
	"""
	Writes a mutation_summary report, where the mutations for each decoy set
	are summarized.
	"""
	# Formatting report content
	template = '{:10s}' + '{:15s}' + '{:10s}' * 3 + '{:18s}' * 2 + '\n'
	
	# Writing header
	if head:
		header = ['Cleaved', 'Sequence', 'Location', 'Mutation', 'Frequency', 
				'delta_E_residue', 'delta_E_interaction']
		with open(name, 'w') as r:
			r.write(template.format(*header))

	# Writing body
	report_lines = []
	for i in range(len(mc.mutations)):
		line = [str(j) for j in [mc.cleaved, mc.sequence, 
				mc.mut_location[i], mc.mutations[i], mc.mut_rate[i], 
				round(mc.mut_res_average_e[i], 3), 
				round(mc.mut_int_average_e[i], 3)]]
		report_lines.append(line)
	report_lines.sort()

	with open(name, 'a') as r:
		for l in report_lines:
			r.write(template.format(*l))


def aggregated_report(mc_set):
	"""
	Whereas the mutations_by_seq_report lists the mutations each sequence 
	experienced, this report summarizes mutation occurrences across the entire
	available dataset.
	"""

	complete_residues_list = []
	for mc in mc_set:
		for i in mc.des_res:
			if i not in complete_residues_list:
				complete_residues_list.append(i)



##############################################################################
def main():
	# Getting user inputs
	args = parse_args()
	folder = args.folder

	# Intitialize Rosetta
	ros_opts = '-mute all'
	init(ros_opts)	

	# Making mutation_collection objects for each sequence
	report_files = sorted(glob(join(folder, '*.fasc')))
	peptide_sequences = [isolate_sequence_from_fasc(i) for i in report_files]
	
	# Simultaneously analyzing and writing to by-sequence report
	seq_mutants = []
	report_name = join(folder.rstrip('/') + '_mutation_summary.txt')
	for n, i in enumerate(peptide_sequences):
		seq_muts = mutation_collection(folder, i)
		seq_mutants.append(seq_muts)

		if i == peptide_sequences[0]:
			mutations_by_seq_report(report_name, seq_muts, head=True)
		else:
			mutations_by_seq_report(report_name, seq_muts)

	# Making cross-sequence aggregated report
	

if __name__ == '__main__':
	main()