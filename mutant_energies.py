#!/usr/bin/python
""" 
PyRosetta4, Python 2.7
Joseph Lubin, 2017
"""
import argparse
from glob import glob
from numpy import mean
from os.path import basename, join, isfile
from pyrosetta import *
from pyrosetta.teaching import EMapVector
from pyrosetta.rosetta import *
from design_protease import *
import cPickle as pickle

def parse_args():
	parser = argparse.ArgumentParser()
	parser.add_argument("folder", help="Pick folder to check")
	parser.add_argument("-des_pep", "--design_peptide", action = "store_true", 
		help="Expand residue selection to peptide, not just protease")
	args = parser.parse_args()
	return args


class mutation_collection:
	"""
	This class object will include information on the mutations in decoys 
	produced by design_protease.py
	"""
	def __init__(self, folder, sequence, design_peptide=False):
		# Identifying decoy location
		self.folder = folder

		# Identifying peptide sequence, and whether it is cleaved
		self.sequence = sequence
		self.seq_and_cleavage()

		# Getting lists of designed and relaxed decoys
		self.collect_models()

		# Defining relevant residues
		self.design_peptide = design_peptide
		self.select_residues()

		# Setting score function
		self.sf = get_fa_scorefxn()

		# Identifying mutations and their frequency and energies
		other_class_values = ['mutations', 'mut_freq', 'mut_location', 
			'mut_res_es', 'mut_int_es']
		for i in other_class_values:
			setattr(self, i, [])

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
		relaxed_search = join(self.folder, '*' + self.sequence + '_relaxed*.pdb*')
		self.relaxed_pdbs = sorted(glob(relaxed_search))

		# Finding designed decoys
		designed_search = join(self.folder, '*' + self.sequence + '_designed*.pdb*')
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
		designable_selector = mutable_residues_selector(self.design_peptide)

		threaded_pose = pose_from_pdb(self.threaded_pdb)
		self.pep_res = selector_to_list(threaded_pose, peptide_selector)
		self.des_res = selector_to_list(threaded_pose, designable_selector)

	##########################################################################
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


class mutations_aggregate:
	def __init__(self, mc_set):
		self.mc_set = mc_set

		# Get full list of mutable residues (differs within the set)
		self.residues_list = self.get_full_residue_list()
		self.mutable_count = len(self.residues_list)

		# Initializing data bins
		self.total_potential_models = 0
		self.actual_potential_models = [0] * self.mutable_count
		self.mutations = []
		while len(self.mutations) < self.mutable_count:
			self.mutations.append({})
		self.mutations_counts = [0] * self.mutable_count

		# Populating bins
		self.starting_residues = self.initial_residues()
		self.collect_mutations_and_frequencies()
		self.count_mutations()

	def get_full_residue_list(self):
		""" 
		Reads through a set of mutation_collection objects and collects a full
		list of all residues that were mutable across the simulations.
		""" 
		complete_residues_list = []
		for mc in self.mc_set:
			for i in mc.des_res:
				if i not in complete_residues_list:
					complete_residues_list.append(i)

		complete_residues_list.sort()
		return complete_residues_list

	def initial_residues(self):
		""" Collect starting residues, before mutation. Sloppy hard-coded. """
		pose_seq = 'GSVVIVGRIILSGRGGPITAYAQQTRGLLGCIITSLTGRDKNQVEGEVQIVSTAAQT'
		pose_seq += 'FLATCINGVCWTVYHGAGTRTIASPKGPVIQMYTNVDQDLVGWPASQGTRSLTPCT'
		pose_seq += 'CGSSDLYLVTRHADVIPVRRRGDSRGSLLSPRPISYLKGSSGGPLLCPAGHAVGIF'
		pose_seq += 'RAAVCTRGVAKAVDFIPVENLETTMRS'

		starting_residues = []
		for i in self.residues_list:
			if i < len(pose_seq):
				starting_residues.append(pose_seq[i-1])

		# Matching length
		while len(starting_residues) < len(self.residues_list):
			starting_residues.append('var')

		return starting_residues

	def mc_name_to_table(self, mc):
		"""
		mutation_collection objects store mutations in the format of N_A_B, 
		where N is the residue number, A is the starting residue, and B is 
		the ending residue. For this application, a list of lists is more 
		useful. This function makes that conversion and appends the frequency 
		of mutation. Returns a list in form [N, B, F], where F is the 
		frequency with which the mutation was observed.
		"""
		converted_mutations_list = [i.split('_') for i in mc.mutations]

		for n, i in enumerate(converted_mutations_list):
			 i[0] = int(i[0])
			 i.pop(1) # Remove original residue
			 i[1] = str(i[1])
			 i.append(mc.mut_freq[n])

		converted_mutations_list.sort()
		return converted_mutations_list

	def collect_mutations_and_frequencies(self):
		""" 
		Reads through the list of mutation_collection objects to extract a 
		list of mutations at each position, and also the total number of 
		decoys both within the set and which have each mutation.
		"""
		for mc in self.mc_set:
			# Getting sum count of all models in all sets
			self.total_potential_models += mc.decoy_count

			# Getting count of models with the potential for mutations, since
			# some model selectors allowed different residues
			for n, i in enumerate(self.residues_list):
				if i in mc.des_res:
					self.actual_potential_models[n] += mc.decoy_count

			# Collecting mutations
			mc_muts = self.mc_name_to_table(mc)
			for i in mc_muts:
				ind = self.residues_list.index(i[0])

				bin_space = self.mutations[ind]
				if i[1] in bin_space:
					bin_space[i[1]] += i[2]
				else:
					bin_space[i[1]] = i[2]


	def count_mutations(self):
		""" 
		Gets count of how many times a residue has mutated across all sets 
		"""
		for n, i in enumerate(self.mutations):
			for m in i:
				self.mutations_counts[n] += i[m]


def aggregated_report(mc_set):
	"""
	Whereas the mutations_by_seq_report lists the mutations each sequence 
	experienced, this report summarizes mutation occurrences across the entire
	available dataset.
	"""
	# Get full list of residues that could have been mutated




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
	
	# Checking if analysis has already been run and pickled
	pickle_name = join(folder, folder.rstrip('/') + '_mutations.pkl')

	if not isfile(pickle_name):
		# Simultaneously analyzing and writing to by-sequence report
		seq_mutants = []
		report_name = join(folder, folder.rstrip('/') + '_mutation_summary.txt')
		for n, i in enumerate(peptide_sequences):
			seq_muts = mutation_collection(folder, i, args.design_peptide)
			seq_mutants.append(seq_muts)

			if i == peptide_sequences[0]:
				mutations_by_seq_report(report_name, seq_muts, head=True)
			else:
				mutations_by_seq_report(report_name, seq_muts)

		# Pickling 
		with open(pickle_name, 'wb') as out:
			pickle.dump(seq_mutants, out, protocol=-1)

	else:
		# Retrieving pickled info
		with open(pickle_name, 'rb') as i:
			seq_mutants = pickle.load(i)

	# Making cross-sequence aggregated report
	

if __name__ == '__main__':
	main()