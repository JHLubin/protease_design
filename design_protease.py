#!/usr/bin/python
""" 
PyRosetta4, Python 2.7
Joseph Lubin, 2017
"""
import argparse
from math import sqrt
from os import makedirs
from os.path import basename, isdir, join
from pyrosetta import *
from pyrosetta.rosetta import *
from pyrosetta.rosetta.core.kinematics import FoldTree
from pyrosetta.rosetta.core.select.residue_selector import \
	AndResidueSelector, ChainSelector, InterGroupInterfaceByVectorSelector,\
	NotResidueSelector, OrResidueSelector, ResidueIndexSelector
from pyrosetta.rosetta.core.pack.task import parse_resfile, TaskFactory
from pyrosetta.rosetta.core.pack.task.operation import \
	OperateOnResidueSubset, PreventRepackingRLT, RestrictToRepackingRLT
from pyrosetta.rosetta.numeric import xyzVector_double_t
from pyrosetta.rosetta.protocols.denovo_design.movers import FastDesign
from pyrosetta.rosetta.protocols.enzdes import ADD_NEW, AddOrRemoveMatchCsts
from pyrosetta.rosetta.protocols.flexpep_docking import FlexPepDockingProtocol
from pyrosetta.rosetta.protocols.relax import FastRelax
from pyrosetta.teaching import MinMover, PackRotamersMover, SimpleThreadingMover
import sys

def parse_args():
	info = "Design a protease around a peptide sequence"
	parser = argparse.ArgumentParser(description=info)
	parser.add_argument("-s", "--start_struct", required=True,
		help="Pick starting PDB")
	parser.add_argument("-od", "--out_dir", required=True,
		help="Name an output directory for decoys")
	parser.add_argument("-cseq", "--cut_peptide_sequence", type=str, 
		action='append', help="List cleaved peptide sequences or provide \
		list file")
	parser.add_argument("-useq", "--uncut_peptide_sequence", type=str, 
		action='append', help="List uncleaved peptide sequences or provide \
		list file")
	parser.add_argument("-cr", "--cat_res", type=int, nargs='+', 
		default=[72, 96, 154], help="The catalytic residues of the protease, \
		excluded from design (defaults are 72, 96, and 154, for HCV)")
	parser.add_argument("-cons", "--constraints", type=str, 
		default='ly104.cst', help="Pick constraints file")
	parser.add_argument("-rad", "--design_rad", type=int, default=8, 
		help="Cutoff for designable residues (Angstroms from peptide, \
		default is 8A)")
	parser.add_argument("-res", "--resfile", type=str,
		help="Pick resfile for design")
	parser.add_argument("-th", "--thread", action = "store_true", 
		help="Option to create threaded models. Use the first time running.")
	args = parser.parse_args()
	return args


def init_opts(cst_file='ly104.cst'):
	""" Produces a list of init options for PyRosetta, including cst file """
	ros_opts = '-mute core -mute protocols -mute basic'
	ros_opts += ' -enzdes::cstfile ' + cst_file
	ros_opts += ' -cst_fa_weight 1.0 -run:preserve_header -out:pdb_gz'
	return ros_opts


def res_ca_cords(pose, res):
	""" Returns the x,y,z coordinates of the A-carbon of a given residue"""
	res_coords = []
	for i in range(3):
		res_coords.append(pose.residue(res).xyz('CA')[i])

	return res_coords


def point_dist(point_1, point_2):
	""" Given two sets of coordinates, determines distance between them """
	sum_difs_squared = 0
	for i in range(3):
		sum_difs_squared += (point_2[i] - point_1[i]) ** 2

	return sqrt(sum_difs_squared)


def res_to_design(pdb, radius=8, exclude_res=[72, 154]):
	""" 
	Determines the chain in the PDB that is the peptide, assuming that the 
	peptide is smallest. Determines the coordinates of all a-carbons in the 
	system, then checks the distances to each atom in the peptide. If the 
	atom is within the cutoff radius, it is added to the list of designable 
	residues. Returns a list of mutable residues.
	"""
	pose = pose_from_pdb(pdb)
	chains = pose.split_by_chain()

	# Determining peptide chain
	pep_chain_no = 1
	for i in range(2, len(chains) + 1):
		if len(chains[i]) < len(chains[pep_chain_no]):
			pep_chain_no = i
	pep_chain = chains[pep_chain_no]
	chains.pop(pep_chain_no) # removes peptide chain from chain list

	# Getting residue number of peptide start
	pep_start = 1
	for i in range(1,pep_chain_no):
		pep_start += chains[i].total_residue()

	# Getting peptide residue CA coordinates:
	pep_coords = []
	for i in range(1, pep_chain.total_residue() + 1):
		pep_coords.append(res_ca_cords(pep_chain, i))

	# Populating the list of designable residues
	mutable_residues = []
	for chain in chains:
		# Note that chains now excludes the peptide chain
		for res in range(1, chain.total_residue() + 1):
			if res in exclude_res:
				# Exclude the catalytic residues from the designable list
				continue
			a_cords = res_ca_cords(pose, res)
			for i in pep_coords:
				# If any res is within radius of any pep res, add to the list
				if point_dist(a_cords, i) <= radius:
					mutable_residues.append(res)
					break

	return mutable_residues  


def selector_to_list(pose, selector):
	""" Converts a selector output vector to a list of selected residues """
	selection_vector = selector.apply(pose)
	selection_list = []
	for i in range(len(selection_vector)): 
		if selection_vector[i+1]==1:
			selection_list.append(i+1)

	return selection_list


def get_seq_list(seq_arg):
	"""
	Takes an argument that can include individual peptide sequences or file(s)
	containing a list of sequences, and returns a list of sequences. 
	Distinguishes between files and sequences by the presence of a dot (.).
	"""
	pep_sequences = []
	for inp in seq_arg:
		if '.' in inp:
			# If input is a file
			with open(inp, 'r') as t:
				lis = t.readlines()
			if len(lis) == 1:
				# If all sequences are listed horizontally on one line
				# rather than one per line, rearrange
				lis = lis[0].split()

			for i in lis:
				pep_sequences.append(i.strip())

		else:
			# Sequence was typed directly into the argument
			pep_sequences.append(inp.strip())

	return pep_sequences
		

def thread_seq(pose, pep_start, pep_length, seq):
	""" Thread a new sequence in for the peptide. """
	tm = SimpleThreadingMover(seq, pep_start)
	tm.apply(pose)
	return pose


def quick_thread(destination, pdb, sequences, cleaved=False, make=False):
	""" 
	Threads a set of sequences onto the peptide portion of the given PDB file,
	outputting a threaded PDB file for each sequence.
	Function is presently hard-coded for this application.
	"""
	pose = pose_from_pdb(pdb)

	thread_files = []

	for seq in sequences:
		# Naming model
		if cleaved:
			pdbname = 'cleaved_ly104_wt_' + seq + '.pdb.gz'
		else:
			pdbname = 'uncleaved_ly104_wt_' + seq + '.pdb.gz'
		out_name = join(destination, pdbname)
		thread_files.append(out_name)

		if make:
			# Threading peptide sequences
			threaded_pose = Pose()
			threaded_pose.assign(pose)
			threaded_pose = thread_seq(threaded_pose, 197, 11, seq)
			threaded_pose.dump_pdb(out_name)

	return thread_files


def mutable_residues_selector():
	"""
	Selects the residues in a shell around the peptide using the 
	InterGroupInterfaceByVectorSelector residue selector
	Presently hard coded for HCV protease.
	Decided to design just around peptide, not pep + cat
	"""
	# Making positive residue selector
	rs = InterGroupInterfaceByVectorSelector()
	rs.group1_selector(ChainSelector("A")) # Protease
	rs.group2_selector(ChainSelector("B")) # Peptide
	rs.nearby_atom_cut(6)
	rs.vector_dist_cut(8)

	# Setting up exclusion of catalytic and peptide residues
	limit_selection = AndResidueSelector()
	not_cat = NotResidueSelector(ResidueIndexSelector('72,96,154')) # Catalytic
	limit_selection.add_residue_selector(not_cat)
	limit_selection.add_residue_selector(ChainSelector("A")) # Exclude peptide
	limit_selection.add_residue_selector(rs)

	return limit_selection


def packable_residues_selector(mutable_selector):
	"""
	Selects the shell of neighbor residues to repack
	Presently hard coded for HCV protease.
	"""
	# Making positive residue selector
	rs = InterGroupInterfaceByVectorSelector()
	rs.nearby_atom_cut(4)
	rs.vector_dist_cut(4)
	rs.group1_selector(NotResidueSelector(mutable_selector))
	rs.group2_selector(mutable_selector)

	# Setting up exclusion of catalytic and mutable residues
	limit_selection = AndResidueSelector()
	not_cat = NotResidueSelector(ResidueIndexSelector('72,96,154')) # Catalytic
	limit_selection.add_residue_selector(not_cat)
	limit_selection.add_residue_selector(rs)
	limit_selection.add_residue_selector(NotResidueSelector(mutable_selector))

	# Add back in the peptide
	expand_selection = OrResidueSelector()
	expand_selection.add_residue_selector(limit_selection)
	expand_selection.add_residue_selector(ChainSelector("B"))

	return expand_selection


def other_residues_selector(mutable_selector, packable_selector):
	""" Selects the residues that are not designable or repackable """
	all_rel_res_sel = OrResidueSelector()
	all_rel_res_sel.add_residue_selector(mutable_selector)
	all_rel_res_sel.add_residue_selector(packable_selector)

	other_res_selector = NotResidueSelector(all_rel_res_sel)

	return other_res_selector


def apply_constraints(pose):
	""" Applies the constraints form the input CST file to a pose """
	cstm = AddOrRemoveMatchCsts()
	cstm.set_cst_action(ADD_NEW)
	cstm.apply(pose)
	return pose


def make_move_map(near_res, pep_start=197, pep_end=208):
	""" 
	Makes a movemap for a protease-peptide system, with all non-peptide 
	residue backbones fixed, and side chains mobile for the peptide and all
	residues in an input list, which is intended to be the nearby residues 
	(8A by default), excluding the catalytic ones. 
	"""
	mm = MoveMap()
	mm.set_bb_true_range(pep_start,pep_end)
	mm.set_chi_true_range(pep_start,pep_end)
	for i in near_res:
		mm.set_chi(i, True)
	
	return mm


def make_fold_tree():
	"""
	Make a fold tree that connects the first catalytic residue to the upstream
	cleaved residue.
	Presently hard-coded for HCV protease
	"""
	ft = FoldTree()
	ft.add_edge(72, 1, -1)
	ft.add_edge(72, 196, -1)
	ft.add_edge(72, 203, 1)
	ft.add_edge(203 ,197, -1)
	ft.add_edge(203 ,207, -1)
	assert ft.check_fold_tree()

	return ft


def make_task_factory(repack_set, other_set):
	""" 
	Makes a TaskFactory with operations that leave the mutable residues 
	designable, restricts the nearby residues to repacking, and prevents 
	repacking of other residues.
	"""
	prevent=PreventRepackingRLT()
	restrict=RestrictToRepackingRLT()

	tf = TaskFactory()
	tf.push_back(OperateOnResidueSubset(prevent, other_set))
	tf.push_back(OperateOnResidueSubset(restrict, repack_set))
	# Everything else left designable

	return tf


def make_pack_task(pose, resfile=None, pack_res=[]):
	""" 
	Makes a packer task for a given pose using an input resfile or list of 
	packable (not designable) residues. 
	"""
	# Packer for protease + peptide\
	task = standard_packer_task(pose)
	if resfile:
		parse_resfile(pose, task, resfile)
	else:
		task.restrict_to_repacking()
		task.temporarily_fix_everything()
		for i in pack_res:
			task.temporarily_set_pack_residue(i, True)

	return task


def fastrelax(pose, score_function, movemap):
	""" 
	Runs the FastRelax protocol on a pose, using given score function and 
	movemap
	"""
	relax = FastRelax()
	relax.set_scorefxn(score_function)
	relax.set_movemap(movemap)

	relax.apply(pose)
	return pose


def fastdesign(pose, score_function, movemap, taskfactory):
	fd = FastDesign()
	fd.set_scorefxn(score_function)
	fd.set_movemap(movemap)
	fd.set_task_factory(taskfactory)

	fd.apply(pose)
	return pose


def minmover(pose, score_function, movemap):
	""" 
	Runs a gradient-base minimization on a pose, using given score function
	and movemap
	"""	
	min_mover = MinMover()
	min_mover.score_function(score_function)
	min_mover.movemap(movemap)

	min_mover.apply(pose)
	return pose	


def design_pack(pose, score_function, task):
	""" Runs packing mover on a pose, using given score function and task """
	pack = PackRotamersMover(score_function, task)
	pack.apply(pose)
	return pose


def res_scores(pose, residues, score_function):
	""" Gets the residue scores for a given set of residues in a pose """
	score_function(pose)
	pose_energies = str(pose.energies()).split('\n') # Table of residue
		# score components, including a header line, so index matches res
	energies_set = []
	for i in residues:
		res_energies = pose_energies[i].split()[1:]
		res_tot_energy = sum([float(j) for j in res_energies])
		energies_set.append(res_tot_energy)

	set_energy = sum(energies_set)
	return set_energy, energies_set


def move_apart(pose, peptide_start, peptide_end):
	""" Moves the peptide a long distance away from the protease """
	# Making displacement vector
	xyz = xyzVector_double_t()
	xyz.x, xyz.y, xyz.z = [100 for i in range(3)]

	# Moving peptide, atom by atom
	for res in range(peptide_start, peptide_end + 1):
		for atom in range(1, pose.residue(res).natoms() + 1):
			pose.residue(res).set_xyz(atom, pose.residue(res).xyz(atom) + xyz)

	return pose


def score_ddg(pose, near_res):
	"""
	Gets a score for the ddG of peptide binding to the protease. This is 
	achieved by taking the peptide and moving each atom a set arbitrary length 
	that is large enough to be far from the protease, then repacking the side
	chains of both the peptide and the mutable protease residues. This 
	function does not take a scorefunction as an input, scoring instead with 
	the default function to ignore the effects of constraints.
	"""
	# Making a new pose to avoid messing up the input
	ddg_pose = Pose()
	ddg_pose.assign(pose)

	sf = get_fa_scorefxn()

	# Score when docked
	dock_score = sf(ddg_pose)

	# Score when separate
	pt = make_pack_task(ddg_pose, pack_res=near_res+range(197,208))
	ddg_pose = move_apart(ddg_pose, 197, 207)
	ddg_pose = design_pack(ddg_pose, sf, pt)
	split_score = sf(ddg_pose)

	ddg = dock_score - split_score

	return [round(i,3) for i in [dock_score, split_score, ddg]]


def ident_mutations(start_pose, end_pose, residues, start_set, end_set):
	"""
	Compares the sequences of a starting pose and ending pose at specified 
	residues and identifies differences. Returns a string listing changes in 
	the format of ANB, where A is the starting residue, N is the residue 
	number, and B is the ending residue.
	"""
	mutations = ''
	for i in residues:
		start_res = start_pose.residue(i).name1()
		end_res = end_pose.residue(i).name1()
		if start_res != end_res:
			r_i = residues.index(i)
			e_dif = '(' + str(round(end_set[r_i] - start_set[r_i], 3)) + ')'
			mut_string = start_res + str(i) + end_res + e_dif
			mutations = ','.join([mutations, mut_string])

	if mutations == '':
		return "NONE"
	else:
		return mutations.lstrip(',')


def set_design(pdb, sf, pep_res, des_res, near_res, num_decoys, resfile, tf):
	"""
	Uses the job distributor to output a set of proteases designed for 
	compatability with a threaded peptide sequence. Outputs a provided number
	of decoys into a directory that is also a required input. Will relax the 
	pdb, then run 20 rounds of design/repacking plus minimization. For all 
	movers, only the residues in the peptide and those within the input list
	are repackable, and only those in the input list are designable. For the 
	relax, the peptide backbone is flexible, and constraints are applied.
	"""
	pose = apply_constraints(pose_from_pdb(pdb))
	ft = make_fold_tree()
	pose.fold_tree(ft)
	#mm = make_move_map(des_res + near_res) # for relax and minimization
	mm = make_move_map(near_res) # for relax and minimization

	dec_name = pdb.replace('.pdb.gz', '_designed')
	jd = PyJobDistributor(dec_name, num_decoys, sf)

	while not jd.job_complete:
		pp = Pose()
		pp.assign(pose)
		# Relaxing
		relax_name = jd.current_name.replace('designed', 'relaxed')
		pp = fastrelax(pp, sf, mm)
		pp.dump_pdb(relax_name)
		rel_pro_set = res_scores(pp, des_res, sf)[1]
		rel_pep_set = res_scores(pp, pep_res, sf)[1]

		# Doing design
		#for i in range(20):
		#	pt = make_pack_task(pp, resfile=resfile)
		#	pp = design_pack(pp, sf, pt)
		#	pp = minmover(pp, sf, mm)
		pp = fastdesign(pp, sf, mm, tf)

		# Getting residue scores, ddG, and mutations list
		prot_res_e, pro_re_set = res_scores(pp, des_res, sf)
		pep_res_e, pep_re_set = res_scores(pp, pep_res, sf)
		dock_split_ddg = score_ddg(pp, des_res + near_res)
		pro_mut = ident_mutations(pose, pp, des_res, rel_pro_set, pro_re_set)
		pep_mut = ident_mutations(pose, pp, pep_res, rel_pep_set, pep_re_set)

		# Making line to add to fasc file
		scores = [prot_res_e, pep_res_e] + dock_split_ddg + [pro_mut, pep_mut]
		temp = "protease_res_scores: {}\tpeptide_res_scores: {}\t"
		temp += "docked_score: {}\tsplit_score: {}\tddG: {}\t"
		temp += "prot_mutations: {}\tpep_mutations: {}"
		score_text = temp.format(*[str(i) for i in scores])
		print score_text, '\n'
		jd.additional_decoy_info = score_text

		jd.output_decoy(pp)


def main():
	# Getting user inputs
	args = parse_args()

	# Initializing PyRosetta
	ros_opts = init_opts(cst_file=args.constraints)
	init(options=ros_opts)

	# Score function
	sf = create_score_function('ref2015_cst')

	# Destination folder for PDB files
	pdb = args.start_struct
	source_pose = pose_from_pdb(pdb)
	dir_nam = args.out_dir
	if not isdir(dir_nam):
		makedirs(dir_nam)

	# Reading inputs for peptide sequences
	cut_seq = get_seq_list(args.cut_peptide_sequence)
	uncut_seq = get_seq_list(args.uncut_peptide_sequence)

	# Creating threaded structures
	make = False
	if args.thread:
		make = True
	t_structs = quick_thread(dir_nam, pdb, cut_seq, cleaved=True, make=make)
	t_structs += quick_thread(dir_nam, pdb, uncut_seq, make=make)

	# Determining peptide part of PDB file, residues near peptide
	pep_res = range(197,208) # Hard code for HCV protease, should update args
	in_shell = args.design_rad
	out_shell = in_shell + 2 # Hard coded 2A, should make a new arg
	cat_res = args.cat_res # Default is for HCV protease
	#des_res = res_to_design(pdb, radius=in_shell, exclude_res=cat_res)
	#inner_res = des_res+cat_res
	#near_res = res_to_design(pdb, radius=out_shell, exclude_res=inner_res)

	# Making residue selectors
	cat_res_selector = ResidueIndexSelector('72,96,154')
	mut_res_selector = mutable_residues_selector()
	pac_res_selector = packable_residues_selector(mut_res_selector)
	oth_res_selector = other_residues_selector(mut_res_selector, 
												pac_res_selector)

	# Converting selectors to lists
	pose = pose_from_pdb(pdb)
	des_res = selector_to_list(pose, mut_res_selector)
	near_res = selector_to_list(pose, pac_res_selector)
	other_res = selector_to_list(pose, oth_res_selector)

	# Making task factory
	tf = make_task_factory(pac_res_selector, oth_res_selector)
	task=tf.create_task_and_apply_taskoperations(pose)

	# Doing design on threaded models
	for struc in t_structs:
		set_design(struc, sf, pep_res, des_res, near_res, 10, args.resfile, tf)


if __name__ == '__main__':
	main()