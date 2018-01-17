#!/usr/bin/env python

from pyrosetta import *
import math, csv, sys, optparse, os, ast, shutil
from pyrosetta.toolbox import *
from transform import *
from SyPRIS_pkgs import *
import rosetta.core.init
import rosetta.core.conformation.symmetry
import rosetta.core.pose.symmetry
import rosetta.core.scoring.symmetry
import rosetta.protocols.simple_moves.symmetry

def add_hetatms(pdb_file, hetatm_lines):
    pdb = read_file(pdb_file)
    new_pdb = pdb + hetatm_lines
    write_file(new_pdb, pdb_file)

def make_flags_file(pdb_file_name,symm_file, res_conversion):
    loop_range = ''
    for pre_ind, residue in res_conversion.iteritems():
        one_range = str(int(residue) -2) + '-' + str(int(residue) +2)
        if not loop_range:
            loop_range += one_range
        else:
            loop_range += ',' + one_range

    with open('%s_chainA.flags' % pdb_file_name[:-4], 'w') as myfile:
        myfile.write('-s ../input_files/%s_chainA.pdb\n-parser:script_vars residue=%s loop_range=%s sym_file=../symmFiles/%s cst_file=../coordFiles/%s_chainA.cst' % (pdb_file_name[:-4], residue, loop_range, symm_file, pdb_file_name[:-4]))

def create_coordCST(pdb_file, res_conversion, hetatm_lines, last_res):
    
    atoms_to_exclude = ['CA', 'CB', 'N', 'O', 'C']
 
    with open(pdb_file, 'r') as myfile:
        pdb = myfile.readlines()

    coordLines = []
    for line in pdb:
        if line[0:6] == 'ATOM  ':
            if line[21] == 'A':
                if line[12:16].strip(' ') not in atoms_to_exclude:
                    coordLine = "CoordinateConstraint " + line[12:16] + " " + res_conversion[line[22:30].strip()] + " CA 1 " + line[30:56] + " HARMONIC 0.0 0.10\n"
                    coordLines.append(coordLine)
    
    if hetatm_lines:
        for line in hetatm_lines:
            coordLine = "CoordinateConstraint " + line[12:16] + " " + str(int(last_res) + 1) + " CA 1 " + line[30:56] + " HARMONIC 0.0 0.10\n"
            coordLines.append(coordLine)
            
    

    with open('%s_chainA.cst' % (pdb_file[:-4]), 'w') as myfile:
        myfile.writelines(coordLines)

def mutate_res(pose_input, rotamer_vals, residue, resType, ncAA):
    print 'res:', residue
    print 'restype:', convert_resname(resType)
    if convert_resname(resType) == None:
        ncAA = rosetta.core.conformation.ResidueFactory.create_residue( ncAA.name_map(resType) )
        pose_input.replace_residue(int(residue), ncAA, True)
    else:
        pose_input = mutate_residue(pose_input, int(residue), '%s' % str(resType[0])) #BPY

    for num, chi in enumerate(rotamer_vals):
        if chi == 'X':
            break
        print 'chi type:', chi
        print 'chi:', num+1
        print 'resind:', int(residue)
        if chi < 0.0:
            chi += 360.
        pose_input.set_chi(num+1, int(residue), float(chi))

    return pose_input

def fast_relax(pose, scorefxn, design_shell, repack_shell, loop_res):
    task_pack = standard_packer_task(pose)
    task_pack.restrict_to_repacking()
    task_pack.temporarily_fix_everything()
    move_map = MoveMap()
    move_map.set_chi(False) #set to fixed rotamer
    move_map.set_bb(False) #set to fixed backbone
    for res in repack_shell + design_shell:
        task_pack.temporarily_set_pack_residue(res, True) #allow residue to repack
        if res in design_shell:
            task_pack.design_residue(res) #allow designable residues to be designable
            move_map.set_bb(res, True) #allow designable residues to also sample backbone position
        move_map.set_chi(res, True) #sample rotamers
    fast_relax_mover = rosetta.protocols.relax.FastRelax(scorefxn, standard_repeats=4)
    fast_relax_mover.set_movemap(move_map)
    monte_carlo_mover = MonteCarlo(pose, scorefxn, 0.8)
    for i in xrange(4):
        fast_relax_mover.apply(pose)
        monte_carlo_mover.boltzmann(pose)
    return
    


def post_SyPRIS_prep(pdb_file_name, scaffold_file, all_red_rot_vals, rotamer_restypes, rotamer_res_inds, symm_file, cof_res_dict, hetatm_lines, last_pdb_res, pose_input, ncAA):
    #need rotamer_vals, pre_resn, rotamer_res all on same index
    print scaffold_file
    pose_from_file(pose_input, scaffold_file) #('chainA_temp.pdb')
    #os.remove('chainA_temp.pdb')
    res_conversion = {}
    for index, pre_resn in enumerate(rotamer_res_inds): 
        residue = pose_input.pdb_info().pdb2pose('A', int(pre_resn))
        print 'rosetta res', residue
        print 'chis to apply', all_red_rot_vals[index]
        print 'restype', rotamer_restypes[index]
        #all_red_rot_vals[index] is the chis with which to apply to new res
        #residue = the pose # from the res number given in rotamer_res_inds
        #rotamer_restypes[index] is the residue type we would like to mutate too
        pose_input = mutate_res(pose_input, all_red_rot_vals[index], residue, rotamer_restypes[index], ncAA)
        #res_conversion will have the chelant model res index as a key and the new pose resi as value
        res_conversion[cof_res_dict[str(pre_resn)]] = str(residue)
    last_res = pose_input.pdb_info().pdb2pose('A', int(last_pdb_res))
    pose_input.dump_pdb('%s_chainA.pdb' % pdb_file_name[:-4])
    #add_hetatms('%s_chainA.pdb' % pdb_file_name[:-4], hetatm_lines) 
    #create_coordCST(pdb_file_name, res_conversion, hetatm_lines, last_res)
    #make_flags_file(pdb_file_name, symm_file, res_conversion)
    
    return 

def main(argv):

    parser = optparse.OptionParser(usage="\n\nTake SyPRIS output and mutate matched residues.")

    parser.add_option('--pdb-path', dest = 'pdb_path',
        help = 'The path to the pdbs listed in the outcsv \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('-o', dest = 'final_out_path',
        help = 'The path to where the output INPUT files should be stored. \n \
                Example: /home/user/my_pdbs/')

    parser.add_option('--rot-lib', dest = 'rot_lib',
        default = '',
        help="Rotamer library dictionary file: ( default = '' )")

    parser.add_option('--cst-path', dest = 'cst_path',
        default = './',
        help="Output cst files path: ( default = '/home/wah49/' )")

    parser.add_option('--flag-path', dest = 'flag_path',
        default = './',
        help="Output flag files path: ( default = '/home/wah49/' )")

    parser.add_option('--params-file', dest = 'params_file',
        default = '',
        help="Specify path to params file when mutating to ncAA.")

    
    (options,args) = parser.parse_args()
   
    ncAA_opt = Pose()
    try:
        ncAA = generate_nonstandard_residue_set( ncAA_opt, ["nBIC.params"] ) #['%s' % options.params_file]  )
    except RuntimeError:
        ncAA = []
    #read off the pdb, and send over
    if options.rot_lib:
        with open(options.rot_lib, 'r') as f:
            s = f.read()
            my_dict = ast.literal_eval(s)
    
    pdb_file_line = options.pdb_path
    pdb_file_line = pdb_file_line.strip('\r\n')
    split_file_line = pdb_file_line.split(',')

    full_scaffold = split_file_line[3]
    pdb_file_name = full_scaffold[:-4] + '_' + '_'.join(split_file_line[1:3]) + '_score_' + split_file_line[-1] + '.pdb'
    chis_to_change = split_file_line[4:-1]

    pdb_obj = Transform(read_file('./new_bis_his_matches/c2_matches/' + pdb_file_name))
    pdb_file_name_split = pdb_file_name.split('_')
    #7 for C_2_up/down
    #4 if the C_2_up wasn't there
    pdb_split_res_types = pdb_file_name_split[ pdb_file_name_split.index('relaxed')+7:\
                                               pdb_file_name_split.index('score') ][::2]
    pdb_split_match_types = pdb_file_name_split[ pdb_file_name_split.index('relaxed')+7:\
                                                 pdb_file_name_split.index('score') ][1::2]
    print '\n'
    print pdb_split_res_types
    print pdb_split_match_types
    print '\n'
    rotamers = []
    rotamer_restypes = []
    rotamer_res_inds = []
    cof_res_dict = {}
    for ind, res_type in enumerate(pdb_split_res_types):
        residue_in_scaff = ''
        chelant_resi = ''
        rot_ind = 0
        for i, chara in enumerate(pdb_split_match_types[ind]):
            try:
                a = int(chara)
                chelant_resi += chara
            except ValueError:
                residue_in_scaff += chara
                if len(residue_in_scaff) == 3:
                    rot_ind = (i+1)
                    break
        print residue_in_scaff
        print pdb_split_match_types[ind]
        pre_rotamer = pdb_split_match_types[ind][rot_ind:]
    
        print pre_rotamer
        rotamer = ''
        bad_indices = []
        for ind, char in enumerate(pre_rotamer):
            print char
            if char in ['a','b']:
                bad_indices.append(ind +1)
            else:
                if ind not in bad_indices:
                    rotamer += char

        rotamers.append(rotamer)
        rotamer_restypes.append(residue_in_scaff)
        rotamer_res_inds.append(res_type[:-3])
        cof_res_dict[res_type[:-3]] = chelant_resi 
    
    if not rotamers:
        sys.exit()
    print 'rotamers:', rotamers
    print 'restypes:', rotamer_restypes
    print 'res indices:', rotamer_res_inds
    print 'cof_res_dict:\n', cof_res_dict
    pdb_obj = Transform(read_file('./new_bis_his_matches/c2_matches/' + pdb_file_name))
    cof_residue_blocks = block_pdb_by_res(read_file('./new_bis_his_matches/c2_matches/' + pdb_file_name))
    cof_residue_blocks = [ x for x in cof_residue_blocks if x[0][17:20].strip().upper() not in rotamer_restypes ]
    hetatm_lines = unblock_pdb_by_res(cof_residue_blocks)
    hetatm_lines = [ x for x in hetatm_lines if x[21] == 'A' ] 
    print 'hetatm_lines', hetatm_lines

    all_red_rot_vals = []
    for indicy, rotamer in enumerate(rotamers):
        chelant_resi_key = rotamer_res_inds[indicy]
        rotamer_vals = get_first_three_chis_of_resnum(copy_transform_object(pdb_obj), [int(cof_res_dict[chelant_resi_key])])
        all_red_rot_vals.append(rotamer_vals[:])
        
    all_red_rot_vals[0] = chis_to_change[:]
    print 'all chis:', all_red_rot_vals
    #print chis_to_change
    #sys.exit()
    scaff_pruned_name = '_'.join(pdb_file_name_split[:(pdb_file_name_split.index('centered')+4)])
    scaffold_file = '/home/wah49/database/pdbs/SCAFF/c2_pdbs/postRelax/' + scaff_pruned_name + '_INPUT.pdb'
    symm_file_name = '_'.join(pdb_file_name_split[:(pdb_file_name_split.index('standard')-1)]) + '.symm'

    scaff_obj = Transform(read_file(scaffold_file))
    last_pdb_res = scaff_obj.get_xyz_info(-1)[-1].split('_')[0]
    
    post_SyPRIS_prep(pdb_file_name, scaffold_file, all_red_rot_vals, rotamer_restypes, rotamer_res_inds, symm_file_name, cof_res_dict, hetatm_lines, last_pdb_res, ncAA_opt, ncAA)

if __name__ == '__main__':
    pyrosetta.init( extra_options='-ignore_zero_occupancy false')# -extra_res /home/wah49/database/params_files/BIC_files/BIC.params' )
    #res_set2 = generate_nonstandard_residue_set( Vector1( ['/home/wah49/params_files/CoO4.params'] ) )
    main(sys.argv[1:])
