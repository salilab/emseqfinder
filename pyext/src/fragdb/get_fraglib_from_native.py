import os
import IMP.atom
import IMP.core
import argparse


def get_SSE_fragments(stride_file):
    '''
    get all the SSE fragments info for the native structure from its
    stride file

    '''
    current_chain = None
    stride_dict = {}
    frag_list = []
    # 1. collect all helix and strand info as a dictionary, where the
    # information is stored based on chain id as well
    with open(stride_file, 'r') as infile:
        for line in infile:
            if line.startswith('ASG '):
                chain_id = line[9]
                resid = int(line[10:15].strip())
                sstype = line[24]

                if current_chain is None or chain_id != current_chain:
                    current_chain = chain_id
                    stride_dict[current_chain] = {}
                    if sstype == 'H' or sstype == 'E':
                        stride_dict[current_chain][resid] = sstype
                else:
                    if sstype == 'H' or sstype == 'E':
                        stride_dict[current_chain][resid] = sstype

    # 2, Loop over each chain to detect the fragments
    for chainid in stride_dict:
        previous_ss = None
        previous_res = None
        current_frag = []
        for resid in stride_dict[chainid]:
            current_ss = stride_dict[chainid][resid]
            if previous_ss is None:
                previous_ss = current_ss
                previous_res = resid
            # only contiguous residues with same SS type should be
            # part of a fragment
            if current_ss == previous_ss and resid - previous_res < 2:
                current_frag.append((resid, current_ss))
                previous_ss = current_ss
                previous_res = resid
            else:
                current_frag.sort(key=lambda y: y[0])
                first_res = current_frag[0][0]
                last_res = current_frag[-1][0]
                length = len(current_frag)
                new_key = (chainid + '_' + current_frag[0][1] + '_'
                           + str(first_res) + '_' + str(last_res) + '_'
                           + str(length))
                frag_list.append(new_key)
                current_frag = [(resid, current_ss)]
                previous_ss = current_ss
                previous_res = resid
        # This step is to collect the last fragment information for each chain
        if len(current_frag) > 1:
            current_frag.sort(key=lambda y: y[0])
            first_res = current_frag[0][0]
            last_res = current_frag[-1][0]
            length = len(current_frag)
            new_key = (chainid + '_' + current_frag[0][1] + '_'
                       + str(first_res) + '_' + str(last_res) + '_'
                       + str(length))
            frag_list.append(new_key)
    # final fragment list
    return frag_list, stride_dict


def get_fragment_coords(pdbname, frag_info, path_to_store_parts):
    '''
    get coords from pdbname chain id and residue ranges
    also save the pdbs of the helical fragments and check if for strands
    there could be multiple strands possibilities
    if 2 strands are possible then save 2-strand pdbs
    '''
    all_pdbs = []
    m = IMP.Model()
    native_h = IMP.atom.read_pdb(pdbname, m, (IMP.atom.ATOMPDBSelector()))
    outdir = path_to_store_parts
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    selected_atom_types = ['N', 'CA', 'C', 'O', 'CB']
    for frag in frag_info:
        # A_E_5_8_4
        chainid, SSEtype, firstres, lastres, length = frag.split('_')
        if int(length) >= 3:
            if SSEtype == 'H':
                # write a pdb now
                h_frag = IMP.atom.Selection(
                    native_h, chain_id=chainid,
                    residue_indexes=range(int(firstres), int(lastres)+1),
                    atom_types=[IMP.atom.AtomType(n)
                                for n in selected_atom_types])
                pdb_name = '_'.join(['h', length, chainid,
                                     firstres, lastres]) + '.pdb'
                IMP.atom.write_pdb(h_frag, outdir + '/' + pdb_name)
                all_pdbs.append(outdir + '/' + pdb_name)
            elif SSEtype == 'E':
                # collect the coordinates first
                s_frag = IMP.atom.Selection(
                    native_h, chain_id=chainid,
                    residue_indexes=range(int(firstres), int(lastres)+1),
                    atom_types=[IMP.atom.AtomType(n)
                                for n in selected_atom_types])
                pdb_name = '_'.join(['s', length, chainid,
                                     firstres, lastres]) + '.pdb'
                IMP.atom.write_pdb(s_frag, outdir + '/' + pdb_name)
                all_pdbs.append(outdir + '/' + pdb_name)
            else:
                print('Something is wrong, should only process strands '
                      'and helix')
    return all_pdbs


def per_chain_segment_coords(pdbname, segment_info, path_to_store_parts):
    '''
    get coords from pdbname chain id and residue ranges
    '''
    all_pdbs = []
    m = IMP.Model()
    native_h = IMP.atom.read_pdb(pdbname, m, (IMP.atom.ATOMPDBSelector()))
    outdir = path_to_store_parts
    if not os.path.exists(outdir):
        os.makedirs(outdir)
    selected_atom_types = ['N', 'CA', 'C', 'O', 'CB']
    for chains in segment_info.keys():
        all_residueids = list(segment_info[chains].keys())
        print(all_residueids)
        # write a pdb now
        h_frag = IMP.atom.Selection(
            native_h, chain_id=chains,
            residue_indexes=all_residueids,
            atom_types=[IMP.atom.AtomType(n) for n in selected_atom_types])
        pdb_name = '_'.join(['all', str(len(all_residueids)), chains,
                             str(all_residueids[0]),
                             str(all_residueids[-1])]) + '.pdb'
        IMP.atom.write_pdb(h_frag, outdir + '/' + pdb_name)
        all_pdbs.append(outdir + '/' + pdb_name)
    return all_pdbs


def mutate_frag_to_ala(pdbname):
    '''
    mutate all residues from a pdb to ALA
    '''
    import mutate_all_ALA
    mutate_all_ALA.mutate_all_ALA(pdbname)
#    os.remove(pdbname)


parser = argparse.ArgumentParser(
    description='Generate a library of parts from native structure, '
                'using the information from STRIDE')
parser.add_argument('ref_pdb', type=str,
                    help='Provide the absolute path of the native pdb '
                         '(we assumed that as reference PDB)')
parser.add_argument('ref_stride', type=str,
                    help='Provide the absolute path of the stride file')
parser.add_argument('path_to_store_parts', type=str,
                    help='path to store the generated parts')

parser.add_argument('--perChain', type=bool, default=False,
                    help='get segments per chain')
args = parser.parse_args()

all_frags, stride_dict = get_SSE_fragments(args.ref_stride)
print(stride_dict)

if not args.perChain:
    all_new_pdbs = get_fragment_coords(args.ref_pdb, all_frags,
                                       args.path_to_store_parts)
else:
    all_new_pdbs = per_chain_segment_coords(args.ref_pdb, stride_dict,
                                            args.path_to_store_parts)

print(all_frags)
