# Example for: selection.mutate()

# This will read a PDB file, change its sequence a little, build new
# coordinates for any of the additional atoms using only the internal
# geometry, and write the mutant PDB file.  It can be seen as primitive
# but rapid comparative modeling for substitution mutants. For more
# sophisticated modeling, see http://salilab.org/modeller/wiki/Mutate%20model
#
# For insertion and deletion mutants, follow the standard comparative
# modeling procedure.

from modeller import Model, Alignment, Selection, Environ
import sys
import os


def mutate_all_ALA(code):
    env = Environ()
    env.io.atom_files_directory = ['../atom_files']

    # Read the topology library with non-hydrogen atoms only:
    env.libs.topology.read(file='$(LIB)/top_heav.lib')
    # To produce a mutant with all hydrogens, uncomment this line:
    # env.libs.topology.read(file='$(LIB)/top_allh.lib')

    # Read the CHARMM parameter library:
    env.libs.parameters.read(file='$(LIB)/par.lib')

    # Read the original PDB file and copy its sequence to the alignment array:
    # pdbpath = sys.argv[1]
    aln = Alignment(env)
    mdl = Model(env, file=code)
    mdl2 = Model(env, file=code)
    aln.append_model(mdl, atom_files=code, align_codes=code)

    # Select the residues to be mutated: in this case all ASP residues:
    # sel = selection(mdl).only_residue_types('ASP')
    sel = Selection(mdl)

    # The second example is commented out; it selects residues '1' and '10'.
    # sel = selection(mdl.residues['1'], mdl.residues['10'])

    # Mutate the selected residues into HIS residues (neutral HIS):
    sel.mutate(residue_type='ALA')

    # Add the mutated sequence to the alignment arrays (it is now the second
    # sequence in the alignment):
    aln.append_model(mdl, align_codes='1fas-1')

    # Generate molecular topology for the mutant:
    mdl.clear_topology()
    mdl.generate_topology(aln['1fas-1'])

    # Transfer all the coordinates you can from the template native structure
    # to the mutant (this works even if the order of atoms in the native PDB
    # file is not standard):
    mdl.transfer_xyz(aln)

    # Build the remaining unknown coordinates for the mutant:
    mdl.build(initialize_xyz=False, build_method='INTERNAL_COORDINATES')
    mdl.res_num_from(mdl2, aln)
    # Write the mutant to a file:
    print(code)
    new_code = os.path.basename(code).split('.')[0]
    mdl.write(file=new_code+'_m.pdb', no_ter=True)


if __name__ == "__main__":
    pdbpath = sys.argv[1]
    mutate_all_ALA(pdbpath)
    print("Done mutation from the main file")
