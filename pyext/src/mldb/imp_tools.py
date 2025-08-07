import IMP
import IMP.atom
import IMP.core
import IMP.pmi
import IMP.pmi.topology
import numpy
from . import config

alpha = IMP.pmi.alphabets.amino_acid
ff = IMP.atom.CHARMMParameters(IMP.atom.get_data_path("top_heav.lib"),
                               IMP.atom.get_data_path("par.lib"))

m = IMP.Model()
rh = IMP.atom.read_pdb(config.refpdb, m)
res = IMP.atom.Residue(IMP.atom.get_by_type(rh, IMP.atom.RESIDUE_TYPE)[0])
CA_ref_coord = IMP.core.XYZ(
    IMP.atom.get_atom(res, IMP.atom.AT_CA)).get_coordinates()
N_ref_coord = IMP.core.XYZ(
    IMP.atom.get_atom(res, IMP.atom.AT_N)).get_coordinates()
C_ref_coord = IMP.core.XYZ(
    IMP.atom.get_atom(res, IMP.atom.AT_C)).get_coordinates()
O_ref_coord = IMP.core.XYZ(
    IMP.atom.get_atom(res, IMP.atom.AT_O)).get_coordinates()
CB_ref_coord = IMP.core.XYZ(
    IMP.atom.get_atom(res, IMP.atom.AT_CB)).get_coordinates()


def get_backbone_particles(pdb):
    m = IMP.Model()
    rh = IMP.atom.read_pdb(pdb, m, IMP.atom.BackbonePDBSelector())
    return IMP.atom.get_leaves(rh)


def compute_ccc_of_particle_set(dmap, ps, resolution, threshold=2.0):

    # Get a cropped map from the particle set
    dc = dmap.get_cropped(ps, threshold**2)
    dc_t = IMP.em.get_threshold_map(dc, 0.0)
    dc_t.std_normalize()

    smap = IMP.em.create_density_map(IMP.em.get_bounding_box(dc),
                                     dmap.get_header().get_spacing())
    smap_c = IMP.em.CoarseL2Norm.get_density_from_particle(
        smap, ps, resolution)
    smap_tc = IMP.em.get_threshold_map(smap_c, 0.0)

    # Normalize the intensities via histogram matching
    IMP.em.CoarseL2Norm.get_normalized_intensities(smap_tc, dc_t, resolution)
    # mm = IMP.em.approximate_molecular_mass(dc_t, 0.0)
    ccc = IMP.em.CoarseCC()
    try:
        return ccc.cross_correlation_coefficient(dc_t, smap_tc, 0.0, True)
    except:  # noqa: E722
        return -1


def make_reference_density_map():
    bbox = IMP.algebra.BoundingBox3D(
        config.res_bounding_points[0], config.res_bounding_points[1])

    dmap = IMP.em.create_density_map(bbox, config.ref_voxel_size)
    return dmap


def get_voxel_coordinates(dmap):
    voxels = {}
    for v in range(dmap.get_number_of_voxels()):
        point = dmap.get_location_by_voxel(v)
        voxels[v] = (point[0], point[1], point[2])
    return voxels


def mutate_residue(residue_particle, new_aa):
    # Given a residue and a new amino acid, change the hierarchy to
    rph = IMP.atom.Hierarchy(residue_particle)
    res = IMP.atom.Residue(residue_particle)
    atoms = rph.get_children()

    # Delete all non backbone atoms
    for a in atoms:
        if IMP.atom.Atom(a).get_atom_type() not in [IMP.atom.AT_CA,
                                                    IMP.atom.AT_N,
                                                    IMP.atom.AT_O,
                                                    IMP.atom.AT_C]:
            p = a.get_particle()
            rph.remove_child(a)
            residue_particle.get_model().remove_particle(p)

    # Get new residue atoms
    res.set_residue_type(IMP.atom.ResidueType(new_aa))
    res.set_name(new_aa)

    # Use ff to add new residue atoms
    topology = ff.create_topology(rph)
    topology.add_atom_types(rph)
    topology.add_missing_atoms(rph)
    topology.add_coordinates(rph)

    bonds = topology.add_bonds(rph)
    ff.create_angles(bonds)
    ff.create_dihedrals(bonds)
    topology.add_impropers(rph)


def parse_configuration(config_file):

    f = open(config_file, "r")
    config = {}
    for line in f.readlines():
        if len(line.strip()) == 0 or line[0] == "#":
            continue

        if "=" in line:
            key = line.split("=")[0].strip()
            value = line.split("=")[1].strip().split(" ")[0].strip()
            value = value.split("#")[0].strip()

            # If there is a period, it's a float
            if "." in value and "/" not in value:
                config[key] = float(value)

            # Else see if it should be an int
            else:
                try:
                    config[key] = int(value)
                except:  # noqa: E722
                    config[key] = value
    return config


def get_psipred_dictionary_from_files(files):
    # Given a set of files .ss2 with name somthing.A.something
    # where A is the chain name
    # Open each up and return a dictionary with key chain_id and
    # a subdictionary of residue numbers and HEC tuples

    psipred = {}
    for pp in files:

        chain = pp.split(".")[-2]
        print("DHSIUDHSA", pp, chain)
        psipred[chain] = open_ss2_file(pp)
    return psipred


def open_ss2_file(filename):
    f = open(filename, "r")

    outdict = {}
    for line in f.readlines():
        if line.strip() == "" or line[0] == "#":
            continue
        fields = line.split()
        resid = int(fields[0])
        outdict[resid] = (float(fields[4]), float(fields[3]),
                          float(fields[5].strip()))

    return outdict


def compare_sequences_and_return_offset(fasta, query):
    # Sequences should be in dictionary form {resnum:res}
    # Returns offset with an integer offset if they are the same
    # Returns False if they are not the same
    # Gaps are ok

    if compare_sequences(fasta, query):
        return 0

    # Available range of offsets i
    minoffset = 1-min(query.keys())
    maxoffset = max(fasta.keys())-max(query.keys())

    for i in range(minoffset, maxoffset+1):
        if compare_sequences(fasta, query, i):
            return i
    return False


def compare_sequences(seq1, seq2, offset=0):
    minres = numpy.min(list(seq1.keys()))
    maxres = numpy.max(list(seq1.keys()))
    for r in range(minres, maxres+1):
        if r+offset not in seq2.keys():
            pass
        elif seq1[r] != seq2[r+offset]:
            return False
    return True


def read_sequences_from_rcsb(sequence_file):
    sequences = IMP.pmi.topology.Sequences(sequence_file)
    seqs = {}
    for s in sequences.sequences.keys():
        s_cs = s.split("_")[-1].split(" ")[0].strip()
        '''
        if "Chains" in s_cs:
            chains = s_cs.split()[1].split(",")
        else:
            chains = s_cs.split()
        for c in chains:
            seqs[c]=sequences[s]
        '''
        seqs[s_cs] = sequences[s]
    return seqs


def read_sequences(sequence_file):
    sequences = IMP.pmi.topology.Sequences(sequence_file)
    seqs = {}
    for s in sequences.sequences.keys():
        s_cs = s.split("|")[-1].strip()
        if "Chains" in s_cs:
            chains = s_cs.split()[1].split(",")
        else:
            chains = s_cs.split()
        for c in chains:
            seqs[c] = sequences[s]
    return seqs


def get_pdb_offset(pdb, sequences):
    # Find the sequence offset between the pdb file and fasta file
    # The fasta file is always the reference

    # Make an IMP hierarchy and compare that sequence to the fasta
    m = IMP.Model()
    root_hier = IMP.atom.read_pdb(pdb, m, IMP.atom.ATOMPDBSelector())
    pdbseq = {}
    pdb_offsets = {}

    # First, get all the residues in each chain.
    for c in IMP.atom.get_by_type(root_hier, IMP.atom.CHAIN_TYPE):
        chain = IMP.atom.Chain(c).get_id()
        pdbseq[chain] = {}

        # Put PDB residues into a dictionary
        for r in IMP.atom.get_by_type(c, IMP.atom.RESIDUE_TYPE):
            res = IMP.atom.Residue(r)
            rt = res.get_residue_type().get_string()
            pdbseq[chain][res.get_index()] = \
                alpha.get_one_letter_code_from_residue_type(rt)

    # Now, compare the fasta to PDB
    for chain in pdbseq.keys():
        offset = compare_sequences_and_return_offset(
            string_to_dict(sequences[chain]), pdbseq[chain])
        pdb_offsets[chain] = offset

    return pdb_offsets


def string_to_dict(string, start=1):
    seqdict = {}
    for i in range(len(string)):
        seqdict[i+start] = string[i]
    return seqdict


def dict_to_string(seqdict):
    lastres = max(seqdict.keys())
    seq = ""
    for i in range(lastres):
        if i+1 in seqdict.keys():
            if len(seqdict[i+1]) == 3:
                seq += alpha.get_one_letter_code_from_residue_type(
                    seqdict[i+1])
            else:
                seq += seqdict[i+1]
        else:
            seq += "-"
    return seq
