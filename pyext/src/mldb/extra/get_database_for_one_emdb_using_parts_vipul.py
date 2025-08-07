import sys
import glob
import os
import IMP
import IMP.em
import IMP.core
import config
import imp_tools
import tools

'''
This script takes in an MRC file and a set of PDB structure_element files and
extracts the density values around the residue as defined in config.py
'''

IMP.set_log_level(IMP.SILENT)

# EMDBID
emdb = sys.argv[1]

# name of the database file to write
database_file = sys.argv[2]

# Name of the EM map and pdb to use (in 0system)
map_file = os.path.join(emdb, "0system", "emdb_normalized_new.map")

# Get all SE pdb files
pdbs = glob.glob(emdb+"/1structure_elements/*.pdb")
print(pdbs)
if len(pdbs) == 0:
    print("No structure elements for this system")
    exit()
try:
    xml = config.get_xml(emdb)  # Get header file
    resolution = tools.get_resolution_from_xml(xml)
except:  # noqa: E722
    # resolution of the map
    resolution = float(sys.argv[3])

# Read the reference pdb
m = IMP.Model()
ref_mh = IMP.atom.read_pdb(config.refpdb, m)

resis = IMP.atom.get_by_type(ref_mh, IMP.atom.RESIDUE_TYPE)
refr = IMP.atom.Residue(resis[0])
psref = [IMP.atom.get_atom(refr, IMP.atom.AT_C),
         IMP.atom.get_atom(refr, IMP.atom.AT_CA),
         IMP.atom.get_atom(refr, IMP.atom.AT_N),
         IMP.atom.get_atom(refr, IMP.atom.AT_CB)]
xyzref = [IMP.core.XYZ(p).get_coordinates() for p in psref]
xyzrefcb = xyzref[-1]

# Make reference density file; the dimension of the bonding box is defined
# in the config.py file
res_bounding_box = IMP.algebra.BoundingBox3D(
    IMP.algebra.Vector3D(config.res_bounding_points[0]),
    IMP.algebra.Vector3D(config.res_bounding_points[1]))

# Open the experimental map
mrw = IMP.em.MRCReaderWriter()
try:
    mrc = IMP.em.read_map(map_file, mrw)
except:  # noqa: E722
    raise Exception("MRC, "+map_file+", cannot be read")

mrc.get_header_writable().set_resolution(resolution)

# resample density if not the same
if mrc.get_header().get_spacing() != config.ref_voxel_size:
    mrc = IMP.em.get_resampled(mrc, config.ref_voxel_size)

##############################################
# Loop over structure elements
of = open(database_file, "w")
for pdb in pdbs:
    # Get STRIDE secondary structure designation
    pdbname = pdb.split('/')[-1].split('.')[0]
    if pdbname.startswith('h'):
        SS = 'H'
    else:
        SS = 'S'

    ilc = tools.is_length_correct_for_parts(pdb)
    # Ignore poly-alanines or if length is not correct
    if not ilc:
        print("SE", os.path.basename(pdb),
              "skipped: | is_length_correct:", ilc)
        continue
    else:
        print("SE", os.path.basename(pdb), "getting residue densities")

    # Open model coordinates
    mh = IMP.atom.read_pdb(pdb, m)
    residues = IMP.atom.get_by_type(mh, IMP.atom.RESIDUE_TYPE)

    skip = False

    # Ensure there are CBs in here, otherwise you get a segfault
    for res in residues:
        ats = [IMP.atom.Atom(a).get_atom_type() for a in res.get_children()]
        print(res, ats)
        if (IMP.atom.AT_CB not in ats
                and IMP.atom.Residue(res).get_residue_type().get_string()
                != "GLY"):
            skip = True

    if skip:
        continue

    # Cycle over all of the residues in the pdb
    for r in residues:
        res = IMP.atom.Residue(r)
        rid = res.get_index()
        chain = IMP.atom.Chain(r.get_parent()).get_id()
        resname = res.get_residue_type().get_string()

        # Mutate GLY to ALA so we have a CB atom to superimpose
        if resname == "GLY":
            imp_tools.mutate_residue(res, "ALA")

        C_atom = IMP.atom.get_atom(res, IMP.atom.AT_C)
        CA_atom = IMP.atom.get_atom(res, IMP.atom.AT_CA)
        N_atom = IMP.atom.get_atom(res, IMP.atom.AT_N)

        psx = [C_atom, CA_atom, N_atom, IMP.atom.get_atom(res, IMP.atom.AT_CB)]

        psbb = [C_atom, CA_atom, N_atom, IMP.atom.get_atom(res, IMP.atom.AT_O)]

        xyzx = [IMP.core.XYZ(p).get_coordinates() for p in psx]
        # Align backbone atoms
        xform1 = IMP.algebra.get_transformation_aligning_first_to_second(
            xyzx, xyzref)
        IMP.atom.transform(mh, xform1)
        m.update()
        xyzx = [IMP.core.XYZ(p).get_coordinates() for p in psx]

        # Superimpose CA residues
        xform2 = IMP.algebra.Transformation3D(xyzref[1]-xyzx[1])
        IMP.atom.transform(mh, xform2)
        m.update()

        xform_tot = xform1*xform2

        res_map = IMP.em.create_density_map(res_bounding_box,
                                            config.ref_voxel_size)
        for v in range(res_map.get_number_of_voxels()):
            point = res_map.get_location_by_voxel(v)

            # Find where this point in the EM map via xform
            xf_point = xform_tot.get_inverse().get_transformed(point)
            # Get the density at this point
            val = IMP.em.get_density(mrc, xf_point)

            # Set the new map value to this
            res_map.set_value(point[0], point[1], point[2], val)

        # Move the structure back to where it came from
        IMP.atom.transform(mh, xform_tot.get_inverse())
        # Get sum and max of densities
        sc_string = tools.get_voxel_value_string(res_map)
        print("current PDB name", pdbname)
        outstring = tools.catstring(
            [emdb, resolution, pdbname, chain, rid,
             resname, SS]) + " " + sc_string
        of.write(outstring+"\n")

of.close()
