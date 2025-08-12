# emdb_home is where you initially download the map, fasta file, pdb file,
# and xml file, after creating a folder with folder name same as EMDBid
emdb_home = "path"

# reference map for voxel data extraction
refpdb = "./reference/ref.pdb"
res_bounding_points = [(-4, -10, -10), (10, 4, 4)]
ref_voxel_size = 1.0


def get_xml(emdb):
    import os
    return os.path.join(emdb_home, str(emdb), "xml")
