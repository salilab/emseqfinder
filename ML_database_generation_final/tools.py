import os
import numpy
import sys

def get_sequence_from_pdb(pdb):
    '''
    from a standard pdb file, return all RESNAME entries 
    from the Calpha atoms
    '''
    resnames = []
    f = open(pdb)
    for l in f.readlines():
        #print(l, l[13:15], l[17:20])
        if l[13:15]=="CA":
            resnames.append(l[17:20])

    f.close()
    return resnames

def is_poly_ala(pdb):
    # Returns true if all residues are alanine
    resis = get_sequence_from_pdb(pdb)
    if all(r == 'ALA' for r in resis):
        return True
    else:
        return False

def is_length_correct(pdb):
    # returns
    # pdb : 
    l = int(os.path.basename(pdb).split("_")[-2])
    resis = get_sequence_from_pdb(pdb)
    if len(resis) == l:
        return True
    else:
        return False


def is_length_correct_for_parts(pdb):
    # returns
    # pdb : 
    resis = get_sequence_from_pdb(pdb)
    if len(resis) >2:
        return True
    else:
        return False


def get_adjacent_values(em, voxel):
    # reutrns a list of the value of all voxels
    # plus one
    loc = em.get_location_by_voxel(v)
    values = []
    vs = em.get_header().get_spacing()

    vecs = [(0,0,1), (0,1,0), (1,0,0), (0,0,-1), (0,-1,0), (-1,0,0)]

    disps = [IMP.algebra.Vector3D(v) for v in vecs]

    for d in disps:
        try:
            values.append(em.get_value(loc+d))
        except:
            pass

    return values

def get_voxel_position_string(dmap):
    for v in range(dmap.get_number_of_voxels()):
        vstring += str(dmap.get_location_by_voxel(v))+" "
    return vstring[0:-1]

def get_voxel_value_string(dmap, const=1):
    vstring = ""
    for v in range(dmap.get_number_of_voxels()):
#        vstring += str(int(dmap.get_value(v)*const*10000)/10000)+" "
        vstring += str(round(dmap.get_value(v), 4))+" "
    return vstring[0:-1]

def get_xml_line(xml, field):
    f = open(xml, "r")

    for l in f.readlines():
        if field in l:
            f.close()
            return l
    f.close()
    return None

def get_pdbid_from_xml(xml):
    l = get_xml_line(xml, "<fittedPDBEntryId>")
    if l is None:
        return l
    return l.strip().split("<fittedPDBEntryId>")[-1][0:4]

def get_threshold_from_xml(xml): 
    l = get_xml_line(xml, "contourLevel")
    if l is None:
        return l
    return float(l.strip().split(">")[1].split("<")[0])

def get_resolution(emdb):
    # get resolution information from EMDB
    from config import get_xml
    return get_resolution_from_xml(get_xml(emdb))
    
def get_resolution_from_xml(xml):
    l = get_xml_line(xml,"<resolutionByAuthor>")
    if l is None:
        return l
    return float(l.strip().split("<resolutionByAuthor>")[-1].split("<")[0])

def get_spacing_from_xml(xml):
    l = get_xml_line(xml, "pixelX")
    return float(l.strip().split(">")[1].split("<")[0])
    
def get_statistics_from_xml(xml):
    # Returns voxel minimum, maximum, average and SD as a tuple
    l = get_xml_line(xml, "<minimum>")
    minimum = float(l.strip().split(">")[1].split("<")[0])
    l = get_xml_line(xml, "<maximum>")
    maximum = float(l.strip().split(">")[1].split("<")[0])
    l = get_xml_line(xml, "<average>")
    average = float(l.strip().split(">")[1].split("<")[0])
    l = get_xml_line(xml, "<std>")
    std = float(l.strip().split(">")[1].split("<")[0])
    return (minimum, maximum, average, std)

def catstring(stuff, delimiter=" "):
    outstring = ""
    for s in stuff:
        outstring+=str(s)+delimiter

    return outstring[0:-1]

def tint(value, n=2):
    return int(value*10**n)/10**n
