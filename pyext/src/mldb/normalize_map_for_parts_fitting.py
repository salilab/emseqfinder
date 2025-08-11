# Script to perform normalization on an EM map
# by histogram matching to a reference

import IMP
import IMP.em
import numpy as np
from statsmodels.distributions.empirical_distribution import ECDF
from statsmodels.distributions.empirical_distribution import StepFunction
import argparse
import os
from . import config
from . import tools


def compute_cdf(dmap, exclude_zero=True):
    # Returns an empirical function for the CDF of voxel intensities
    # If exclude_zero = True, then ignore zeros (as they are cropped pixels)
    vals = []
    for i in range(dmap.get_number_of_voxels()):
        val = dmap.get_value(i)
        if not exclude_zero or val != 0:
            vals.append(val)
    if len(vals) == 0:
        raise ValueError("The input density map contains no valid "
                         "voxel values.")
    f = ECDF(vals)  # Return the Empirical CDF of an array as a step function
    return f


def compute_inverse_cdf(dmap, rn=(0, 1), exclude_zero=True):
    # Returns a function that maps a CDF to the intensity value
    # f = compute_cdf(exclude_zero=exclude_zero)
    f = compute_cdf(dmap, exclude_zero=exclude_zero)
    x = np.arange(rn[0], rn[1], 0.0001)
    y = f(x)
    fi = StepFunction(y, x)
    return fi


def compare_stats_hist(ref, exp):
    # Function to debug
    expcdf = compute_cdf(exp)
    refinvcdf = compute_inverse_cdf(ref)
    refcdf = compute_cdf(ref)
    emax = exp.get_max_value()
    emin = exp.get_min_value()
    print("Val ecdf refval refcdf")
    for i in np.arange(emin, emax, 0.01):
        cc = expcdf(i)
        refval = refinvcdf(cc)
        print(i, cc, refval, refinvcdf(refval))

    print("Val ecdf refval refcdf")
    for i in np.arange(0, 1, 0.01):
        print(i, refcdf(i), refinvcdf(i))


def convert_map_to_ndarray(dmap):
    vals = []
    for i in range(dmap.get_number_of_voxels()):
        val = dmap.get_value(i)
        vals.append(val)
    map_array = np.array(vals)
    return map_array


def hist_match(source_map, template_map):
    """
    Adjust the pixel values of a Map such that its histogram
    matches that of a target Map

    Arguments:
    -----------
        source: IMP.map object
            Map to transform; the histogram is computed over the flattened
            array
        template: IMP.map object
            Template map; can have different dimensions to source
    """

    source = convert_map_to_ndarray(source_map)
    template = convert_map_to_ndarray(template_map)

    # get the set of unique pixel values and their corresponding indices and
    # counts
    s_values, bin_idx, s_counts = np.unique(source, return_inverse=True,
                                            return_counts=True)
    t_values, t_counts = np.unique(template, return_counts=True)

    # take the cumsum of the counts and normalize by the number of pixels to
    # get the empirical cumulative distribution functions for the source and
    # template maps (maps pixel value --> quantile)
    s_quantiles = np.cumsum(s_counts).astype(np.float64)
    s_quantiles /= s_quantiles[-1]
    t_quantiles = np.cumsum(t_counts).astype(np.float64)
    t_quantiles /= t_quantiles[-1]

    of = open("hmatch_new.dat", "w")

    for i in range(source_map.get_number_of_voxels()):
        expval = source_map.get_value(i)
        if expval != 0:
            id_ = np.where(s_values == expval)
#            print(s_quantiles[id_])
            # interpolate linearly to find the pixel values in the template map
            # that correspond most closely to the quantiles in the source map
            interp_t_value = np.interp(s_quantiles[id_], t_quantiles, t_values)
            of.write(str(i)+" "+str(expval)+" "+str(interp_t_value[0])+"\n")
            source_map.set_value(i, interp_t_value)


def histogram_match_map(ref, exp, em_map, df=None):
    '''
    Match the histogram of exp to ref
    '''
    # First, compute CDF for both ref and exp:
    expcdf = compute_cdf(exp)  # compute the CDF for experimental map
    # compute a mapping function that maps between the frequencies to
    # pixel values for the ref map
    refinvcdf = compute_inverse_cdf(ref)
    refcdf = compute_cdf(ref)  # compute the CDF for the ref map

    if df is not None:
        of = open(os.path.join(os.path.dirname(em_map), "hmatch.dat"), "w")
    for i in range(exp.get_number_of_voxels()):
        expval = exp.get_value(i)
        # get how frequent that pixel value (expval) occurs in
        # the experimental map
        expc = expcdf(expval)
        # get what pixel value occurs with same frequency in the ref map
        refval = refinvcdf(expc)
        if expval != 0:
            if df is not None:
                of.write(str(i)+" "+str(expval)+" "+str(expc)+" "
                         + str(refval)+" "+str(refcdf(refval))+"\n")
            exp.set_value(i, refval)


def get_percentile_of_value(dmap, value, thresh=-100000):

    # For the given value, compute its percentile over the values
    # in the density map

    dmap_vals = []
    for i in range(dmap.get_number_of_voxels()):
        val = dmap.get_value(i)
        if val > thresh:
            dmap_vals.append(val)
    pct = (np.array(dmap_vals) < value).sum() / len(dmap_vals)
    return pct


parser = argparse.ArgumentParser(
    description='Given an EMDB id, normalize the density map using '
                'histogram matching')
parser.add_argument('emdb', metavar='emdb', type=str,
                    help='Provide an EMDB id')
parser.add_argument('--thresh', metavar='thresh', type=str,
                    help='map threshold value')
parser.add_argument(
    "--database_home", dest="database_home", type=str,
    help="Directory containing all data files used in the protocol",
    default=".")

args = parser.parse_args()
curr_dir = os.getcwd()
emdb = args.emdb

# Filepaths
em_map = os.path.join(args.database_home, str(emdb), "0system", "emdb.map")
norm_em_map = os.path.join(args.database_home, str(emdb), "0system",
                           "emdb_normalized_new.map")
pdb_file = os.path.join(args.database_home, str(emdb), "0system",
                        "native.pdb")
data_file = os.path.join(args.database_home, str(emdb), "0system",
                         "normalization_new.dat")

# If there's no EM map, then there's nothing to do!
if not os.path.exists(em_map):
    print("NO em map; check path", em_map)
    exit()

df = open(data_file, "w")
df.write("voxel exp_val exp_cdf ref_val ref_cdf")

dmap = IMP.em.read_map(em_map, IMP.em.MRCReaderWriter())
ref_dmap = IMP.em.read_map(config.reference_map, IMP.em.MRCReaderWriter())

# Get the coordinates of the PDB Calpha atoms
m = IMP.Model()
rh = IMP.atom.read_pdb(pdb_file, m, IMP.atom.CAlphaPDBSelector())
particles = IMP.atom.Selection(rh).get_selected_particles()

# Crop map coordinates greater than 14 angstroms from all CA atoms
dmap_c = dmap.get_cropped(particles, 14.0)

try:
    xml = config.get_xml(emdb)
    thresh = tools.get_threshold_from_xml(xml)
except:  # noqa: E722
    mass = IMP.atom.get_mass_from_number_of_residues(len(particles))
    thresh = float(args.thresh)
    print("thresh value {} for the emdb {}".format(str(thresh), emdb))

# Then, threshold the map at the user-defined value
dmap_t = IMP.em.get_threshold_map(dmap_c, thresh)

# Match histogram of thresholded map to reference
histogram_match_map(ref_dmap, dmap_t, em_map, df=df)

# Write out normalized map
IMP.em.write_map(dmap_t, norm_em_map, IMP.em.MRCReaderWriter())
