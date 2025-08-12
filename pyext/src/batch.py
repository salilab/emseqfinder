import subprocess
import shutil
import contextlib
import pathlib
import tempfile
import gzip
import os
import sys
from . import get_data_path
from IMP import ArgumentParser
from IMP.emseqfinder.compute_dynamic_threshold import compute_threshold
from IMP.emseqfinder.calculate_seq_match_batch import calculate_seq_match


__doc__ = "Perform all steps of the emseqfinder protocol."


def process_pdb(pdbfile, resolution, database_home, reference_map,
                final_output_file):
    base = pathlib.Path(pdbfile.stem)
    mapfile = pathlib.Path('cryoem_maps') / base.with_suffix('.map')
    fastafile = pathlib.Path('fasta_files') / base.with_suffix('.fasta')

    # Check file existence
    if not mapfile.exists():
        print(f"[WARNING] Missing map file for {base!r}. Skipping...")
        return
    if not fastafile.exists():
        print(f"[WARNING] Missing FASTA file for {base!r}. Skipping...")
        return

    print("===============================================")
    print(f"[INFO] Processing {base}")

    threshold = compute_threshold(mapfile)
    print(f"[INFO] Using threshold: {threshold:.4f}")

    # Clean and recreate project directory
    shutil.rmtree(base, ignore_errors=True)
    (base / '0system').mkdir(parents=True)
    (base / '1structure_elements').mkdir()

    # STRIDE assignment
    stride_file = base / base.with_suffix('.stride')
    stride_file.unlink(missing_ok=True)
    p = subprocess.run(["stride", pdbfile, "-rA", "-f" + str(stride_file)])
    if p.returncode != 0:
        print(f"[ERROR] STRIDE failed for {base}. Skipping...")
        return

    # Create unique fragment directory
    frag_dir = pathlib.Path(f"fragments_{base}")
    shutil.rmtree(frag_dir, ignore_errors=True)
    frag_dir.mkdir(parents=True)
    p = subprocess.run(
        [sys.executable, '-m',
         'IMP.emseqfinder.fragdb.get_fraglib_from_native', pdbfile,
         stride_file, frag_dir])
    if p.returncode != 0:
        print(f"[ERROR] Fragment generation failed. Skipping {base}.")
        return

    # Copy system files
    shutil.copy(pdbfile, base / "0system" / "native.pdb")
    shutil.copy(mapfile, base / "0system" / "emdb.map")
    shutil.copy(fastafile, base / "0system")
    for frag in frag_dir.glob("*.pdb"):
        shutil.copy(frag, base / "1structure_elements")

    # Normalize map
    p = subprocess.run(
        [sys.executable, '-m',
         'IMP.emseqfinder.mldb.normalize_map_for_parts_fitting', base,
         '--thresh', str(threshold),
         '--database_home', database_home,
         '--reference_map', reference_map])
    if p.returncode != 0:
        print(f"[ERROR] Normalization failed. Skipping {base}.")
        return

    # Remove old database files if exist
    for p in (f"{base}_ML_side.dat", f"{base}_ML_side.pkl",
              f"{base}_ML_side_ML_prob.dat", f"{base}_ML_output.txt"):
        pathlib.Path(p).unlink(missing_ok=True)

    # Create database
    p = subprocess.run(
        [sys.executable, '-m',
         'IMP.emseqfinder.mldb.get_database_for_one_emdb_using_parts', base,
         '--database_home', database_home,
         f"{base}_ML_side.dat", str(resolution)])
    if p.returncode != 0:
        print(f"[ERROR] ML DB generation failed. Skipping {base}.")
        return

    # Convert to pickle
    p = subprocess.run(
        [sys.executable, '-m', 'IMP.emseqfinder.convert_MLDB_topkl',
         f"{base}_ML_side.dat", f"{base}_ML_side"])
    if p.returncode != 0:
        print(f"[ERROR] PKL conversion failed. Skipping {base}.")
        return

    # Prediction
    p = subprocess.run(
        [sys.executable, '-m', 'IMP.emseqfinder.final_ML_predict',
         f"{base}_ML_side.pkl", '10000'])
    if p.returncode != 0:
        print(f"[ERROR] Prediction failed. Skipping {base}.")
        return

    # Evaluate prediction
    p = subprocess.run(
        [sys.executable, '-m', 'IMP.emseqfinder.evaluate_output_database',
         f"{base}_ML_side_ML_prob.dat", f"{base}_ML_output.txt"])
    if p.returncode != 0:
        print(f"[ERROR] Evaluation failed. Skipping {base}.")
        return

    # Calculate and append sequence match
    calculate_seq_match([f"{base}_ML_output.txt"], final_output_file)

    # Cleanup
    shutil.rmtree(frag_dir, ignore_errors=True)

    print(f"[INFO] Finished processing {base} ✅")
    print("===============================================")
    print("")


def parse_args():
    parser = ArgumentParser(
        description="Perform all steps of the emseqfinder protocol")
    parser.add_argument(
        "--db-resolution", dest="resolution", type=float,
        help="Resolution used for database generation", default=4.0)
    parser.add_argument(
        "--database_home", dest="database_home", type=str,
        help="Directory containing all data files used in the protocol",
        default=".")
    parser.add_argument(
        "--reference_map", dest="reference_map", type=str,
        help="reference map for voxel data extraction",
        default=get_data_path('reference/ref.mrc.gz'))

    return parser.parse_args()


@contextlib.contextmanager
def get_reference_map(args):
    if args.reference_map.endswith('.gz'):
        # IMP's mrc reader cannot read compressed files, so decompress
        # to a temporary location and return that instead
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpfile = os.path.join(tmpdir, 'ref.mrc')
            with gzip.open(args.reference_map, 'rb') as fh_in:
                with open(tmpfile, 'wb') as fh_out:
                    shutil.copyfileobj(fh_in, fh_out)
            yield tmpfile
    else:
        yield args.reference_map


def main():
    args = parse_args()

    # Output file for overall results
    final_output_file = pathlib.Path("batch_matching_results.txt")

    with get_reference_map(args) as reference_map:
        for pdbfile in pathlib.Path("pdb_files").glob("*.pdb"):
            process_pdb(pdbfile, args.resolution, args.database_home,
                        reference_map, final_output_file)


if __name__ == '__main__':
    main()
