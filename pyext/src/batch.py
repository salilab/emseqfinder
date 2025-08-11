import subprocess
import shutil
import pathlib
import sys
from IMP import ArgumentParser
from IMP.emseqfinder.compute_dynamic_threshold import compute_threshold


__doc__ = "Perform all steps of the emseqfinder protocol."


def process_pdb(pdbfile):
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
         '--thresh', str(threshold)])
    if p.returncode != 0:
        print(f"[ERROR] Normalization failed. Skipping {base}.")
        return


def parse_args():
    parser = ArgumentParser(
        description="Perform all steps of the emseqfinder protocol")
    return parser.parse_args()


def main():
    args = parse_args()

    # Resolution used for database generation
    resolution = 4

    # Output file for overall results
    final_output_file = pathlib.Path("batch_matching_results.txt")

    # Add header if file doesn’t exist
    if not final_output_file.exists():
        with open(final_output_file, 'w') as fh:
            print("Result_File Total_Percentage_Matched "
                  "Total_Abs_Percentage_Matched", file=fh)

    # Loop through all PDB files
    for pdbfile in pathlib.Path("pdb_files").glob("*.pdb"):
        process_pdb(pdbfile)


if __name__ == '__main__':
    main()
