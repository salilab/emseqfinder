\brief Assignment of sequence to backbone fragments traced in a cryo-EM map 

EMSequenceFinder is a method for assigning amino acid residue sequence to
backbone fragments traced in an input cryo-electron microscopy (cryo-EM) map.
EMSequenceFinder relies on a Bayesian scoring function for ranking 20 standard
amino acid residue types at a given backbone position, based on the fit to a
density map, map resolution, and secondary structure propensity. The fit to a
density is quantified by a convolutional neural network that was trained on
5.56 million amino acid residue densities extracted from cryo-EM maps at
3–10 Å resolution and corresponding atomic structure models deposited in
the Electron Microscopy Data Bank (EMDB). For more information, see
[Mondal et al, 2025](https://pubmed.ncbi.nlm.nih.gov/40719420/).

# Dependencies:

In addition to IMP's own dependencies, EMSequenceFinder also requires

 - The [STRIDE](https://webclu.bio.wzw.tum.de/stride/install.html)
   command line tool for secondary structure prediction
 - The `mrcfile`, `scipy`, `scikit-learn`, `statsmodels`, `pandas`,
   and `tensorflow` Python packages

One way to get these dependencies is via
[conda-forge](https://conda-forge.org/). In order for TensorFlow prediction
to work correctly on GPUs with `libdevice`, you may have to run
`export XLA_FLAGS=--xla_gpu_cuda_data_dir=$CONDA_PREFIX`

# emseqfinder: command line tool to run emseqfinder protocol {#emseqfinder}

The protocol is typically run using the `emseqfinder` command line tool
in the following fashion:

 1. Add input files. Create directories `pdb_files`, `cryoem_maps` and
    `fasta_files` containing input files in `.pdb`, `.map` and `.fasta` format
    respectively. Name all three files for a given run with the same stem
    (e.g. `EMD-8637.pdb`, `EMD-8637.map`, `EMD-8637.fasta`).
 2. Run the pipeline with `emseqfinder batch`.

Upon successful execution, the following output files will be generated:
 - `*_ML_side_ML_prob.dat` contains fragment-wise sequence scores.
 - `batch_matching_results.txt` contains overall sequence matching accuracy
   per structure.

# Info

_Author(s)_: Dibyendu Mondal, Vipul Kumar

_Maintainer_: `benmwebb`

_License_: [LGPL](https://www.gnu.org/licenses/old-licenses/lgpl-2.1.html)
This library is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2 of the License, or (at your option) any later version.

_Publications_:
 - D. Mondal, V. Kumar, T. Satler, R. Ramachandran, D. Saltzberg, I. Chemmama, K.B. Pilla, I. Echeverria, B.M. Webb, M. Gupta, K.A. Verba, A. Sali. Recognizing amino acid sidechains in a medium resolution cryo-electron density map. Prot Sci 34, e70217, 2025.
 - See [main IMP papers list](@ref publications).
