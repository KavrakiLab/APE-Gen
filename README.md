# APE-Gen: Anchored Peptide-MHC Ensemble Generator

## Description:

- APE-Gen is a tool that generates multiple clash-free conformations of a peptide bound to an MHC

- All that is required for input is the sequence of the peptide and the MHC allotype name (if supported)

- Minimal example: `python APE_Gen.py QFKDNVILL HLA-A*24:02`

- For a peptide with "n" residues, APE-Gen will use a template of another peptide with "n" residues to place the anchor residues in the correct pocket of the MHC
  - to change the template, modify `n-mer-templates.txt` for the corresponding n-mer, then add the template pdb file (peptide and MHC) to the `templates/` folder

- A receptor model is also needed to run APE-Gen
  - if the allotype of the MHC is supported (found in `receptor-class-templates.txt`), then simply specify the name as input into APE-Gen
  - if the allotype is not supported, but one has the PDB file of the receptor
    - place the PDB file inside the `templates/` folder
    - modify `receptor-class-templates.txt` with an appropriate name
    - then call APE-Gen with the chosen name
  - One way to obtain a model of the receptor is to use the provided scripts in `modeller_scripts/` for homology modelling (see below)

- PDB File preparation:
  - For PDB files to be used as input anywhere in the script, the chains must be labelled in a particular way
  - chain A is the heavy alpha chain, chain B is the light beta-immunoglobulin chain, chain C is the peptide

- Specifying receptor degrees-of-freedom
  - modify `flex_res.txt` to add/remove residues that are allowed to be flexible during the SMINA minimization step

- `dunbrack.bin` and `loco.score` are files required for the RCD step
- `align.py` is required for using PYMOL for alignment of PDB files
- `minimize.py` is needed for further minimization using OpenMM


## Output

- Each round of APE-Gen is saved within a folder with the index of the round (counting from 0)

- `full_system_confs/` contains the ensemble of conformations of peptide and MHC after energy refinement and filtering
  - each conformation is named by the index of the loop as generated by RCD
    - not every conformation makes it past the energy refinement and filtering step
- `peptide_confs.pdb` contains only the peptide conformations
- `filtered_energies.npz` is a numpy file that contains the energies of the ensemble according to the SMINA scoring function
  - it contains two arrays of the same size: `filtered_indices` which contains indices of each conformation and `filtered_energies` which contains the corresponding energies


## Helper scripts

- `get_pMHC_pdb.py`
  - usage: `python get_pMHC_pdb.py <pdb code>`
  - assumes pdb code is of a peptide-MHC structure
  - adds missing atoms/residues, removes all waters and ions, labels chains as A,B,C where chain C is the peptide

- `mutate.py`
  - Usage: `~/pymol/bin/pymol -qc mutate.py <pdb> <selection> <new_residue in 3-letter code> <name of new pdb file>`
  - example: `~/pymol/bin/pymol -qc mutate.py 0.pdb C/1/ ALA 0_mutated.pdb`

## Options help:

```
usage: APE-Gen.py [-h] [-n NUM_CORES] [-l NUM_LOOPS] [-t RCD_DIST_TOL] [-r]
                  [-d] [-p] [-a ANCHOR_TOL] [-o] [-g NUM_ROUNDS]
                  [-b {receptor_only,pep_and_recept}] [-s]
                  peptide_input receptor_class

Anchored Peptide-MHC Ensemble Generator

positional arguments:
  peptide_input         Sequence of peptide to dock or pdbfile of crystal
                        structure
  receptor_class        Class descriptor of MHC receptor. Use REDOCK along
                        with crystal input to perform redocking.

optional arguments:
  -h, --help            show this help message and exit
  -n NUM_CORES, --num_cores NUM_CORES
                        Number of cores to use for RCD and smina computations.
                        (default: 8)
  -l NUM_LOOPS, --num_loops NUM_LOOPS
                        Number of loops to generate with RCD. (Note that the
                        final number of sampled conformations may be less due
                        to steric clashes. (default: 100)
  -t RCD_DIST_TOL, --RCD_dist_tol RCD_DIST_TOL
                        RCD tolerance (in angstroms) of inner residues when
                        performing IK (default: 1.0)
  -r, --rigid_receptor  Disable sampling of receptor degrees of freedom
                        specified in flex_res.txt (default: False)
  -d, --debug           Print extra information for debugging (default: False)
  -p, --save_only_pep_confs
                        Disable saving full conformations (peptide and MHC)
                        (default: False)
  -a ANCHOR_TOL, --anchor_tol ANCHOR_TOL
                        Anchor tolerance (in angstroms) of first and last
                        backbone atoms of peptide when filtering (default:
                        2.0)
  -o, --score_with_openmm
                        Rescore full conformations with openmm (AMBER)
                        (default: False)
  -g NUM_ROUNDS, --num_rounds NUM_ROUNDS
                        Number of rounds to perform. (default: 1)
  -b {receptor_only,pep_and_recept}, --pass_type {receptor_only,pep_and_recept}
                        When using multiple rounds, pass best scoring
                        conformation across different rounds (choose either
                        'receptor_only' or 'pep_and_recept') (default:
                        receptor_only)
  -s, --min_with_smina  Minimize with SMINA instead of the default Vinardo
                        (default: False)
```


## Installation instructions:

1) install miniconda 
	- https://conda.io/miniconda.html

2) using conda, install smina, pdbfixer, numpy, mdtraj, and openMM
	- `conda install -c bioconda smina`
	- `conda install -c omnia pdbfixer`
	- `conda install -c conda-forge mdtraj`
  - `conda install -c schrodinger pymol`
  - (optional) `conda install -c omnia -c conda-forge openmm`

3) install RCD
	- http://chaconlab.org/modeling/rcd/rcd-download
	- make sure RCD is added to path so that `rcd` is a command in the terminal
  - intel mkl may be needed (`conda install -c intel mkl`) and added to library path

4) install vina
	- `wget vina.scripps.edu/download/autodock_vina_1_1_2_linux_x86.tgz`
	- `tar -xvzf autodock_vina_1_1_2_linux_x86.tgz`
	- then add bin into path so that `vina_split` is a command in the terminal


## Using Modeller scripts

- Usage: `python model_receptor.py <fasta file containing alpha chain seq> <template PDB>`
- Make sure all files are in the same directory
- Will create 2 models and choose the best model according to MODELLER's DOPE score
  - can change this parameter inside `model_receptor.py`
- Model with the best DOPE score is found in `best_model.pdb`
- Example: `python model_receptor.py P01892.fasta 3I6L.pdb`
  - models HLA-A*02:01 using 3I6L as a template
    - 3I6L contains a model of HLA-A*24:02

## Using Docker file

- Run: `docker build -t apegen .`
- Minimal example: `python /APE-Gen/APE_Gen.py LLWTLVVLL HLA-A*02:01 -o`
- Modeller still needs to be installed (along with obtaining license key) 

