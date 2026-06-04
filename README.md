# brain-controllability-stimulation
Brain stimulation and network controllability analysis using MATLAB and Julia

## Computational Pipeline

The analysis consists of a two-stage workflow combining Julia and MATLAB implementations.

### 1. Julia simulation

The first step is to run the Julia simulation code:

This code requires the following input files:

- `StructuralMatrix30.mat`
- `IntCon30.mat`

These files contain the structural connectivity matrix and initial conditions for the 30-subject dataset.

The output of this step includes the following stimulation-derived measures for each subject:

- FC_before_stim
- FC_during_stim
- FC_post_stim

The results are saved and used as input for the MATLAB analysis pipeline.

---

### 2. MATLAB analysis

After running the Julia simulation, navigate to the folder containing the generated output files and run the MATLAB scripts.

The MATLAB code is used to:

- Load simulation outputs (FC_before_stim, FC_during_stim, FC_post_stim)
- Compute statistical relationships with controllability metrics (AC and MC)
- Generate figures and correlation analyses
- Perform subject-level and group-level comparisons

Users can select which relationships and figures to compute depending on the desired analysis (e.g., AC–SEDB, MC–FEPB, etc.).

---

### 3. Notes

- Ensure that all required `.mat` files are in the correct directory before running the Julia code.
- MATLAB scripts assume that the working directory is set to the output folder of the Julia simulation.
- Results can be reproduced by running the pipeline in the specified order.
