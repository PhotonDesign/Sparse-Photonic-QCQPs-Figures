Files in this repository and their function:
- inverse_designMMA: this completes one optimization, given the design-region parameters at the top of the code. It relies on the MMA code in the folder "GCMMA-MMA-code-1.5"
- createHundredResults: this runs inverse_designMMA 100 times (or as many as desired) and stores the relevant data
- fdfd_1d: this runs the FDFD simulation. "Helper" codes include chi_from_alpha.m, set_fdfd_grid.m, get_M_ABC.m, and sim_dir_adj.m (simulates the direct and adjoint fields)
- plotHundredResults: this code can be used to recreate the figures in the main text. Simply uncomment the desired plot

Sub-folders:
- MMAFiftyResultsNearZeroStart: 50 runs of MMA with nearly zero material at the start (starting "from scratch")
- MMAFiftyResultsRandStart: 50 runs of MMA with random layer configurations
- MMAHundredResults: compilation of the two different data sets of 50 runs each
- GCMMA-MMA-code-1.5: MMA code downloaded from https://www.smoptit.se/
- Layer-Parameterization-Old: previous inverse design approach using parameterized layers