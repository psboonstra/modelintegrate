This repository accompanies "A review of some existing and novel methods for
incorporating historical models to improve estimation of coefficients" by Philip
S. Boonstra and Pedro Orozco del Pino (2023). It allows for reproducing the
simulation study presented in that paper. 

There are five main parts to this

1. `run_sims.Rmd` / `run_sims.html` / `run_sims.R` provide the top-level
commands for conducting the simulation study. Read through the documented code
on `run_sims.html` (which is the knitted version of `run_sims.Rmd`), and when
you are ready to run it your self, use `run_sims.R` (which is the result of
`knitr::purl(run_sims.Rmd)`). Each instance of the script can generate and
analyze an arbitrary number of simulated datasets for one of the 96 data
generating mechanisms presented in the manuscript. However, for some of the data
generating mechanisms it takes many minutes or even hours to analyze a single
simulated dataset. Thus, although you can run this code locally on your own
computer, if you wish to conduct the full simulation study from the manuscript,
you will need to use an high-performance compute cluster and distribute
multiple instances of this script in parallel. The file `slurm_template.txt`
gives the skeleton of the SLURM batch script we used; you will need to fill in
the particulars before you can use it.

2. Having successfully run `run_sims.R` and saved all of the outputs into a
folder called `out` (which will be done automatically above), read through and
run `process_main_results.Rmd` / `process_main_results.html` /
`process_main_results.R` to turn the raw results into the tables and figures in
the manuscript. You can do this step on your local computer, but you obviously
first need to download `out`.

3. The folder `methods` contains R scripts that implement the statistical methods
discussed in the manuscript. It does not need to be directly called by the user. 

4. The folder `sim_functions` contains R scripts pertaining to running the
simulation study.  It does not need to be directly called by the user. 

5. The folder `aux_functions` contains other relevant R scripts with helper
functions.  It does not need to be directly called by the user. 
