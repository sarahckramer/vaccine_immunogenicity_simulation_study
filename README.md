# vaccine_immunogenicity_simulation_study
The purpose of this repository is to test a method to infer parameter values from antibody kinetic data using nonlinear mixed effects models.

All code is contained in the src directory.

Synthetic antibody kinetic data with bi- or monoexponential decay patterns can first be generated using the python scripts, "generate_synthetic_antibody_data_BIEXP.py" and "generate_synthetic_antibody_data_MONOEXP.py." All required functions for generating data are contained within "functions_python.py."

The R scripts, "fit_model_to_synth_data_BIEXP.R" and "fit_model_to_synth_data_MONOEXP.R," can then be run to fit bi- and monoexponential models, respectively, to the synthetic data, using the R package [nlme][1]. The code loops over n_participants to include varying numbers of "individuals" when fitting the models, and fits to data at the timepoints specified in select_timepoints. Both of these vectors can be edited to explore additional combinations of participant numbers and blood sampling timepoints, but only those included in the code have been rigorously tested. Note that fitting the biexponential model typically requires at least four blood sampling timepoints; the monoexponential model can handle as few as two post-vaccination timepoints. The parameter "noisy" controls whether results and diagnostic plots will be printed as the code runs.

Finally, the R script "compare_to_truth.R" generates plots comparing the estimated model parameters to the "true" values initially used to generate the data in the python code.

This repository uses the R package [renv][2] to handle package version control. To load the versions of R packages used in this repository, you can run "renv::restore()" before running any of the R scripts.

[1]: https://www.rdocumentation.org/packages/nlme/versions/3.1-152
[2]: https://rstudio.github.io/renv/articles/renv.html
