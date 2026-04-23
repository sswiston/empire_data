EMPIRE: The Ellipse Model for Phylogenetic Inference of Range Evolution
---

This dataset is associated with Swiston, McHugh, & Landis 2026 (DOI TBD). It contains two major elements. First is a simulation study performed using R associated with the validation of the EMPIRE model. Second is an empirical analysis of the Australian Sphenomorphinae under the EMPIRE model. Here, we provide the version of the EMPIRE R package used in this study, scripts for simulation and empirical studies, empirical data including the phylogenetic tree and present-day species ranges for the Australian Sphenomorphinae, and output from the empirical analysis with the MultiFIG model. We also provide figures and supplemental figures associated with the manuscript.

Description of the data and file structure

Overview: The **package** .zip file contains the version of the EMPIRE R package used in this study (for maintained package, see https://github.com/sswiston/empire). The **code** .zip file contains all of the scripts necessary for generating and analyzing simulations for the simulation study, running the empirical analysis of the Australian Sphenomorphinae, and plotting the results of both analyses. The **data** .zip file contains the empirical dataset for the analysis of the Australian Sphenomorphinae. The **output** .zip file contains all of the unprocessed and processed output from the empirical analysis of the Australian Sphenomorphinae (**skinks**), the processed output from the simulation study (**sim**), as well as all plots and supplemental figures from the associated manuscript (**plots**).

* `package` contains the version of the EMPIRE R package used in this study - for more information, you can download the package at https://github.com/sswiston/empire and examine the vignette (with a worked example) and documentation for objects and functions
    * `README.md` contains basic package information
    * `LICENSE` contains license information
    * `DESCRIPTION`, `NAMESPACE` provide basic package structure
    * `R` contains .R files for all package objects and functions
    * `man` contains documentation for all package objects and functions
* `code` contains R code for EMPIRE simulation study and empirical analysis of the Australian Sphenomorphinae, as well as plotting scripts.
    * `sim.R` the script for simulating standard datasets and analyzing those simulated datasets under EMPIRE
    * `sim_sensitivity.R` the script for simulating datasets with missing taxa and analyzing those simulated datasets under EMPIRE
    * `sim_noisy.R` the script for simulating noisy datasets and anlyzing those simulated datasets under EMPIRE
    * `sim_rase.R` the script for simulating standard datasets under EMPIRE and analyzing those simulated datasets using RASE
    * `skinks.R` the script for analyzing the Australian Sphenomorphinae under EMPIRE
    * `plots.Rmd` the script for processing and plotting output
    * `plots.html` the knitted document produced by `plots.Rmd`
* `data` contains data for the empirical analysis of the Australian Sphenomorphinae
    * `australia` contains data associated with the Australian continent and environment
        * `australia.cpg`, `australia.dbf`, `australia.shp`, `australia.shx` comprise the shapefile for plotting the outline of Australia
        * `Ma0_raster.csv`, `Ma20_raster.csv`, `Ma40_raster.csv` are precipitation rasters from Pohl et al. 2022 for time slices 0Ma, 20Ma, and 40Ma
        * `bom_raster.csv` is the precipitation raster from the Australian Government Bureau of Meteorology for the present (0Ma)
    * `ranges` contains range data associated with the Australian Sphenomorphinae
        * `australian_skinks.dbf`, `australian_skinks.prj`, `australian_skinks.shp`, `australian_skinks.shx` comprise the shapefile for the ranges of each extant species of Australian Sphenomorphine, obtained from Roll et. al. 2017 and Caetano et al. 2022
        * `ellipse_data.csv` contains the range ellipses for each extant species of Australian Sphenomorphinae
    * `spheno.tre` contains a cropped phylogeny from Title et al. 2024 for the Australian Sphenomorphinae
* `output` contains output of simulated and empirical analyses
    * `sim` contains output of simulation study
        * `sim_results.tsv` contains the output for standard simulated datasets:
            * sim = simulation #
            * param = parameter name
            * true = true simulating value
            * est = estimate (posterior mean)
            * hpd_low = 95% HPD lower bound
            * hpd_high = 95% HPD upper bound
            * covered = whether the 95% HPD interval covers the truth
            * ess = effective sample size
            * tree_size = number of taxa in the tree
        * `sensitivity_results.tsv` contains the output for simulated datasets with missing taxa:
            * sim = simulation #
            * param = parameter name
            * true = true simulating value
            * est = estimate (posterior mean)
            * hpd_low = 95% HPD lower bound
            * hpd_high = 95% HPD upper bound
            * ess = effective sample size
        * `noisy_results.tsv` contains the output for simulated datasets with added noise:
            * sim = simulation #
            * param = parameter name
            * true = true simulating value
            * est = estimate (posterior mean)
            * hpd_low = 95% HPD lower bound
            * hpd_high = 95% HPD upper bound
            * ess = effective sample size
        * `rase_results.tsv` contains the output for simulated datasets analyzed with RASE:
            * sim = simulation #
            * param = parameter name
            * true = true simulating value
            * empire = estimate from EMPIRE (posterior mean)
            * rase = estimate from RASE (posterior mean)
            * ess = effective sample size
            * mean_a = mean log area across simulated ancestral ellipses
            * ratio = the ratio of the absolute difference between RASE estimate and the true value divided by the true value
    * `skinks` contains output of empirical analysis of Australian Sphenomorphinae
        * `*_log.tsv` contain the MCMC traces from an EMPIRE MCMC analysis
            * `model` contains model parameters
            * `d` contains daughter configuration at each internal node
            * `m` contains cladogenetic mode at each internal node
            * `c` contains concentric circle of splitting at each internal node
            * `h` contains direction line of splitting at each internal node
            * `x` contains x coordinate of centroid location at each internal node
            * `y` contains y coordinate of centroid location at each internal node
            * `r` contains x coordinate of tilt point (controls elongation) at each internal node
            * `s` contains y coordinate of tilt point (controls elongation) at each internal node
            * `a` contains log area at each internal node
        * `node_reconstructions.tsv` contains reconstructed ellipse data for each ancestral range ellipse:
            * node = internal node number
            * r/s/a/x/y = posterior means
            * scenario_d/m/c/h = the values associated with the most common scenario (combination of d/m/c/h)
            * scenario_r/s/a/x/y = posterior means of continuous characters across iterations with the selected scenario
            * scenario_sample_size = the number of iterations with the selected scenario
            * scenario_freq = the relative frequency of the selected scenario
            * scenario_d/m/c/h_freq = the relative frequency of the values associated with the selected scenario
            * scenario_r/s/a/x/y_low = the lower bound of the 95% HPD across iterations with the selected scenario
            * scenario_r/s/a/x/y_high = the upper bound of the 95% HPD across iterations with the selected scenario
            * scenario_r/s/a/x/y_range = the range of the 95% HPD across iterations with the selected scenario
            * node_age = the age of the internal node
            * selection = the direction line, grouping symmetrical directions
            * true_x = the x component of the true (projected) splitting direction
            * true_y = the y component of the true splitting direction
            * true_angle = the angle of the true splitting direction
            * true_selection = the true splitting direction binned according to the closest direction line
        * `aridity_data.csv` contains information about the aridification experienced by each ancestral range ellipse:
            * node = internal node number
            * reconstructed = mean aridification across the reconstructed range ellipse in mm/day/MY
            * mean = mean aridification across all possible daughters from the reconstructed range ellipse
            * clade = which clade the node belongs to
        * `processed_directions.csv` contains information about the true (projected) splitting direction at each internal node for a sample of the MCMC iterations:
            * h = the direction line of splitting
            * r = the x coordinate of tilt point (controls elongation)
            * s = the y coordinate of tilt point (controls elongation)
            * x = the x coordinate of centroid location
            * y = the y coordinate of centroid location
            * selection = the direction line, grouping symmetrical directions
            * true_x = the x component of the true splitting direction
            * true_y = the y component of the true splitting direction
            * true_angle = the angle of the true splitting direction
            * true_selection = the true splitting direction binned according to the closest direction line
            * node = internal node number
    * `plots` Contains plots associated with Swiston, McHugh, & Landis 2026 (DOI TBD)

Sharing/Access information

Links to other publicly accessible locations of the data:

* https://bitbucket.org/sswiston/empire_data/src/main/

Data was derived from the following sources:
* Title, P. O., Singhal, S., Grundler, M. C., Costa, G. C., Pyron, R. A., Colston, T. J., Grundler, M. R., Prates, I., Stepanova, N., Jones, M. E. H., Cavalcanti, L. B. Q., Colli, G. R., Di-Poï, N., Donnellan, S. C., Moritz, C., Mesquita, D. O., Pianka, E. R., Smith, S. A., Vitt, L. J., & Rabosky, D. L. (2024). The macroevolutionary singularity of snakes. Science, 383(6685), 918–923. https://doi.org/10.1126/science.adh2449
* Roll et. al. 2017 The global distribution of tetrapods reveals a need for targeted reptile conservation. Nature Ecology & Evolution 1:1677-1682
* Caetano, et al. 2022. Automated assessment reveals that the extinction risk of reptiles is widely underestimated across space and phylogeny. PLoS Biology, 20(5): e3001544.
* Pohl, A., Wong Hearing, T., Franc, A., Sepulchre, P., & Scotese, C. R. (2022). Dataset of Phanerozoic continental climate and Köppen–Geiger climate classes. Data in Brief, 43(108424), 108424. https://doi.org/10.1016/j.dib.2022.108424
* http://www.bom.gov.au/climate/maps/averages/rainfall/
* Swiston, S. K., McHugh, S. W., & Landis, M. J. (2026). EMPIRE: The Ellipise Model for Phylogenetic Inference of Range Evolution. BioRxiv. DOI TBD

Code/Software

* The EMPIRE R package: https://github.com/sswiston/empire
* The dataset contains .R and .Rmd files for performing both simulation and empirical analyses under the EMPIRE model. These files were designed to be run using R version 4.4.0.

Version Changes:
NONE