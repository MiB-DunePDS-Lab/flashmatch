# flashmatch
A maximum likelihood-approach to match charge and light signals in DUNE-FD

## "fm" macros
You should run in sequence the following macros:
```
fm_calibrator.cpp
fm_distributions.cpp
fm_parametrizer.cpp
fm_ana.cpp
```

All the parameters are defined in **config.json**.
Here is a brief description of the parameters:
`max_nfiles`: the maximum number of files to be processed. Set it to a large number to process all the files in the input directory.
`n_opdet`: number of optical detector (it should be removed).
`input_dir`: the directory where the input files are located.

**fm_caibrator:** This script takes SolarNuAna outputs to study the relationship between the charge
of a signal TPC cluster, the neutrino energy and the drift time. It is fundamental to learn how to
estimate the neutrino energy to build the likelihood function for the flashmatch algorithm.
The most probable values of the *Q/E vs Drift Time* plot are fitted with an exponential function, while
the *Q~corr~ vs E~true~* plot is fitted with a linear function (even though it shows more complex structure).


**fm_distributions:** The script produces the *Reconstruction Probability vs Expcted #PE* plot, which
describes the probability of registering a hit and clusterise that hit in a Flash given the expected number of photoelectrons.
In addition, it constructs alsto the *Expected #PE vs Reco #PE* plot, which describes the relationship between the expected and reconstructed number of photoelectrons in a flash.
For convenience, it creates a TTree with the Cluster-Flash infos to be used in the next step of the analysis.

**fm_parametrizer:** it fits the *Reconstruction Probability vs Expcted #PE* plot with a sigmoid function, and the *Expected #PE vs Reco #PE* plot with a sigmoid-like
function, and the projection of the *Expected #PE vs Reco #PE* plot with a lognormal fognormal function.

**fm_ana:** This macro complutes the loglikelihood on (true)Cluster-(true)Flash pairs and compare it with a random
matching given by another CLuster with the same Flash.
