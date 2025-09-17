# Supplementary code repository for: YAP1 is a key regulator of EWS::FLI1-dependent malignant transformation upon IGF-1 mediated reprogramming of bone mesenchymal stem cells

Rahil Noorizadeh, Barbara Sax, Tahereh Javaheri, Branka Radic-Sarika, Valerie Fock, Veveeyan Suresh, Maximilian Kauer, Aleksandr Bykov, Danijela Kurija, Michaela Schlederer, Lukas Kenner, Gerhard Weber, Wolfgang Mikulits, Florian Halbritter, Richard Moriggl, Heinrich Kovar

St. Anna Children's Cancer Research Institute (CCRI), Vienna, Austria

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.14780468.svg)](https://doi.org/10.5281/zenodo.14780468)

## Abstract

Ewing sarcoma (EwS) is an aggressive cancer of adolescents in need of effective treatment. Insulin-like growth factor (IGF) 1 is an autocrine growth factor for EwS, but only 10% of patients responded to IGF-1 receptor blockade. Although presumed to originate from mesenchymal progenitors during bone development, targeting of the EwS driver oncogene EWS::FLI1  to the mesenchymal lineage in a mouse model does not result in tumor formation but to skeletal malformations and perinatal death. We report that transient exposure to IGF-1 concentrations mimicking serum levels during puberty reprogrammes limb-derived mesenchymal cells of EWS::FLI1-mutant mice to stable transformation and tumorigenicity. We identify a modular mechanism of IGF-1-driven tumor promotion in the early steps of EwS pathogenesis, in which Yap1 plays a central role. Pharmacologic Yap1/Tead inhibition reverses the transformed phenotype of EWS::FLI1 expressing cells. Our data provide a rationale for combined IGF-1R and YAP/TEAD inhibition in the treatment of EwS patients.

## Repository structure

* `project.Dockerfile` defines the environment used to carry out all experiments
* `config.yaml` is used to set paths 
* `R/` holds R function definitions and misc utility scripts
* `Rmd/` holds R markdown documents for the individual steps of the project
* `bash/` holds shell scripts to build and run the docker image, and to parse the config file
* `metadata/` holds custom geneset definitions

## Reproducing the results

The file `R/knit.R` calls all `Rmd/*.rmd` files in order to reproduce the analysis.

Paths in the `config.yaml` file starting with "/path/to/" will have to be set.

## Links

**Paper:** [10.1016/j.celrep.2025.115381](https://doi.org/10.1016/j.celrep.2025.115381)

**Data (GEO):** [GSE269007](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE269007)
