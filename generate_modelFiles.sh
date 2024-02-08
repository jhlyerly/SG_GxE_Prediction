#! /bin/bash

rm -R Intermediate_Outputs

Rscript s1_process_phenoGeno.R

Rscript s2_getEnvironmentalData.R
