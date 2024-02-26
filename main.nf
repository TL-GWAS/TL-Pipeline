#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// Misc Parameters
params.VERBOSITY = 0
params.RNG = 123
params.TRAITS_DATASET = "You need to provide a Traits dataset."
params.UKB_ENCODING_FILE = "NO_UKB_ENCODING_FILE"
params.CALL_THRESHOLD = 0.9
params.POSITIVITY_CONSTRAINT = 0.01
params.MAF_THRESHOLD = 0.01
params.COHORT = "UKBB"
params.UKB_CONFIG = "${projectDir}/assets/ukbconfig.yaml"
params.UKB_WITHDRAWAL_LIST = "${projectDir}/assets/NO_WITHDRAWAL_LIST"
params.OUTDIR = "${launchDir}/results"

// Confounding adjustment by PCA
params.NB_PCS = 6
params.QC_FILE = "${projectDir}/assets/NO_QC_FILE"
params.LD_BLOCKS = "${projectDir}/assets/NO_LD_BLOCKS"
params.FLASHPCA_EXCLUSION_REGIONS = "${projectDir}/assets/exclusion_regions_hg19.txt"

// Estimands Generation 
params.BATCH_SIZE = 400

// CUSTOM
params.STUDY_DESIGN = "CUSTOM"
params.ESTIMANDS_FILE = "NO_ESTIMANDS_FILE"

// FROM_ACTORS
params.EXTRA_CONFOUNDERS = "${projectDir}/assets/NO_EXTRA_CONFOUNDER"
params.EXTRA_COVARIATES = "${projectDir}/assets/NO_EXTRA_COVARIATE"
params.ENVIRONMENTALS = "${projectDir}/assets/NO_EXTRA_TREATMENT"
params.ORDERS = "1,2"

// SVP Parameters
params.SVP = false
params.GRM_NSPLITS = 100
params.NB_SVP_ESTIMATORS = 100
params.MAX_SVP_THRESHOLD = 0.9
params.ESTIMATOR_KEY = "TMLE"

// TMLE Parameters
params.KEEP_IC = params.SVP == true ? true : false
params.PVAL_THRESHOLD = 0.05
params.TL_SAVE_EVERY = 100
params.ESTIMATOR_FILE = "glmnet"

// Outputs Parameters
params.ARROW_OUTPUT = "dataset.arrow"
params.JSON_OUTPUT = "NO_JSON_OUTPUT"
params.HDF5_OUTPUT = "results.hdf5"

// Negative Control Parameters
params.PERMUTATION_HDF5_OUTPUT = "permutation_${params.HDF5_OUTPUT}"
params.PERMUTATION_JSON_OUTPUT = params.JSON_OUTPUT == "NO_JSON_OUTPUT" ? params.JSON_OUTPUT : "permutation_${params.JSON_OUTPUT}"
params.MAX_PERMUTATION_TESTS = ""
params.MAX_PERMUTATION_ATTEMPTS = 1
params.PERMUTATION_ORDERS = "1"
params.MAF_MATCHING_RELTOL = 0.05
params.VARIANTS_TO_RANDOMIZE = "NO_VARIANT_TO_RANDOMIZE"
params.N_RANDOM_VARIANTS = 10

include { TARGENE } from './workflows/targene.nf'
include { PERMUTATION_TEST; RANDOMIZATION_TEST } from './workflows/negative_control.nf'
include { PCA } from './workflows/pca.nf'
include { MAKE_DATASET } from './workflows/dataset.nf'
// Loco functions/processes
include { LocoMergeBEDS: RunPCALoco } from './modules/loco.nf'

log.info """\
         ${workflow.manifest.name} v${workflow.manifest.version}
         ==========================
         Cohort Type  : ${params.COHORT}
         Study Design : ${params.STUDY_DESIGN}
         --
         run as       : ${workflow.commandLine}
         started at   : ${workflow.start}
         config files : ${workflow.configFiles}
         container    : ${workflow.containerEngine}
         ==========================
         """
         .stripIndent()

workflow  {
    TARGENE()
}
