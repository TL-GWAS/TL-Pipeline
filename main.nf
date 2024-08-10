#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// BGEN files
params.BGEN_FILES = "${projectDir}/assests/NO_BGEN_FILE"
params.BQTLS = "${projectDir}/assests/NO_BQTLS"
params.TRANS_ACTORS = "${projectDir}/assests/NO_TRANS_ACTORS"

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
params.BATCH_SIZE = 200

// CUSTOM
params.STUDY_DESIGN = "CUSTOM"
params.ESTIMANDS_FILE = "NO_ESTIMANDS_FILE"

// SVP Parameters
params.SVP = false
params.GRM_NSPLITS = 100
params.NB_SVP_ESTIMATORS = 100
params.MAX_SVP_THRESHOLD = 0.9
params.ESTIMATOR_KEY = "TMLE"

// TMLE Parameters
params.KEEP_IC = params.SVP == true ? true : false
params.PVAL_THRESHOLD = 0.05
params.TL_SAVE_EVERY = params.BATCH_SIZE
params.ESTIMATOR_FILE = "glmnet"

// Outputs Parameters
params.ARROW_OUTPUT = "dataset.arrow"
params.JSON_OUTPUT = "NO_JSON_OUTPUT"
params.HDF5_OUTPUT = "results.hdf5"

// Simulations
params.TRAIN_RATIO = 6
params.SAMPLE_GA_HITS = true
params.GA_MAX_VARIANTS = 50
params.GA_DISTANCE_THRESHOLD = 1000000
params.GA_PVAL_THRESHOLD = 1e-5
params.GA_MAF_THRESHOLD = 0.01
params.N_REPEATS = 1
params.RNGS = [1]
params.SAMPLE_SIZES = [100]
params.GA_TRAIT_TABLE = "${projectDir}/assets/Traits_Table_GeneATLAS.csv"
params.MIN_FACTOR_LEVEL_OCCURENCES = 10
params.MAX_SAMPLING_ATTEMPTS = 1000
params.NSAMPLES_FOR_TRUTH = 1000000

include { TARGENE } from './workflows/targene.nf'
include { PCA } from './workflows/pca.nf'
include { MAKE_DATASET } from './workflows/dataset.nf'
// Loco functions/processes
include { LOCOGWAS } from './workflows/loco.nf'

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
    if (params.STUDY_DESIGN == "GWAS")
        LOCOGWAS()
    else
        TARGENE()
}
