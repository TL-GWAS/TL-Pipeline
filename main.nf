#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

params.SNPS_EXCLUSION_LIST = "NO_FILE"
params.QUERIES_MODE = "given"
params.CALL_THRESHOLD = 0.9
params.MINOR_CAT_FREQUENCY = 0.001
params.SAVE_FULL = false
params.PHENOTYPES_BATCH_SIZE = 1
params.GRM_NSPLITS = 100
params.MAF_THRESHOLD = 0.01
params.NB_PCS = 6
params.NB_VAR_ESTIMATORS = 0
params.MAX_TAU = 0.8
params.PVAL_SIEVE = 0.05
params.OUTDIR = "$launchDir/results"

include { IIDGenotypes } from './modules/genotypes.nf'
include { generatePCs; MergeExtraCovariatesAndPCs } from './modules/covariates.nf'
include { FromASBxTransActors; FromGivenQueries } from './modules/queries.nf'
include { UKBFieldsList; UKBConv; TraitsFromUKB } from './modules/ukb_traits.nf'
include { TMLE as TMLEContinuous; TMLE as TMLEBinary} from './modules/tmle.nf'
include { PhenotypesBatches as ContinuousPhenotypesBatches; PhenotypesBatches as BinaryPhenotypesBatches} from './modules/tmle.nf'
include { GRMPart; AggregateGRM } from './modules/grm.nf'
include { SieveVarianceEstimation } from './modules/sieve_variance.nf'
include { Summary } from './modules/summary.nf'


def NbPhenotypes() {
    if (params.PHENOTYPES_LIST != "NO_FILE") {
        reader = file(params.PHENOTYPES_LIST).newReader()
        int lines = 0
        while (reader.readLine() != null) { 
            lines++
        }
        return lines
    }
    else {
        binReader = file(params.BINARY_PHENOTYPES).newReader()
        nbBin = binReader.readLine().split(" ").size()

        contReader = file(params.CONTINUOUS_PHENOTYPES).newReader()
        nbCont = contReader.readLine().split(" ").size()
        // Remove twice the FID and IID columns
        return nbBin + nbCont - 4
    }
}

NB_PHENOTYPES = NbPhenotypes()

workflow extractTraits{
    traits_config = Channel.value(file("$params.TRAITS_CONFIG"))
    encrypted_dataset = Channel.value(file("$params.ENCRYPTED_DATASET"))
    encoding_file = Channel.value(file("$params.ENCODING_FILE"))
    withdrawal_list = Channel.value(file("$params.WITHDRAWAL_LIST"))

    UKBFieldsList(traits_config)
    UKBConv(UKBFieldsList.out, encrypted_dataset, encoding_file)
    TraitsFromUKB(UKBConv.out, traits_config, withdrawal_list)

    emit:
        sample_ids = TraitsFromUKB.sample_ids
        binary_phenotypes = TraitsFromUKB.binary_phenotypes
        continuous_phenotypes = TraitsFromUKB.continuous_phenotypes
        confounders = TraitsFromUKB.confounders
        covariates = TraitsFromUKB.covariates
}

workflow generateIIDGenotypes {
    take:
        sample_ids

    main:
        qc_file = Channel.value(file("$params.QC_FILE"))
        flashpca_excl_reg = Channel.value(file("$params.FLASHPCA_EXCLUSION_REGIONS"))
        ld_blocks = Channel.value(file("$params.LD_BLOCKS"))
        bed_files_ch = Channel.fromFilePairs("$params.UKBB_BED_FILES", size: 3, checkIfExists: true){ file -> file.baseName }

        IIDGenotypes(flashpca_excl_reg, ld_blocks, bed_files_ch, qc_file, sample_ids)

    emit:
        IIDGenotypes.out
}

workflow generateQueriesAndGenotypes{
    take:
        sample_ids

    main:
        bgen_files_ch = Channel.fromPath("$params.UKBB_BGEN_FILES", checkIfExists: true)
        excluded_snps = Channel.fromPath(file("$params.SNPS_EXCLUSION_LIST"))
        if (params.QUERIES_MODE == "ASBxTransActors") {
            asb_snp_ch = Channel.fromPath("$params.ASB_FILES", checkIfExists: true)
            trans_actors = Channel.fromPath("$params.TRANS_ACTORS_FILE", checkIfExists: true)
            outputs = FromASBxTransActors(bgen_files_ch.collect(),
                                                asb_snp_ch.collect(), 
                                                trans_actors, 
                                                excluded_snps,
                                                sample_ids)
        }
        else if (params.QUERIES_MODE == "given"){
            query_files = Channel.fromPath("$params.QUERY_FILES", checkIfExists: true).collect()
            outputs = FromGivenQueries(bgen_files_ch.collect(), query_files, excluded_snps, sample_ids)
        }

    emit:
        genotypes = outputs.genotypes
        queries = outputs.queries

}

workflow generateEstimates {
    take:
        genotypes_file
        queries_files
        continuous_phenotypes_file
        binary_phenotypes_file
        covariates_file

    main:
        estimator_file = Channel.value(file("$params.ESTIMATORFILE", checkIfExists: true))

        // compute TMLE estimates for continuous targets
        ContinuousPhenotypesBatches(continuous_phenotypes_file)
        queries_to_continuous_phenotype_batches = queries_files.combine(ContinuousPhenotypesBatches.out.flatten())
        TMLEContinuous(genotypes_file, continuous_phenotypes_file, covariates_file, estimator_file, queries_to_continuous_phenotype_batches, "Real")
        
        // compute TMLE estimates for binary targets
        BinaryPhenotypesBatches(binary_phenotypes_file)
        queries_to_binary_phenotype_batches = queries_files.combine(BinaryPhenotypesBatches.out.flatten())
        TMLEBinary(genotypes_file, binary_phenotypes_file, covariates_file, estimator_file, queries_to_binary_phenotype_batches, "Bool")

        hdf5_files = TMLEContinuous.out.flatten()
                        .concat(TMLEBinary.out.flatten())
                        .map { it -> [it.getName().split("_batch")[0], it]}
                        .groupTuple()

    emit:
        hdf5_files
}


workflow generateSieveEstimates {
    take:
        snps_tmle_files
        iid_genotypes
    
    main:
        if (params.NB_VAR_ESTIMATORS != 0){
            // Build the GRM
            grm_parts = Channel.from( 1..params.GRM_NSPLITS )
            GRMPart(iid_genotypes.collect(), params.GRM_NSPLITS, grm_parts)
            AggregateGRM(GRMPart.out.collect())
            // Sieve estimation
            sieve_estimates = SieveVarianceEstimation(snps_tmle_files, AggregateGRM.out.grm_ids, AggregateGRM.out.grm_matrix)
        }
        else {
            sieve_estimates = snps_tmle_files.map(it -> [it[0], "NO_FILE"])
        }
    emit:
        sieve_estimates
}

workflow generateSummaries {
    take:
        tmle_files
        sieve_files
    
    main:
        // joining on the prefix which corresponds to a tuple of SNPS
        Summary(tmle_files.join(sieve_files))
}

workflow {
    // Extract traits
    extractTraits()

    // Generate queries
    generateQueriesAndGenotypes(extractTraits.out.sample_ids)

    // Generate IID Genotypes
    generateIIDGenotypes(extractTraits.out.sample_ids)


    // generate estimates
    generateEstimates(
        generateQueriesAndGenotypes.out.genotypes.first(),
        generateQueriesAndGenotypes.out.queries.flatten(),
        extractTraits.out.continuous_phenotypes,
        extractTraits.out.binary_phenotypes,
        extractTraits.out.confounders
        extractTraits.out.covariates
    )

    // generate sieve estimates
    generateSieveEstimates(generateEstimates.out, generateIIDGenotypes.out)

    // generate Summaries
    generateSummaries(generateEstimates.out, generateSieveEstimates.out)
}