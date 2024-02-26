include { SVP; AggregateGRM; GRMPart } from '../modules/svp.nf'

workflow SVPWorkflow {
    take:
        hdf5_result
        iid_genotypes
    
    main:
        grm_parts = Channel.from( 1..params.GRM_NSPLITS )
        GRMPart(iid_genotypes.collect(), params.GRM_NSPLITS, grm_parts)
        AggregateGRM(GRMPart.out.collect())
        // Sieve estimation
        SVP(
            hdf5_result.collect(), 
            AggregateGRM.out.grm_ids, 
            AggregateGRM.out.grm_matrix,
        )
}