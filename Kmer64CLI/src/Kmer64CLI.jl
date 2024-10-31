module Kmer64CLI

using Kmer64

function show_usage()
    println(stderr, "Usage: Kmer64CLI reads1.fastq reads2.fastq query.fasta out1_path out2_path [kmer_length=40 [check_rc=false]]")
    exit(1)
end

function julia_main()::Cint
    try
        # Need at least 5 arguments
        if length(ARGS) < 5
            show_usage()
        end
        
        # Parse required arguments
        reads1_path = ARGS[1]
        reads2_path = ARGS[2]
        query_path = ARGS[3]
        out1_path = ARGS[4]
        out2_path = ARGS[5]
        
        # Parse optional arguments
        k = length(ARGS) >= 6 ? parse(Int, ARGS[6]) : 40
        check_rc = length(ARGS) >= 7 ? parse(Bool, ARGS[7]) : false
        
        # Call the filter function
        hits = Kmer64.filter_paired_reads(
            reads1_path, reads2_path,
            query_path,
            out1_path, out2_path;
            k=k, check_rc=check_rc
        )
        
        println("Wrote $(length(hits)) read pairs containing query kmers")
        
        return 0
    catch e
        println(stderr, "error: $e")
        return 1
    end
end

end # module