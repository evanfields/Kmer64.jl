module Kmer64CLI

using Kmer64

function show_usage()
    println(stderr, "Usage: Kmer64CLI reads1.fastq reads2.fastq query.fasta [kmer_length=40] [check_rc=false]")
    exit(1)
end

function julia_main()::Cint
    try
        # Need at least 3 arguments
        if length(ARGS) < 3
            show_usage()
        end
        
        # Parse required arguments
        reads1_path = ARGS[1]
        reads2_path = ARGS[2]
        query_path = ARGS[3]
        
        # Parse optional arguments
        k = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 40
        check_rc = length(ARGS) >= 5 ? parse(Bool, ARGS[5]) : false
        
        # Call the filter function
        hits = Kmer64.filter_paired_reads(reads1_path, reads2_path, query_path; 
                                        k=k, check_rc=check_rc)
        
        println("Found $(length(hits)) read pairs containing query kmers")
        
        return 0
    catch e
        println(stderr, "error: $e")
        return 1
    end
end

end # module