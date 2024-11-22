module Kmer64CLI

using Kmer64
using ArgParse

"""
    parse_commandline()

Parse command-line arguments for using Kmer64 as a kmer filter, returning the parsed arguments as a Dict.
"""
function parse_commandline()
    
    s = ArgParseSettings(
        description = "Filter paired-end FASTQ files to find reads containing specific k-mers.",
        version = "0.0.1",
        add_version = true
    )

    @add_arg_table! s begin
        "--reads1"
            help = "Path to first FASTQ file of paired-end reads"
            required = true
        "--reads2"
            help = "Path to second FASTQ file of paired-end reads"
            required = true
        "--query", "-q"
            help = "Path to FASTA file containing query sequence"
            required = true
        "--out1", "-o"
            help = "Output path for filtered first reads"
            required = true
        "--out2", "-O"
            help = "Output path for filtered second reads"
            required = true
        "--kmer-length", "-k"
            help = "Length of k-mers to use for matching"
            arg_type = Int
            default = 40
        "--check-rc"
            help = "Also check reads for reverse complement of query kmers"
            action = :store_true
        "--single-thread"
            help = "Force single-threaded execution even if multiple threads are available"
            action = :store_true
        "--force", "-f"
            help = "Overwrite output files if they already exist"
            action = :store_true
    end

    return parse_args(s)
end

"""
    validate_args(args::Dict)

Validate the parsed command-line arguments. May throw an error.
"""
function validate_args(args::Dict)
    # Check that input files exist
    for (arg, name) in [
        ("reads1", "first reads file"),
        ("reads2", "second reads file"),
        ("query", "query sequence file")
    ]
        if !isfile(args[arg])
            error("Cannot find $name at '$(args[arg])'")
        end
    end

    # Check that k-mer length is reasonable
    if args["kmer-length"] < 1 || args["kmer-length"] > 64
        error("k-mer length must be between 1 and 64")
    end

    # Check output paths
    for outfile in [args["out1"], args["out2"]]
        # Check that output directory exists
        outdir = dirname(abspath(outfile))
        if !isdir(outdir)
            error("Output directory '$outdir' does not exist")
        end
        
        # Check if output file exists (unless --force is set)
        if !args["force"] && isfile(outfile)
            error("Output file '$outfile' already exists. Use --force to overwrite.")
        end
    end
end

"""
    julia_main()::Cint

Main entry point for the CLI application. Returns 0 on success, 1 on error.
"""
function julia_main()::Cint
    try
        args = parse_commandline()
        validate_args(args)
        
        # Determine whether to use threaded version
        using_threads = !args["single-thread"] && Threads.nthreads() > 1
        if using_threads
            println("Using $(Threads.nthreads()) threads")
            filter_fun = Kmer64.filter_paired_reads_threaded
        else
            println("Running in single-threaded mode")
            filter_fun = Kmer64.filter_paired_reads
        end

        # Do the filtering
        filter_time = @elapsed begin
            n_hits = filter_fun(
                args["reads1"], args["reads2"],
                args["query"],
                args["out1"], args["out2"];
                k = args["kmer-length"],
                check_rc = args["check-rc"],
            )
        end
        
        # Report results
        println("Filter completed in $(round(filter_time; digits=2)) seconds")
        println("Wrote $n_hits read pairs containing query k-mers")
        
        return 0
    catch e
        if e isa InterruptException
            println("\nInterrupted by user")
            return 1
        end
        println(stderr, "ERROR: $(sprint(showerror, e, catch_backtrace()))")
        return 1
    end
end

end # module
