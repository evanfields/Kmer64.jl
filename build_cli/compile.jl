using Pkg

# Activate a temporary environment for building
cd(@__DIR__)
Pkg.activate(".")

# Add PackageCompiler and your CLI package
Pkg.add("PackageCompiler")
Pkg.develop(path="../Kmer64CLI")  # Adjust path as needed

using PackageCompiler

# Create the app
create_app(
    "../Kmer64CLI", "Kmer64CLI_app";
    incremental=false,
    filter_stdlibs=true,
    force=true,
    precompile_statements_file = joinpath(@__DIR__, "precompile_statements.jl")
)

println("Build complete! Executable is at Kmer64CLI_app/bin/Kmer64CLI")