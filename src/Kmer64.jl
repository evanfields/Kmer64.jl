module Kmer64

using BioSymbols: DNA, DNA_A, DNA_G, DNA_C, DNA_T, NucleicAcid
using BioSequences: LongDNA, DNASeq, BioSequences, isambiguous, hasambiguity

using FASTX
import Base: hash

include("kmer_filters.jl")
include("parallelize.jl")

struct Kmer
    length::UInt8
    data::UInt128
end
Base.hash(k::Kmer, h::UInt) = hash(k.data, hash(k.length, h))
# == already correctly falls back to ===

function Kmer(seq::DNASeq)
    length(seq) <= 64 || error("Can't build Kmer with more than 64 bases.")
    x = zero(UInt128)
    for nuc in seq
        x = (x << 2) | _nuc_to_bits(nuc)
    end
    return Kmer(length(seq), x)
end

"""Compute the kmer formed by taking the last k-1 bases from `kmer` and
concatenating `next_base`."""
function next_kmer(kmer::Kmer, next_base::DNA)
    empty_bits = 0x80 - (kmer.length << 0x1)
    return Kmer(
        kmer.length,
        ((kmer.data << (empty_bits + 0x2)) >> empty_bits) | _nuc_to_bits(next_base)
    )
end
function _next_data(kdata::UInt128, next_base::DNA, relevant_bits_mask)
    return ((kdata << 0x2) & relevant_bits_mask) | _nuc_to_bits(next_base)
end


const n2ba = Vector{UInt8}(undef, 8)
n2ba[reinterpret(UInt8, DNA_A)] = 0b00
n2ba[reinterpret(UInt8, DNA_T)] = 0b11
n2ba[reinterpret(UInt8, DNA_G)] = 0b01
n2ba[reinterpret(UInt8, DNA_C)] = 0b10

"""Convert a non-ambiguous DNA (G/A/T/C) to UInt8 representation.
!!! Unsafe! Will produce wrong answers or crashes if you call with any other DNA symbol."""
@inline function _nuc_to_bits(nuc::DNA)::UInt8
    return @inbounds n2ba[reinterpret(UInt8, nuc)]
end

@inline function int_to_nuc(i)
    i == 0 && return DNA_A
    i == 3 && return DNA_T
    i == 1 && return DNA_G
    i == 2 && return DNA_C
    error("may not pass int outside [0,3]")
end

function to_dna(kmer::Kmer)
    letters = DNA[]
    x = kmer.data
    for _ in 1:kmer.length
        push!(letters, int_to_nuc(x & 0b11))
        x = x >> 2
    end
    return reverse(letters)
end

function reverse_2bit_chunks(x::UInt128)
    # Swap bits within each 2-bit chunk
    x = ((x & 0x55555555555555555555555555555555) << 1) | 
        ((x & 0xaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa) >> 1)
    # Then reverse, returning each chunk to its original orientation
    return bitreverse(x)
end

function reverse_complement(kmer::Kmer)
    # slide significant bits to the left side so they'll be on the right after reversing
    n = kmer.length
    sig_bits = 2 * n
    slid = kmer.data << (128 - sig_bits)
    reversed = reverse_2bit_chunks(slid)
    # we've now reversed the sequence, take the complement
    negated_mask = typemax(UInt128) << sig_bits
    return Kmer(kmer.length, ~(reversed | negated_mask))
end

end # module Kmer64
