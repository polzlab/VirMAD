# using CairoMakie
# using Gtk
struct Contig
    name::String
    length::Int
end
struct Interval
    block::Int
    genome::String
    contigInd::Int
    s::Int
    e::Int
    strand::Char
end
sfx(x::Vector{T}) where {T<:AbstractString} = split.(x, '\t');
getgen(x::AbstractString) = replace(x, r"_contig.+$" => "")
rng(x) = x[1]+1:x[2]
rng(x,i) = x[i,1]+1:x[i,2]
rng2(x,i) = x[i,1]+2:x[i,2]
nrow(x) = size(x,1)
flatten(x) = reduce(vcat,x)
blocklen(x) = x[3]-x[2]+1
function get_best_gen(block_desc)
    gens = [getgen(x[3]) for x in block_desc]
    lend = Dict(g => 0 for g in gens)
    for x in block_desc
        k = getgen(x[3])
        lend[k] += 1
    end
    o = sort(collect(lend), by=x -> x[2])
    return o[1][1]
end

# read block_coords file
p = "test/500/blocks_coords.txt"
lns = readlines(p) |> sfx;
breaks = [length(x)==1 && startswith("---")(x[1]) for x in lns] |> findall
lns_blocks = [[1; breaks .+ 1] [breaks .- 1; length(lns)]]


# get genome with smallest number of contigs and its IDs
refi = rng(lns_blocks,1)
ref_gen = get_best_gen(lns[refi])
ref_ids = [parse(Int,x[1]) for x in lns[refi] if startswith(ref_gen)(x[3])]
ref = lns[refi];
cl = parse(Int,ref[end][1])
contigs = Vector{Contig}(undef,cl)
for x in ref
    ci = parse(Int,x[1])
    contigs[ci] = Contig(x[3],parse(Int,x[2]))
end

# extract all intervals and build ref_ints which is 
# the intervals in the ref genome, ordered by contig
# length
intervals = Vector{Interval}();
for i in 2:size(lns_blocks,1)
    for x in lns[rng2(lns_blocks,i)]
        ci = parse(Int,x[1])
        gen = getgen(contigs[ci].name)
        co = parse.(Int, x[3:4])
        i1 = min(co...)
        i2 = max(co...)
        push!(intervals, Interval(i-1,gen,ci,i1,i2,x[2][1]))
    end
end;
sort!(intervals,by=x->(x.genome,x.contigInd,x.s));
ref_ints = filter(x -> x.genome==ref_gen, intervals);
coo = sortperm(contigs[ref_ids], by=x->-x.length)
co = ref_ids[coo]
cis = [x.contigInd for x in ref_ints];
io = reduce(vcat,[findall(cis.==i) for i in co]);
ref_ints = ref_ints[io];

spots = Vector{Vector{Interval}}(undef,length(ref_ints));
for (i,r) in enumerate(ref_ints)
     spots[i] = filter(x -> x.block==r.block, intervals)
end
