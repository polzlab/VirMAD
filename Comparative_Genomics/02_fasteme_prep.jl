using StringViews
using Mmap
using LinearAlgebra
using SparseArrays
using DataFrames
make_range(I, i) = (I[i-1]+1):(I[i]-1)
function split_lines(x)
    I = findall(x .== UInt8('\n'))
    [StringView(x[make_range(I, i)]) for i in 2:length(I)]
end
function rem_ext(x::T)::T where {T<:AbstractString}
    join(split(x, '.')[1:end-1], '.')
end
function process_raw(x, g2i)
    lns = split_lines(x)
    lnss = split.(lns, '\t')
    g1vec = basename.(getindex.(lnss, 1)) |> x -> rem_ext.(x)
    g2vec = basename.(getindex.(lnss, 2)) |> x -> rem_ext.(x)
    ani = parse.(Float64, getindex.(lnss, 3))
    af1 = parse.(Float64, getindex.(lnss, 4))
    af2 = parse.(Float64, getindex.(lnss, 5))
    af = min.(af1,af2)
    edges =  DataFrame(g1 = g1vec, g2 = g2vec, ani = ani, af = af)
    edges.i = [g2i[k] for k in edges.g1]
    edges.j = [g2i[k] for k in edges.g2]
    edges.x = [1 - (x.ani * x.af)*1/(100^2) for x in eachrow(edges)]
    return edges
end


length(ARGS) != 5 && error("need exactly 5 input arguments")
data_dir = ARGS[1]
out_path = ARGS[2]
skaniS = parse(Int, ARGS[3])
skaniM = parse(Int, ARGS[4])
skaniC = parse(Int, ARGS[5])


# get gens
gens = rem_ext.(readdir("$data_dir/contigs_fixed_hdrs")) |> sort
ngen = length(gens)
g2i = Dict(k => i for (i, k) in enumerate(gens))

# read raw data
begin
    local edge_path = "$(data_dir)/skani_s$(skaniS)_m$(skaniM)_c$(skaniC).edges"
    local edge_raw = open(edge_path, "r") do f mmap(f) end;
    edges = process_raw(edge_raw, g2i)
end;

# make output matrices
ani = ones(Float16, ngen, ngen);
for e in eachrow(edges)
    ani[e.i, e.j] = ani[e.j, e.i] = e.x
end
for i in axes(ani,1); ani[i,i] = 0; end
ani2 = round.(ani, digits=19);

# write output
open(out_path, "w") do io
    println(io, ngen)
    for i in axes(ani2, 1)
        println(io, gens[i], '\t', join(ani2[:, i], '\t'))
    end
end
