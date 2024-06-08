using StringViews
using Mmap
using Base.Threads
using NaturalSort
using PhyloTools
function split_lines(fn)
    raw = read(fn) |> StringView
    raw_lines = split(raw, r"[\n\r]+")
    raw_lines_split = split.(raw_lines, '\t')
    filter!(x -> length(x) > 1, raw_lines_split)
    return raw_lines_split
end
gid(x)::String = basename(x)[1:end-4]
function readgff(f)
    l = filter(x -> !startswith("#")(x) && contains("\tCDS\t")(x), readlines(f))
    return split.(l, '\t')
end
function make_hdr(x, g)
    x0 = split(x[end], ";")
    gi = contains("gene=").(x0)
    xs = replace.(sum(gi) > 0 ? x0[[1, 2, findfirst(gi)]] : x0[1:2], r"(^[^=]+=)" => "")
    xs2 = replace.(xs, "%2C" => ",")
    xs3 = replace.(xs2, r"\s+" => "_")
    o = [g; x[[1, 4, 5, 7]]; xs3]
    return join(o, "#")
end
gff2hdr(gff, g) = make_hdr.(gff, g)
function genes2hdrs(g, gffFls, genes)
    f2 = filter(endswith("/$g.gff3"), gffFls)[1]
    hdrs = gff2hdr(readgff(f2), g)
    k2h = Dict(k => filter(contains("#$(split(k,' ')[1])#"), hdrs)[1] for k in genes)
    return k2h
end
hdrchk(x) = eachmatch(r"#contig_\d+#",x , overlap=true) |> collect |> length

# get file paths
begin
    local defdir = "../../output/defensefinder"
    dfdirs = readdir(defdir, join=true)
    local annodir = "../../data/vibrio/annotations"
    local gffdir = "$annodir/gff"
    gffFls = readdir(gffdir, join=true)
end

Threads.@threads for d in dfdirs
# for d in dfdirs
    g = basename(d)
    fn = filter(endswith("defense_finder_genes.tsv"), readdir(d, join=true))[1]
    raw_lines_split = split_lines(fn)
    ln1 = popfirst!(raw_lines_split)
    genes = getindex.(raw_lines_split, 2)
    hdrs = genes2hdrs(g,gffFls,genes)
    @assert all([hdrchk(x)==1 for x in values(hdrs)]) "broke headers in $d"
    open("$d/genes.tsv", "w") do io
        println(io, join(ln1,'\t'))
        for i in eachindex(raw_lines_split)
            l = raw_lines_split[i]
            println(io, join([l;hdrs[l[2]]],'\t'))
        end
    end
end
