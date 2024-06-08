using PhyloTools
using NaturalSort
using ProgressMeter
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
function main1(g, faaFls, gffFls, io, lk)
    f1 = filter(contains(g), faaFls)[1]
    f2 = filter(contains(g), gffFls)[1]
    faa = readfasta(f1)
    hdrs = gff2hdr(readgff(f2), g)
    k2h = (k => filter(contains("#$(split(k,' ')[1])#"), hdrs)[1] for k in keys(faa)) |> collect
    sort!(k2h, by=x -> x[2], lt=NaturalSort.natural)
    lock(lk) do
        for (k, h) in k2h
            v = replace(faa[k], r"[\n\r]" => "")
            length(v) == 0 && @error "$k"
            println(io, '>', h, '\n', v)
        end
    end
end
function phagehdr(x, g)
    xs = split(x, "_CDS_")
    inds = match(r"(\d+)\.\.(\d+)", xs[2])
    ori = contains("complement")(xs[2]) ? "-" : "+"
    xx = [g; xs[1]; inds[1]; inds[2]; ori]
    return join(xx, "#")
end
function main2(g, faaFls, io, lk)
    f1 = filter(contains(g), faaFls)[1]
    faa = readfasta(f1)
    k2h = (k => phagehdr(k, g) for k in keys(faa)) |> collect
    sort!(k2h, by=x -> x[2], lt=NaturalSort.natural)
    lock(lk) do
        for (k, h) in k2h
            v = replace(faa[k], r"[^A-Za-z]" => "")
            length(v) == 0 && @error "$k"
            println(io, '>', h, '\n', v)
        end
    end
end
function main()
    organism = length(ARGS) > 0 ? lowercase(ARGS[1]) : "vibrio"
    faaFls = readdir("../../data/$organism/annotations/faa", join=true) |> x -> filter(endswith(".faa"), x)
    gffFls = readdir("../../data/$organism/annotations/gff", join=true) |> x -> filter(endswith(".gff3"), x)
    gens = gid.(faaFls)
    outdir = "../../output/protclu/$organism"
    !isdir(outdir) && mkpath(outdir)
    lk = ReentrantLock()
    io = open("$outdir/prots.faa", "w")
    prg = ProgressUnknown(; desc="concatenating proteomes")
    Threads.@threads for g in gens
        organism == "vibrio" ? main1(g, faaFls, gffFls, io, lk) : main2(g, faaFls, io, lk)
        next!(prg)
    end
    finish!(prg)
    close(io)
end;

main()



