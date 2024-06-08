using NaturalSort
using .Threads: @threads
using DataFrames
using PhyloTools
using PhyloTools: degree
using Colors
function getclu(path)
    l1 = readlines(path)
    l2 = split.(l1, '\t')
    clu = Dict(k => String[] for k ∈ first.(l2) |> unique)
    for x ∈ l2
        k = x[1]
        v = x[2]
        push!(clu[k], v)
    end
    return clu
end
genname(x::T) where {T<:AbstractString} = split(x, "#")[1];
ori(x::T) where {T<:AbstractString} = split(x, "#")[end];
# gene filtering functions
o1(pcdf, i1, i2) = (i1 .< pcdf.i1) .& (pcdf.i1 .< i2 .<= pcdf.i2) # starts before, ends inside
o2(pcdf, i1, i2) = (i1 .< pcdf.i1) .& (pcdf.i2 .<= i2) # starts before, ends after
o3(pcdf, i1, i2) = (pcdf.i1 .<= i1 .<= pcdf.i2) .& (pcdf.i1 .<= i2 .<= pcdf.i2) # starts inside, ends inside
o4(pcdf, i1, i2) = (pcdf.i1 .<= i1 .<= pcdf.i2) .& (pcdf.i2 .<= i2) # starts inside, ends after
f(pcdf, i1, i2) = o1(pcdf, i1, i2) .| o2(pcdf, i1, i2) .| o3(pcdf, i1, i2) .| o4(pcdf, i1, i2);
function countmap(x::Vector{T}) where {T}
    o = Dict{T,Int}()
    for v in x
        if haskey(o, v)
            o[v] += 1
        else
            o[v] = 1
        end
    end
    return o
end
function xfmt(x::Vector{T}) where {T<:Number}
    i0 = x .== 0
    o = ["" for _ in eachindex(x)]
    o[i0] .= "0"
    scaled = log10.(x[x.>0])
    expo = floor.(Int, scaled)
    rem = scaled .- expo
    scalex = round.(10 .^ rem, digits=2)
    scale = [string(x)[end-1:end] == ".0" ? string(Int(x)) : string(x) for x in scalex]
    o[.!i0] = scale .* "⋅10" .* Makie.UnicodeFun.to_superscript.(expo)
    return o
end
function colf(i::Int)
    # cols = [:lightgrey, :orange, :green, :skyblue]
    [colorant"lightgrey", colorant"orange", colorant"lightgrey", colorant"lightgrey"][i+1]
end
function scf(i::Int)
    [:transparent, :black, :transparent, :transparent][i+1]
end


tr = readnw("../../output/trees/tre.phage_10_10_10") |> midpoint_root |> mad;
PhyloTools.fixdist!(tr);
begin
    # find nonredundant genomes to keep
    local th = 5e-3
    tkp = String[]
    local md = Dict(id(n) => 0.0 for n in postwalk(tr))
    for n in postwalk(tr)
        isleaf(n) && continue
        md[id(n)] = [md[id(c)] + distance(c) for c in children(n)] |> maximum
    end
    # find edge nodes and push their leaves to tkp
    for n in prewalk(tr)
        isleaf(n) && continue
        if (md[id(n)] >= th) && (md[id(n[1])] < th)
            push!(tkp, leafnames(n)[1])
        end
    end
end


# read prot clu to create pcdf_anti DataFrame
# columns of pcdf_anti: prot, clu, i1, i2
begin
    local idt = 99
    local ppath = "../../output/protclu/phages/prots_$(idt).clu"
    local pclu = getclu(ppath)
    local pc1 = [v for (_, vs) ∈ pclu for v ∈ vs]
    local pc1s = split.(pc1, "#")
    local c1 = parse.(Int, getindex.(pc1s, 3))
    local c2 = parse.(Int, getindex.(pc1s, 4))
    local pc2 = [k for (k, vs) ∈ pclu for v ∈ vs]
    pcdf = DataFrame(gen=genname.(pc1), prot=pc1, clu=pc2, i1=c1, i2=c2)
    # local torem = filter(row -> row.i1 .>= 1.2e5, pcdf).gen |> Set
    # pcdf = filter(row -> !(row.gen in torem) && (row.gen in tkp), pcdf)
    # read mariusz's predictions and extract protein clusters flagged as antidefense
    local mar = readlines("../../data/model_predictions_20240529_with_clu.tsv") |> x -> split.(x, '\t')
    local mh = popfirst!(mar)
    local trgt = getindex.(mar, 13)
    local pns = [split(x, " - ")[2] for x in getindex.(mar, 2)]
    local i1 = parse.(Float64, getindex.(mar, 7)) |> x -> Int.(x)
    local i2 = parse.(Float64, getindex.(mar, 8)) |> x -> Int.(x)
    local x = zip(pns, i1, i2, trgt) |> collect
    local li = [findfirst(startswith(pns[i]).(pcdf.prot) .& f(pcdf, i1[i], i2[i])) for i in eachindex(pns)]
    pcdf[!, :target] .= ""
    pcdf[li, :target] = trgt
    anti_clu = pcdf.clu[li] |> unique
    pcdf_anti = filter(row -> row.clu ∈ anti_clu, pcdf)
    for k in anti_clu
        i = pcdf_anti.clu .== k
        trgts = unique(pcdf_anti.target[i]) |> x -> filter(y -> y != "", x)
        if length(trgts) == 1
            pcdf_anti.target[i] .= trgts[1]
        end
    end
end;

ag = genname.(pcdf_anti.prot) |> unique |> Set;
pcdf_ag = filter(x -> genname(x.prot) ∈ ag, pcdf);
sort!(pcdf_ag, :prot, lt=natural);
pcdf_ag.anti = pcdf_ag.clu .∈ Ref(anti_clu);
aggens = genname.(pcdf_ag.prot) |> unique
torem = filter(row -> row.i1 .>= 1.2e5, pcdf).gen |> Set
gens_all = setdiff(intersect(aggens, tkp), torem);


# prune tree
tr = PhyloTools.extract(tr, gens_all)
ln2i = Dict(k => i for (i, k) in enumerate(leafnames(tr)));


# find structural genes
begin
    local rawl = readlines("../../output/pharokka/pharokka_proteins_summary_output.tsv")
    local rs = split.(rawl, '\t')
    ctg = getindex.(rs, 5)
    local ci = .!(ctg .∈ Ref(Set(["", "other", "unknown function"])))
    # local ci = (ctg .== "DNA, RNA and nucleotide metabolism")
    ctg = [getindex.(rs, 1)[ci] ctg[ci]]
    local i = genname.(ctg[:, 1]) .∈ Ref(Set(gens_all))
    ctg = ctg[i, :]
end
early = [
    # "DNA, RNA and nucleotide metabolism",
    # "moron, auxiliary metabolic gene and host takeover",
    "integration and excision",
    # "transcription regulation"
] |> Set;
late = setdiff(unique(ctg[:, 2]), early) |> Set

# calculate distribution of antidefenses relative to phage center
ppd = Dict{AbstractString,Vector{Any}}();
mdist = Dict{AbstractString,Vector{Vector{Number}}}()
# TODO
# calculate distance from early vs late distributions as well
for g in gens_all
    ddf = filter(x -> genname(x.prot) == g, pcdf_ag)
    cf = ctg[genname.(ctg[:, 1]).==g, :]
    cfd = Dict(zip(cf[:, 1], cf[:, 2]))
    rs = Set(cf[:, 1])
    ddf.tp = Int.(ddf.anti)
    tpi = (ddf.prot .∈ Ref(rs)) .&& (ddf.tp .!= 1)
    for i in findall(tpi)
        v = cfd[ddf.prot[i]] ∈ early ? 2 : 3
        ddf.tp[i] = v
    end
    # calculate distance from mid
    begin
        mid = (ddf.i2 |> maximum) / 2
        gmid = ddf.i1 .+ ((ddf.i2 - ddf.i1) ./ 2)
        md = (gmid .- mid) ./ mid
        obs = md[ddf.anti]
        rnd = rand(md, sum(ddf.anti) * 1000)
        mdist[g] = [obs, rnd]
    end
    o = ori.(ddf.prot)
    m = [ddf.i1 ddf.i2 o ddf.anti]
    go = gene.(eachrow(m), ln2i[g])
    clr = [colf(i) for i in ddf.tp]
    sc = [scf(i) for i in ddf.tp]
    ppd[g] = [go, clr, sc, mdist]
end

# leafnames positions
tpd = PhyloTools.treepositions(tr);
lfs = getleaves(tr);
lbls = getindex.(split.(name.(lfs), '_'), 1);
lbly = [tpd[k][2] for k in id.(lfs)];
xmax = [tpd[k][1] for k in id.(lfs)] |> maximum;
lblx = ones(length(lbls)) .* xmax;

# highlight short phages clade
sps = ["1.011.O.", "1.141.A.", "1.044.O.", "1.095.O.", "1.057.O.", "1.020.O."];
tr2 = PhyloTools.extract(deepcopy(tr), sps);
ln2i2 = Dict(k => i for (i, k) in enumerate(leafnames(tr2)));
df = filter(row -> row.gen ∈ sps, pcdf);
trgts = df.target |> unique |> x -> filter(y -> y != "", x);
clrs2 = WGLMakie.Colors.distinguishable_colors(length(trgts));
t2c = Dict(trgts[i] => clrs2[i] for i in eachindex(trgts));
ppd2 = deepcopy(filter(x -> x[1] ∈ sps, ppd));
for g in sps
    dfg = filter(row -> row.gen==g, df)
    sort!(dfg, :i1);
    ti = dfg.target .!= ""
    ppd2[g][2][ti] .= [t2c[k] for k in dfg.target[ti]]
    o = ori.(dfg.prot)
    m = [dfg.i1 dfg.i2 o]
    ppd2[g][1] = gene.(eachrow(m), ln2i2[g])
end


WGLMakie.activate!(px_per_unit=2, scalefactor=2)
prnt = true
fsize = (800, 500)
begin
    # prepare
    fig = Figure(fontsize=10, size=fsize)
    gl1 = fig[1, 1] = GridLayout()
    ax1 = Axis(gl1[1, 1])
    ax2 = Axis(gl1[1, 2], xticks=[Int(x * 1e4) for x in 0:13], xtickformat=xfmt, xgridcolor=:lightgrey, xgridwidth=0.25)
    ax3 = Axis(fig[1, 1], xlabel="Viral Coordinate", ylabel="Density", title="antidefens vs rest",
        aspect=1, width=150, height=150, halign=0.9, valign=0.9, tellwidth=false, tellheight=false)
    axs = [ax1; ax2]
    linkyaxes!(ax1, ax2)
    # plot tree
    plot!(ax1, tr)
    text!(ax1, Point2f.(zip(lblx, lbly .- 0.5)); text=lbls, fontsize=8)
    # plot gene tracks
    for g in gens_all
        go, clr, sc, _ = ppd[g]
        plot!(ax2, go, color=clr, strokewidth=0.5, strokecolor=sc)
    end
    # plot antidefense position density plot
    obs = reduce(vcat, first.(values(mdist))) |> Vector{Float64}
    rnd = reduce(vcat, getindex.(values(mdist), 2)) |> Vector{Float64}
    dnsty = 0.07
    density!(ax3, rnd, color=(:lightgrey, 0.3), strokecolor=:lightgrey, strokewidth=5, boundary=(-1, 1), bandwidth=dnsty)
    density!(ax3, obs, color=(:orange, 0.1), strokecolor=:orange, strokewidth=5, boundary=(-1, 1), bandwidth=dnsty)
    # design
    limits!(ax1, 0, 1.3, nothing, nothing)
    limits!(ax2, 0, nothing, nothing, nothing)
    limits!(ax3, -1, 1, 0, nothing)
    hidespines!.(axs)
    hideydecorations!.(axs)
    hidexdecorations!(axs[1])
    colsize!(gl1, 1, Relative(0.1))
    colgap!(gl1, 1, 0)
end;
if prnt
    using CairoMakie
    CairoMakie.activate!()
    od = "../../output/pdfs/may2024"
    isdir(od) && mkpath(od)
    CairoMakie.save("$od/phage_genomes_marius_antidefense_highlight.pdf", fig, size=fsize)
end


# plot small subtree of short phages
begin
    fsize2=(300,100)
    fig2 = Figure(size=fsize2)
    ax1 = Axis(fig2[1,1])
    ax2 = Axis(fig2[1,2])
    linkyaxes!(ax1,ax2)
    plot!(ax1, tr2)
    text!(ax1, Point2f.(zip(fill(.45,ntip(tr2)), (1:ntip(tr2)) .- 0.4)); text=leafnames(tr2), fontsize=8)
    for g in sps
        go, clr, sc, _ = ppd2[g]
        plot!(ax2, go, color=clr, strokewidth=0.5, strokecolor=sc)
    end
    hidedecorations!.([ax1,ax2])
    hidespines!.([ax1,ax2])
    colsize!(fig2.layout,1,Relative(.2))
    colgap!(fig2.layout,1,0)
    limits!(ax1,-.1,1.5,nothing,nothing)
    limits!(ax2,0,nothing,nothing,nothing)
end
if prnt
    using CairoMakie
    CairoMakie.activate!()
    od = "../../output/pdfs/may2024"
    isdir(od) && mkpath(od)
    CairoMakie.save("$od/phage_genomes_marius_antidefense_highlight_subtree.pdf", fig2, size=fsize2)
end
