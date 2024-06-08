using CairoMakie
CairoMakie.activate!()
# using WGLMakie
# WGLMakie.activate!(scalefactor=2)
using NaturalSort
using Graphs
using GraphMakie
using GraphMakie.NetworkLayout
using SparseArrays
include("helper_functions.jl")

# read prot clu
begin
    idt = ARGS[1]
    local vpath = "../../output/protclu/vibrio/prots_$(idt).clu"
    local ppath = "../../output/protclu/phages/prots_$(idt).clu"
    idt = parse(Float64, idt)
    vclu = readclu(vpath)
    pclu = readclu(ppath)
    pp2clu = Dict(v => k for (k, vs) in pclu for v in vs)
    vp2clu = Dict(v => k for (k, vs) in vclu for v in vs)
end;
# find close phages that infect a host with close relatives
begin
    thp = 5e-2
    thv = 5e-2
    close_vibs_all = rlf("../../output/close_relatives/by_phage_$(thv)_vib_10_10_10.clu")
    close_phages_all = rlf("../../output/close_relatives/by_vib_$(thp)_phage_10_10_10.clu")
    pairs = triangulate(close_phages_all, close_vibs_all)
end;


# loop over pairs to build putative interaction graph
global edges_txt = Matrix{AbstractString}(undef, 0, 2)
for p in pairs
    vibs = p[1]
    phages = p[2]
    #
    # reconstruct interactions
    cff = filter(x -> contains(phages[1])(x[2]), close_phages_all) # get lines from close_phages_all that have current phage group
    cvf = filter(x -> contains(vibs[1])(x[2]), close_vibs_all) # get lines from close_vibs_all that have current vib group
    filter!(x -> x[1] in vibs, cff)
    filter!(x -> x[1] in phages, cvf)
    # make interaction matrix
    v2i = Dict(x => i for (i, x) in enumerate(vibs))
    p2j = Dict(x => i for (i, x) in enumerate(phages))
    m = zeros(Bool, length(vibs), length(phages))
    for k in vibs
        i = v2i[k]
        x0 = filter(x -> x[1] == k, cff)
        if length(x0) == 1
            x1 = x0[1][2]
            pos, neg = split(x1, ';')
            for x2 in split(pos, ',')
                if haskey(p2j, x2)
                    j = p2j[x2]
                    m[i, j] = true
                end
            end
        end
    end
    mp = unique(m, dims=2)
    mpcnt = []
    for mpc in eachcol(mp)
        t1 = findall(mpc .== 1)
        if length(t1) > 0
            push!(mpcnt, t1)
        end
    end
    #
    # read defense system files
    v2gfam = Dict(v => Set{AbstractString}() for v in vibs)
    for vib in vibs
        lraw = split.(readlines("../../output/defensefinder/$vib/genes.tsv"), '\t')
        hdr = popfirst!(lraw)
        # v2syst[vib] = Set(getindex.(lraw,3))
        # gfs = [haskey(vp2clu,x[end]) ? vp2clu[x[end]] : x[end] for x in lraw] |> Set
        gfs = [vp2clu[x[end]] for x in lraw] |> Set
        v2gfam[vib] = gfs
    end
    core_defense = reduce(intersect, values(v2gfam)) |> collect |> sort
    var_defense = setdiff(reduce(union, values(v2gfam)), core_defense) |> collect |> sort
    # putative defenses are core defenses and variable defenses that match the infection pattern
    a2j = Dict(a => j for (j, a) in enumerate(var_defense))
    avm = zeros(Bool, length(vibs), length(var_defense))
    for (i, v) in enumerate(vibs)
        va = intersect(v2gfam[v], var_defense)
        if length(va) > 0
            for a in va
                j = a2j[a]
                avm[i, j] = true
            end
        end
    end
    avmi = zeros(Bool, size(avm, 2))
    for (i, c) in enumerate(eachcol(avm))
        for cnst in mpcnt
            if all(c[cnst] .== 1)
                avmi[i] = true
            end
        end
    end
    defenses = [core_defense; var_defense[avmi]] # final putative defense list
    #
    # get phage pangenome
    proteomes = [filter(x -> startswith(p)(x[1]), pp2clu) |> values |> unique |> Set for p in phages]
    core_phage = reduce(intersect, proteomes)
    var_phages = [setdiff(x, core_phage) for x in proteomes]
    var_phage = reduce(union, var_phages) |> collect
    var2i = Dict(var => i for (i, var) in enumerate(var_phage))
    pvm = zeros(Bool, length(var_phage), length(phages))
    for (j, phage) in enumerate(phages)
        pv = var_phages[j]
        i = [var2i[k] for k in pv]
        pvm[i, j] .= true
    end
    # find rows (protein families) in pvm that correspond to infection outcome
    pvmi = zeros(Bool, size(pvm, 1))
    for (i, r) in enumerate(eachrow(pvm))
        if r in eachrow(m)
            pvmi[i] = true
        end
    end
    pipfs = var_phage[pvmi]
    if length(pipfs) > 0
        for ad in pipfs
            for d in defenses
                global edges_txt = vcat(edges_txt, [ad d])
            end
        end
    end
end;

# build graph
nodes = [sort(unique(edges_txt[:, 1]), lt=natural); sort(unique(edges_txt[:, 2]), lt=natural)];
nad = edges_txt[:, 1] |> unique |> length
nd = edges_txt[:, 2] |> unique |> length
n2i = Dict(n => i for (i, n) in enumerate(nodes));
m = spzeros(Int, length(nodes), length(nodes))
g = SimpleGraph(length(nodes))
for e in eachrow(edges_txt)
    i1 = n2i[e[1]]
    i2 = n2i[e[2]]
    m[i1, i2] += 1
    m[i2, i1] = m[i1, i2]
    add_edge!(g, i1, i2)
end;
# lout = Stress()(g);


f= Figure(size=(400,400), fontsize=8)
ax1 = Axis(f[1:3, 1])
hidedecorations!(ax1);
hidespines!(ax1);
ax2 = Axis(f[3, 2], title="Putative viral anti-defenses\nnode degree distribution", xlabel="degree", ylabel="# of antidefenses")
nc = [:black for _ in 1:nv(g)]
nc[1:nad] .= :red
# gp = graphplot!(ax1, g, layout=lout,
gp = graphplot!(ax1, g,
    edge_color=:lightgrey, node_color=nc, node_size=5, edge_width=.1)
deg = degree(g)
hist!(ax2, deg[1:nad], bins=0:5:maximum(deg[1:nad]), strokewidth=1)
limits!(ax2, 0, 45, 0, nothing)
colsize!(f.layout, 1, Relative(0.8))
pdfod = "../../output/pdfs/association_networks"
!isdir(pdfod) && mkpath(pdfod)
CairoMakie.save("$pdfod/Extended_Data_Fig1_B_network_protclu_id_$idt.pdf", f, pt_per_unit=1)





