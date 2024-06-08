using StringViews
using Mmap
using PhyloTools
using SparseArrays
# using WGLMakie
using CairoMakie
CairoMakie.activate!()
import CairoMakie: save as cairosave
rng(x) = x[1]:x[2]
function mat2rng(x)
    o = Vector{AbstractUnitRange}(undef, size(x, 1))
    Threads.@threads for i in axes(x, 1)
        o[i] = rng(x[i, :])
    end
    return o
end
splt(x) = split(x, '\t');
fix_name(x) = replace(x, r"_.+$" => "")
function readtree(p, tkp)
    tkp2 = fix_name.(tkp)
    tmp = readnw(p) |> midpoint_root
    l = filter(x -> !(fix_name(name(x)) in tkp2), getleaves(tmp))
    tree = PhyloTools.prune!(tmp, l)# |> mad;
    for n in postwalk(tree)
        n.data.distance = abs(n.data.distance)
    end
    return mad(tree)
end
function xf(x, d)
    xs = split(x, r"[,;]")
    xi = [d[k] for k in xs]
    return extrema(xi)
end

begin # read data
    # read phage and vib isolate lists
    vibs = readlines("../../data/vib_isolates.txt")
    phages = readlines("../../data/phage_isolates.txt")
    #
    # read vibrio tree
    vib_tree = readtree("../../output/trees/tre.vib_10_10_10", vibs)
    vibs = leafnames(vib_tree)
    vib2ind = Dict(k => i for (i, k) in enumerate(vibs))
    vtp = PhyloTools.treepositions(vib_tree)
    vib2coord = Dict(name(x) => vtp[id(x)][2] for x in getleaves(vib_tree))
    #
    # read phage tree
    phage_tree = readtree("../../output/trees/tre.phage_10_10_10", phages)
    phages = leafnames(phage_tree) |> x -> replace.(x, r"_.+$" => "")
    phage2ind = Dict(k => i for (i, k) in enumerate(phages))
    ptp = PhyloTools.treepositions(phage_tree)
    phage2coord = Dict(fix_name(name(x)) => ptp[id(x)][2] for x in getleaves(phage_tree))
    #
    # read interaction matrix
    p = "../../data/infection_matrix.tsv"
    raw = open(p, "r") do f
        mmap(f)
    end
    nls = findall(raw .== UInt8('\n'))
    c1 = [1; nls .+ 1]
    c2 = [nls .- 1; length(raw)]
    intervals = [c1 c2]
    rngs = mat2rng(intervals)
end;

# fill matrix
m = zeros(Bool, length(vibs), length(phages));
for ind in 3:size(intervals, 1)
    d = StringView(raw[rngs[ind]]) |> splt
    if haskey(vib2ind, d[2]) && haskey(phage2ind, d[1])
        i = vib2ind[d[2]]
        j = phage2ind[d[1]]
        m[i, j] = true
    end
end

# diversity params and figsize
sz = (400, 400)
v_th = 5e-2
p_th = 5e-2
begin
    # initialize figure
    fig = Figure(size=sz, fontsize=8)
    ax1 = Axis(fig[1, 1]) # vib tree
    ax2 = Axis(fig[1, 2], aspect=1) # cross-infection matrix
    ax3 = Axis(fig[2, 2]) # phage tree
    # actual plotting
    sm = sparse(m')
    sx, sy, _ = findnz(sm)
    markers = [Rect(i-.5, j-.5, 1, 1) for (i, j) in zip(sx, sy)]
    poly!(ax2, markers, color=(:black,.5), strokecolor=:black, strokewidth=.1)
    vt2 = deepcopy(vib_tree)
    plot!(ax1, vt2; linewidth=0.2)
    plot!(ax3, phage_tree, up=true, linewidth=0.2)
    # axes
    linkyaxes!(ax2, ax1)
    linkxaxes!(ax2, ax3)
    hidedecorations!.([ax1; ax2; ax3])
    hidespines!.([ax1; ax2; ax3])
    rowgap!(fig.layout, 1, 0.0)
    rowsize!(fig.layout, 2, Relative(0.1))
    colsize!(fig.layout, 1, Relative(0.1))
    colgap!(fig.layout, 1, 0.0)
    tightlimits!.([ax1; ax2; ax3])
    # highlight close phages that differentially infect a host
    ppath = "../../output/close_relatives/by_vib_$(p_th)_phage_10_10_10.clu";
    yrawaf = splt.(readlines(ppath));
    yraw = filter!(x -> haskey(vib2ind, x[1]), yrawaf)
    y = [vib2coord[x[1]] for x in yraw];
    x = [xf(x[2], phage2coord) for x in yraw];
    rects = Vector{Rect}();
    for i in eachindex(y)
        x0 = x[i][1] - 0.5
        y0 = y[i] - 0.5
        w = x[i][2] - x[i][1] + 1
        h = 1
        push!(rects, Rect(x0, y0, w, h))
    end;
    lines!.(ax2, rects, linewidth=.5, color=:red);
    # highlight close hosts that are differentially protected from a virus
    vpath = "../../output/close_relatives/by_phage_$(v_th)_vib_10_10_10.clu";
    yrawaf = splt.(readlines(vpath));
    yraw = filter!(x -> haskey(phage2ind, x[1]), yrawaf)
    x = [phage2coord[x[1]] for x in yraw];
    y = [xf(x[2], vib2coord) for x in yraw];
    rects = Vector{Rect}();
    for i in eachindex(y)
        x0 = x[i] - 0.5
        y0 = y[i][1] - 0.5
        w = 1
        h = y[i][2] - y[i][1] + 1
        push!(rects, Rect(x0, y0, w, h))
    end;
    lines!.(ax2, rects, linewidth=.5, color=:orange);
end;
#
# save plot
ofn = "../../output/pdfs/association_networks/Extended_Data_Fig1_A_nahant_v$(v_th)_p$(p_th).pdf"
cairosave(ofn, fig; size=sz, pt_per_unit=1)




