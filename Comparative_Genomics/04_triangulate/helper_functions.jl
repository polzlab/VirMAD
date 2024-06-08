using StringViews
using NaturalSort

# structs
struct Res
    vibs::AbstractVector
    phages::AbstractVector
    m::AbstractArray
    defense_var::AbstractVector
    phageprot_var::AbstractVector
    putative_defense::AbstractVector
    defense_mat::AbstractArray
    putative_phageprot::AbstractVector
    phageprot_mat::AbstractArray
    core_defenses::AbstractVector
end


# oneliners
rlf(x) = split.(readlines(x), '\t')
ngen(x) = first.(split.(x, "#")) |> unique |> length
gs(x) = split(x, r"[,;]") |> sort;
ppf(x, pp) = [any(startswith(p).(x[2])) for p in pp] |> all
pnf(x, pn) = [any(startswith(p).(x[2])) for p in pn] |> any

function triangulate(close_phages_all, close_vibs_all)
    o = Set{Tuple{Vector{String},Vector{String}}}()
    for p in close_phages_all
        vib = p[1]
        x = p[2]
        vs = filter(y -> vib in gs(y[2]), close_vibs_all)
        if length(vs) > 0
            for vibs in vs
                push!(o, (gs(vibs[2]), gs(x)))
            end
        end
    end
    return o
end

function readclu(path)
    l1 = readlines(path)
    l2 = split.(l1, '\t')
    clu = Dict(k => String[] for k in first.(l2) |> unique)
    for x in l2
        k = x[1]
        v = x[2]
        push!(clu[k], v)
    end
    return clu
end


function check_dim(xm, m, dim)
    f(x, y) = all(collect(x) .== collect(y))
    od = dim == 1 ? 2 : 1
    tkp = zeros(Bool, size(xm, dim))
    for (i, x) in enumerate(eachslice(xm, dims=dim))
        tkp[i] = [(f(y, x) | f(y, 1 .- x)) for y in eachslice(m, dims=dim)] |> any
    end
    return tkp
end


function extend_interaction_matrix(defenses, pipfs, vclu, pclu, infections)
    dprots = [vclu[k] for k in defenses]
    pprots = [pclu[k] for k in pipfs]
    bg = reduce(union, [first.(split.(prots, "#")) for prots in dprots])
    pg = reduce(union, [first.(split.(prots, "#")) for prots in pprots])
    ii = (first.(infections) .∈ Ref(pg)) .&& (getindex.(infections, 2) .∈ Ref(bg))
    ifx2 = infections[ii,:]
end

stringview_lines(p) = split(read(p) |> StringView, r"[\n\r]+")

function read_infection_matrix(p)
    raw = stringview_lines(p)
    rs = split.(raw, '\t')
    popfirst!(rs)
    popfirst!(rs)
    return [x[1:2] for x in rs if length(x)>5 && parse(Int,x[6])>1]
end

function push2dict!(d::Dict{KT,Vector{VT}},k::KT,v::VT) where {KT,VT}
    if haskey(d,k)
        push!(d[k],v)
    else
        d[k] = v
    end
    return nothing
end

