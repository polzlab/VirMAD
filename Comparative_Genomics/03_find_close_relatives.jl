using Leiden
using NaturalSort
using SparseArrays
using StringViews
f(x, th) = @inbounds findall(parse.(Float16, x[2:end]) .< th);
function Base.count(x::Vector{T}) where {T}
    d = Dict(k => 0 for k in x)
    for k in x
        d[k] += 1
    end
    d
end
function push2dict!(d, k, v)
    if haskey(d, k)
        push!(d[k], v)
    else
        d[k] = [v]
    end
end

function main()

    # read infection matrix
    begin
        local infpath = "../data/infection_matrix.tsv"
        local infraw = read(infpath) |> StringView
        local i1 = split(infraw, '\n')
        # ifxndesc = split(popfirst!(i1), '\t')
        # ifxnhdr = split(popfirst!(i1), '\t')
        popfirst!(i1)
        popfirst!(i1)
        ifxn = [[x[1:2]; parse(Int, x[6])] for x in split.(i1, '\t')]
        pid = Dict{String,Vector{String}}()
        vid = Dict{String,Vector{String}}()
        for x in ifxn
            if x[3] > 1
                push2dict!(pid, x[1], x[2])
                push2dict!(vid, x[2], x[1])
            end
        end
        @assert all([phage in vid[vib] for (phage,vibs) in pid for vib in vibs])
        @assert all([vib in pid[phage] for (vib,phages) in vid for phage in phages])
    end

    # read raw distance matrix
    fls = readdir("../output/trees/dist_mats", join=true)
    outdir = "../output/close_relatives"
    !isdir(outdir) && mkpath(outdir)

    for anipath in fls

        # read distance matrix
        bn = replace(basename(anipath), ".dist" => "")
        println(bn)
        pflag = startswith("phage")(bn)
        raw = read(anipath) |> StringView
        rs_pre = split(raw, '\n')       
        popfirst!(rs_pre)
        pop!(rs_pre)
        rs = split.(rs_pre, '\t')
        ii = (first.(rs) .∈ (pflag ? Ref(keys(pid)) : Ref(keys(vid)))) |> findall
        jj = [1;ii .+ 1]
        rs = [x[jj] for x in rs[ii]]
        gens = first.(rs)
        ngen = length(gens)

        # make graph adjacency matrix
        for dist_th in [1e-3, 2e-3, 5e-3, 1e-2, 2e-2, 5e-2, 1e-1, 1]

            println("\t", dist_th)
            iis = f.(rs, dist_th)
            m = spzeros(Bool, ngen, ngen)
            for i in eachindex(iis)
                @inbounds ii = iis[i]
                @inbounds m[ii, ii] .= true
            end
            lraw = filter(x->length(x)>1,Leiden.leiden(m; resolution=.1)[2])
            clu = [gens[i] for i in lraw]

            # create lines to write
            towrite = String[]
            for ks in clu
                tchk = reduce(union, [pflag ? pid[k] : vid[k] for k in ks]);
                for k in tchk
                    tchk2 = filter(x -> x ∈ ks, pflag ? vid[k] : pid[k])
                    if length(tchk2)!=length(ks)
                        pos = eltype(ks)[]
                        neg = eltype(ks)[]
                        for kk in ks
                            kk ∈ tchk2 ? push!(pos,kk) : push!(neg,kk)
                        end
                    if all(length.([pos,neg]).>0) 
                        ps = join(pos, ',')
                        ns = join(neg, ',')
                        push!(towrite, "$k\t$ps;$ns")
                    end
                    end
                end
            end

            # write output
            open("$outdir/$(pflag ? "by_vib_" : "by_phage_")$(dist_th)_$(bn).clu", "w") do io
                println(io, pflag ? "vib\tinfects;protected" : "phage\tinfected;protected")
                for l in towrite
                    println(io, l)
                end
            end

        end
    end
end

main()
