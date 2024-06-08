using ProgressMeter
using .Threads: @threads

# functions
bn(x) = basename(x)[1:end-3]
function phanotate(p1,od)
    cmd1 = `phanotate.py -o $od/gbk/$(bn(p1)).gbk -f genbank -l 30 $p1`
    cmd2 = `genbank.py $od/gbk/$(bn(p1)).gbk -f gff -o $od/gff/$(bn(p1)).gff`
    cmd3 = `genbank.py $od/gbk/$(bn(p1)).gbk -f faa -o $od/faa/$(bn(p1)).faa`
    cmd4 = `genbank.py $od/gbk/$(bn(p1)).gbk -f fna -o $od/fna/$(bn(p1)).fna`
    run(cmd1)
    run(cmd2)
    run(cmd3)
    run(cmd4)
    return nothing
end

# get file paths
gd = "../../data/phages/contigs_fixed_hdrs"
fls = readdir(gd, join=true)

# define output dirs
od = "../../data/phages/annotations"
gbd = "$od/gbk"
gffd = "$od/gff"
fnad = "$od/fna"
faad = "$od/faa"
ods = [od;gbd;gffd;fnad;faad]
[!isdir(x) && mkpath(x) for x in ods];

# run phanotate
prg = Progress(length(fls); desc="running phanotate");
@threads for f in fls
    phanotate(f,od)
    next!(prg)
end
finish!(prg)
