using FASTX
using LinearAlgebra
using Statistics
using Plots
using Random
using LaTeXStrings
using UMAP
using DelimitedFiles

io = open("proteins.fasta", "r")
reads = [(description(rec), sequence(rec)) for rec in FASTAReader(io)]
headers = first.(reads)
sequences = last.(reads)
n_seqs = length(sequences)
species = [occursin("YEAST", h) ? -1 : 1 for h in headers]
amino_acids = mapreduce(unique, union, sequences) |> sort!

const N = 50_000
hdv() = rand((-1,1), N)
shift(x, k=1) = circshift(x, k)
import Base.bind
bind(xs::Vector{Int}...) = reduce(.*, xs)
const aa_hdvs = Dict(aa=>hdv() for aa in amino_acids)
bind([hdv(), hdv(), hdv()]...)

const trimer_hdvs = Dict(aa1 * aa2 * aa3 => 
					bind(aa_hdvs[aa1], shift(aa_hdvs[aa2]), shift(aa_hdvs[aa3], 2))
			for aa1 in amino_acids for aa2 in amino_acids for aa3 in amino_acids)


function embed_sequences(sequences)
	# preallocate an empty matrix
	w=3
	hdvs = zeros(Int, length(sequences), N)
	for (i, seq) in enumerate(sequences)
		v = @view hdvs[i,:]  # ref to hdv i
		for pos in 1:length(seq)-(w-1)
			trimer = seq[pos:pos+(w-1)]
			v .+= trimer_hdvs[trimer]
			#v .+= bind([aa_hdvs[seq[j]] for j in pos:pos+w]...)
			#v .+= bind([shift(aa_hdvs[trimer[j]],j) for j in 1:length(trimer)]...)
			#println(trimer)
			#println([seq[j] for j in pos:pos+w])
		end
		v .= sign.(v)
	end
	return hdvs
end

Xseq = embed_sequences(sequences);

# Select D dimensions that are more variable
D=20000
m = mean(Xseq,dims=1)[:]
s = std(Xseq, dims=1)[:]
fig = scatter(m,s, xlabel=L"\mu", ylabel=L"\sigma", label="", xtickfontsize=15, ytickfontsize=15, guidefontsize=20)

use = sortperm(abs.(m))[1:D]
X = Xseq[:, use]

shuffled = shuffle(1:999)
X = X[shuffled, :]
y = species[shuffled]

writedlm("proteins_X", X)
writedlm("proteins_y", y)

# now we create a umap embedding for visualization using UMAP
K = Xseq[shuffled,:]*Xseq[shuffled,:]' / D
embedding = UMAP.umap(1 .- K, 2; n_neighbors=15, min_dist=0.1, metric=:precomputed)'
embedding = (embedding .- minimum(embedding, dims=1))./(maximum(embedding, dims=1).-minimum(embedding, dims=1))
scatter(embedding[:,1], embedding[:,2], markersize=10, c=y, alpha=0.3)
writedlm("proteins_UMAP", embedding)