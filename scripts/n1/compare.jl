using DrWatson
DrWatson.quickactivate(
	"CavityEquilibria")

using Artifacts, LazyArtifacts

using Plots
using JLD2

Δo, zo = let
	@load artifact"old-n1-data/data.jld2" M
	M.Δ, M.z
end

Δ, z = let
	@load datadir("P1", "latest.jld2") M
	M.ω, M.z
end

p = plot(;
	 xlabel="Δ/2π kHz",
	 ylabel="z/μm",
	 yaxis=:log,
)
scatter!(p, Δo, zo .*1e6, label="reference simulation")
scatter!(p, Δ, z[1, :] .* 1e6, label="NLsolve, N=1")

savefig(p, joinpath(mkpath(plotsdir("n1", "comparison")), "latest.png"))

if isinteractive()
    gui(p)
end
