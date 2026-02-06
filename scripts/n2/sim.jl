using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Random

using Artifacts, LazyArtifacts

using JLD2
using Plots
if isinteractive()
    plotlyjs()
end
using ProgressBars

using CavityEquilibria
using CavityEquilibria.Parameters
using CavityEquilibria.Laser
using CavityEquilibria.RootFinding
using CavityEquilibria.Util


const PREFIX = "n2"
const params = DEFAULT_PARAMS
const N = 2
const κ = 18e4*2π
const ϕ = 0 

M = let
	@load datadir("P1", "latest.jld2") M
    M
end

Δs = unique(M.Δ)
function zsp(Δ)
	#M.z[M.Δ .== Δ]
    exp10.(range(-2, 1, length=5)) |> collect
end

z = []
ω = []
φ = []

convergence_failed = 0
outofbounds = 0
	
for Δ in ProgressBar(Δs)
	l = params.λ0 * (1- Δ/params.ω0)
	vϕ = [0, ϕ]
	vd = [-l/2, l/2]
	for zg in Iterators.product([zsp(Δ) for _ in 1:N]...)
		try
			r = find_roots(
				collect(zg), vd, vϕ; 
				Δ = Δ, κ = κ, params = params
			)
		
			if (length(z) == 0 || !any(v-> isapprox(v, r), z) && all(abs.(r) .< 1))
				push!(z, r)
				push!(ω, Δ)
                push!(φ, ϕ)
            elseif any(abs.(r) .> 1)
                global outofbounds
                outofbounds += 1
			end
		catch err
			if !(err isa ConvergenceError)
				rethrow(err)
			end
            global convergence_failed
            convergence_failed += 1
            @error "convergence error with"
		end
	end
end
@info "got $(outofbounds) solutions that were out of bounds"
@info "got $(convergence_failed) solutions that did not converge"

@show z
z = Iterators.reduce(hcat, z)

M = (Δ=ω, z=z, ϕ=φ)
let
    p = mkpath(datadir(PREFIX))
    @save joinpath(p, "latest.jld2") M
    @save joinpath(p, "$(now_nodots()).jld2") M
end

p = plot(;
  xaxis=:log,
  xlims=(10^-2, 50),
  xlabel="z1/μm",
  yaxis=:log,
  ylims=(10^-2, 50),
  ylabel="z2/μm",
  zlabel="Δ/ 2πkHz"
)
scatter!(p,
         z[1, :] .* 1e6, z[2, :] .* 1e6, ω /(2π*1e3)
)

if isinteractive()
    gui(p)
end
