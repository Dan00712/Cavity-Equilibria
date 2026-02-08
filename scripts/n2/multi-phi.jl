using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Random
using Logging

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


const params = DEFAULT_PARAMS
const PREFIX = "n2/multi-phi"
const N = 2
const κ = 18e4*2π
const Φ = [0, π/2, -π/2, π]  
const Φs = ["0", "π/2", "-π/2", "π"]

M = let
	@load datadir("P1", "latest.jld2") M
    M
end

Δs = unique(M.Δ)
function zsp(Δ)
    M.z #[Δ .== M.Δ]
end

z = []
ω = []
φ = []
lc = ReentrantLock()

@info "running Φ values"
for ϕ in Φ
    @info "ϕ = $(round(ϕ, digits=2))"
    Threads.@threads for Δ in ProgressBar(Δs)
        z_ = []
        ω_ = []
        φ_ = []

        l = params.λ0 * (1- Δ/params.ω0)
        vϕ = [0, ϕ]
        vd = [-l/2, l/2]
        for zg in Iterators.product([zsp(Δ) for _ in 1:N]...)
            try
                r = find_roots_log(
                    collect(zg), vd, vϕ;
                    Δ = Δ, κ = κ, params = params
                )

                if (length(z) == 0 || !any(v-> isapprox(v, r), z_) && all(abs.(r) .< 1))
                    push!(z_, r)
                    push!(ω_, Δ)
                    push!(φ_, ϕ)
                end
            catch err
                if !(err isa ConvergenceError)
                    rethrow(err)
                end
            end
            lock(lc) do
                append!(z, z_)
                append!(ω, ω_)
                append!(φ, φ_)
            end
        end
    end
end

@info "reformatting data"
z = Iterators.reduce(hcat, z)
z = Float64.(z)

ω = Float64.(ω)
φ = Float64.(φ)

@info "storing data"
M = (Δ=ω, z=z, ϕ=φ)
let
    p = mkpath(datadir(PREFIX))
    @save joinpath(p, "latest.jld2") M
    @save joinpath(p, "$(now_nodots()).jld2") M
end

@info "plotting data"
p = plot(;
  xaxis=:log,
  xlims=(10^-2, 50),
  xlabel="z1/μm",
  yaxis=:log,
  ylims=(10^-2, 50),
  ylabel="z2/μm",
  zlabel="Δ/ 2πkHz"
)

for (l, ϕ) in ProgressBar(zip(Φs, Φ))
    scatter!(p,
             z[1, φ .== ϕ] .* 1e6, z[2, φ .== ϕ] .* 1e6, ω[φ .== ϕ] /(2π*1e3),
             label=l
    )
end

savefig(p, joinpath(
                    mkpath(plotsdir(PREFIX)), "latest.html"
))
savefig(p, joinpath(
                    mkpath(plotsdir(PREFIX)), "$(now_nodots()).png"
))

if isinteractive()
    gui(p)
end
