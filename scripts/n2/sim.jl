#!/usr/bin/env julia
using DrWatson
HOME = ENV["HOME"]
DrWatson.quickactivate(
    joinpath(HOME, "projects", "bac", "CavityEquilibria"),
    "CavityEquilibria",
)

using Random

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

const Φ = [0]

const N = 2
const sz = let
    VALS = 15
    Random.seed!(0)

    #range is [-2, 2]
    logsz = rand(VALS) .* 4 .- 2
    1e-6 .* vcat(exp10.(logsz), -1 .* exp10.(logsz))
end
const vz = Iterators.product([sz for _ = 1:N]...)
const d = 1.5e-6    # m === 1.5 μm
const κ = 18e4 * 2π
const Δ = let
    # kc * d/2 = 2 π
    kc = 4π/d
    # c = ωc/kc
    ωc = params.c * kc
    
    ωc - params.ω0
end
const vd = [-d/2, d/2]


ϕ, z = let
    ϕ = []
    z = []
    for ϕ_ in ProgressBar(Φ)
        for zg in vz
            try
                r = find_roots(collect(zg), vd, [0, ϕ_]; Δ = Δ, κ = κ, params = params)

                if (length(z) == 0 || !any(v -> isapprox(v, r), z)) && all(abs.(r) .< 1)
                    push!(z, r)
                    push!(ϕ, ϕ_)
                end
            catch err
                if !(err isa ConvergenceError)
                    rethrow(err)
                end
            end
        end
    end

    z = Iterators.reduce(hcat, z)

    ϕ, z
end

const PREFIX="n2-single"
let
    p = joinpath(datadir(PREFIX) |> mkpath, "latest.jld2")
    M = (ϕ=ϕ, z = z)
    @save p M
end

p = plot(;
#    zlabel = "Δ/ 2π kHz",
    xlabel = "z1/μm",
    ylabel = "z2/μm",
    #yscale=:log,
    #xscale=:log,
)

scatter!(
    p,
    z[1, :] .* 1e6, # x
    z[2, :] .* 1e6, # y
#    ϕ;
#    markersize = 5,
)
savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "latest.html"))
savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "latest.pdf"))
savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "$(now_nodots()).png"))

if isinteractive()
    gui(p)
end
