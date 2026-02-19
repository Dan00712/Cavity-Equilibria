using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Logging
using LinearAlgebra

using JLD2
using Plots
gr()
default(show= false)
using ProgressBars

using CavityEquilibria.Parameters
using CavityEquilibria.RootFinding
using CavityEquilibria.Util

const PREFIX = "n2-red"
const N = 2
const φ = 0
markers = [
          :xcross,
          :cross,
          :diamond,
          :dtriangle,
          :utriangle
]
colors = [
          :blue,
          :red,
          :green,
          :orange,
          :purple,
          :yellow
]


ω, ϕ, z = let
	@load datadir("n2", "latest.jld2") M

    vd(Δ; params) = begin
        l = params.λ0 * (1-Δ/params.ω0)
        [-l, l]
    end
    params=DEFAULT_PARAMS
    f(z1, z2, Δ) = norm(vL([z1, z2], vd(Δ, params=params), [0, φ], Δ=Δ, κ=18e4*2π, params=params))

    vs = [f(z1i, z2i, Δi) for (z1i, z2i, Δi) in zip(M.z[1,:], M.z[2,:], M.Δ)]
    sel = vs .== 0

    M.Δ[sel], M.ϕ[sel], M.z[:, sel]
end

Δs = unique(ω)
Δs = Δs[round.(Int, range(1, length(Δs), length=4))]

p = plot(;
         xlabel="z1 /μm",
         ylabel="z2 /μm"
)

for (Δ, m, c) in zip(Δs, markers, colors)
    scatter!(p,
             z[1, ω .== Δ],
             z[2, ω .== Δ],
             marker=m, color=c,
             label="Δ=$(round(Δ/2π/1e3; digits=2)) 2π kHz"
   )
end

savefig(p, joinpath(
                    mkpath(plotsdir(PREFIX)), "latest.pdf"
))
savefig(p, joinpath(
                    mkpath(plotsdir(PREFIX)), "$(now_nodots()).png"
))

if isinteractive()
    gui(p)
end
