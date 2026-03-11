using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Random
using LinearAlgebra

using Artifacts, LazyArtifacts

using JLD2
using Plots
plotlyjs()
using ProgressBars

using CavityEquilibria
using CavityEquilibria.Parameters
using CavityEquilibria.Laser
using CavityEquilibria.RootFinding
using CavityEquilibria.Util


const PREFIX = "n2-heat"
const N = 2
const κ = 18e4*2π
const ϕ = 0
const params = DEFAULT_PARAMS

const z1 = let
    μ = 5.156e-6
    extending = 5e-6
    range(μ-extending, μ+extending, 1000)
end
const z2 = let
    μ = -5.533e-6
    extending = 5e-6
    range(μ-extending, μ+extending, 1000)
end
const Δ =  174.3030303030303 *2π

vd(Δ; params) = begin
    l = params.λ0 * (1-Δ/params.ω0)
    [-l, l]
end
vL2(z1, z2, Δ, params) = vL(
                           [z1, z2], 
                           vd(Δ; params=params),
                           [0, ϕ]; Δ=Δ, κ=κ, params=params)

L(z1, z2, Δ, params) = norm(vL2(z1, z2, Δ, params))

function main()
    f(z1, z2) = L(z1, z2, Δ, params)

    vs = Matrix{Float64}(undef, length(z1), length(z2))
    Threads.@threads for j in ProgressBar(axes(vs, 2))
        for i in axes(vs, 1)
            @inbounds vs[i, j] = f(z1[i], z2[j])
        end
    end

    p = heatmap(z1.*1e3, z2.*1e3, vs, label="Case: z₁ ∼ z₂")
    savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "latest.pdf"))
    savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "$(now_nodots()).png"))

    if isinteractive()
        gui(p)
    end

    p, vs
end

f = main()
