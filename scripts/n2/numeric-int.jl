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


const PREFIX = "n2-plain"
const params = DEFAULT_PARAMS
const N = 2
const κ = 18e4*2π
const ϕ = 0 

function smallest_distance(z1, z2)
    min([norm([z1, z2]- col) for col in eachcol(M.z)]...)
end

function integral_curved(zg; Δ)
    R = smallest_distance(zg...)*9/10
    ϕ = range(0, 2π, 100)

    vs = map(ϕ) do p
        x = zg[1] + R*cos(p)
        y = zg[2] + R*sin(p)

        u = f(x, y, Δ)[1]
        v = f(x, y, Δ)[2]
        return (u, v)
    end

    total_angle = 0.0
    for i in 1:(length(vs)-1)
        v1 = vs[i]
        v2 = vs[i+1]

        # Signed angle between v1 and v2: atan(det, dot)
        # This handles the "cyclic shift" and branch cuts automatically
        det = v1[1]*v2[2] - v1[2]*v2[1]
        dot = v1[1]*v2[1] + v1[2]*v2[2]
        total_angle += atan(det, dot)
    end

    return round(total_angle / (2π))
end

vd(Δ; params) = begin
   l = 2*params.λ0 * (1-Δ/params.ω0)
   [-l, l]
end
vp = [0, ϕ]

f(z1, z2, Δ) = vL(
                           [z1, z2], 
                           vd(Δ; params=params),
                           vp;
                           Δ=Δ, κ=κ, params=params
              )

@load "data/n2/latest.jld2" M

function main()
    near_misses = []
    ω = Float64[]
    l = ReentrantLock()
    Threads.@threads for i in ProgressBar(collect(eachindex(M.Δ)))
        W = integral_curved(M.z[:, i]; Δ=M.Δ[i])

        if W == 0
            lock(l) do
                push!(near_misses, M.z[:, i])
                push!(ω, M.Δ[i])
            end
        end
    end
    @info "got $(length(near_misses)) near misses\n\t$(length(near_misses)/length(M.Δ))%"

    (near_misses, ω)
end

foo = main()
