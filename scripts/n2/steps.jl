using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

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


const PREFIX = "steps"
const params = DEFAULT_PARAMS
const N = 2
const κ = 18e4*2π
const ϕ = 0 
const dt = 1e-4


z1, z2, Δ = let
    i = ARGS[1] |> x-> parse(Int, x)
    fname = datadir("n2", "latest.jld2") 
    @load  fname M
    M.z[:, i]..., M.Δ[i]
end

l = params.λ0 * (1- Δ/params.ω0)

function main()
    η = [z1, z2, 0, 0, α_c_eq([z1, z2], [-l/2, l/2], [0,0];
                              Δ=Δ, κ=κ, params=params)]
    t = Float64[]
    q = [η[1:2]]
    for _ in range(start=0, step=dt, length=1_000)
        η += dηodt(η; Δ=Δ, κ=κ, params=params) .* dt
        push!(q, η[1:2])
    end

    q = real.(Iterators.reduce(hcat, q))

    p = scatter(q[1, 2:end] .* 1e6, q[2, 2:end] .* 1e6;
        label="curve of [z1=$(round(z1 .* 1e6; digits=2))μm, z2=$(round(z2 .* 1e6; digits=2))μm]",
         xlabel="z1/μm",
         ylabel="z2/μm",
         xlims=(-30, 30),
         ylims=(-30, 30)
    )
    scatter!(p, [z1*1e6], [z2*1e6], label="initial starting point")

    savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "latest.png"))
    if isinteractive()
        gui(p)
    end
    p, q
end

f = main()
