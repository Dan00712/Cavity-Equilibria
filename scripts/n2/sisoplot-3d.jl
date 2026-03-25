using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Logging
using LinearAlgebra

using JLD2
using PlotlyJS
using ProgressBars


using CavityEquilibria.Parameters
using CavityEquilibria.RootFinding
using CavityEquilibria.Util

const params = DEFAULT_PARAMS
const PREFIX = "n2-iso"
const N = 2
const κ = 18e4 * 2π

ω, ϕ, z = let
    path = if length(ARGS) >= 1
        ARGS[1]
    else
        datadir(PREFIX, "latest.jld2")
    end

	@load path M

    M.Δ, M.ϕ, M.z
end

#=trace1 = scatter(
         x=z[1, :] .* 1e6, y=z[2, :] .* 1e6, z=ω /(2π*1e3),
    label="equilibria-positions",
    #xaxis=:log,
#    xlims=(10^-2, 50),
    xlabel="z1/μm",
    #yaxis=:log,
#    ylims=(10^-2, 50),
    ylabel="z2/μm",
    zlabel="Δ/ 2πkHz"
) =#

z1 = (z[1, begin:40:end]) |> collect
z2 = (z[2, begin:40:end]) |> collect
ωs = unique(ω) |> collect

vals = [
    begin 
        l =  params.λ0 * (1 - Δ/params.ω0)
        norm(vL([z1i, z2i], [-l/2, l/2], [0, 0]; Δ=Δ, κ=κ, params=params))
    end
    for z1i in z1, z2i in z2, Δ in ωs
]

trace2 = isosurface(x = z1 .* 1e6, y = z2 .* 1e6, z = ω, 
                    value=vals,
                    isomin=-.1, isomax = 1e-5,
                    opacity=.6, surface_count=1
)
layout = Layout(
                title="surface",
                scene=(
                       xaxis=attr(title="z1/μm"),
                       yaxis=attr(title="z1/μm"),
                       zaxis=attr(title="z1/μm")
                )
)
p = plot(trace2, layout)

savename = if length(ARGS) >= 2
    ARGS[2]
else
    "latest"
end
  
savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "$(savename).html"))
savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "$(now_nodots()).html"))

if isinteractive()
    display(p)
end
