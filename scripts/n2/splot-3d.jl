using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Logging

using JLD2
using Plots
plotlyjs()
default(show= false)
using ProgressBars

using ScatteredInterpolation

using CavityEquilibria.Util

const PREFIX = "n2"
const N = 2
const φ = 0

ω, ϕ, z = let
	@load datadir(PREFIX, "latest.jld2") M

    selector = M.z[1,:] .> 0 .&& M.z[2, :] .> 0
    M.Δ[selector], M.ϕ[selector], M.z[:, selector]
    #M.Δ, M.ϕ, M.z
end

itp = interpolate(Polyharmonic(2), z, ω)
z1_ = (range((extrema(z[1, :]))..., 200))
z2_ = (range((extrema(z[2, :]))..., 200))
Ds = [evaluate(itp, [z1i, z2i])[1] for z1i in z1_, z2i in z2_]

p = plot(;
    xaxis=:log,
#    xlims=(10^-2, 50),
    xlabel="z1/μm",
    yaxis=:log,
#    ylims=(10^-2, 50),
    ylabel="z2/μm",
    zlabel="Δ/ 2πkHz"
)
scatter!(p,
         z[1, :] .* 1e6, z[2, :] .* 1e6, ω /(2π*1e3),
         label="equilibria-positions"
)
#=surface!(p, z1_ .* 1e6, z2_ .* 1e6, Ds ./(2π*1e3),
       label="interpolated surface",
       colorbar=false
      )
      =#

  
savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "latest.html"))
savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "$(now_nodots()).html"))

if isinteractive()
    gui(p)
end
