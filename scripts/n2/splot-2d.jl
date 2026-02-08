using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Logging

using JLD2
using Plots
gr()
default(show= false)
using ProgressBars

using CavityEquilibria.Util

const PREFIX = "n2"
const N = 2
const φ = 0

ω, ϕ, z = let
	@load datadir(PREFIX, "latest.jld2") M
    M.Δ, M.ϕ, M.z
end

p = plot(;
    xlabel="Δ /2πkHz",
    yaxis=:log,
    ylims=(10^-2, 50),
    ylabel="z/μm",
    legend=:bottomleft,
)

scatter!(p,
         ω[φ .== ϕ] /(2π*1e3),
         z[1, φ .== ϕ] .* 1e6,  
         label="z1"
)
scatter!(p,
         ω[φ .== ϕ] /(2π*1e3),
         z[2, φ .== ϕ] .* 1e6,  
         label="z2"
)

savefig(p, joinpath(
                    mkpath(plotsdir(PREFIX)), "latest.pdf"
))
savefig(p, joinpath(
                    mkpath(plotsdir(PREFIX)), "$(now_nodots()).png"
))

if isinteractive()
    gui(p)
end
