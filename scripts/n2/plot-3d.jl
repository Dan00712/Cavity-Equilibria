using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Logging

using JLD2
using Plots
if isinteractive()
    plotlyjs()
end
using ProgressBars

using CavityEquilibria.Util

const PREFIX = "n2/multi-phi"
const N = 2
const Φ = [0, π/2, -π/2, π]  
const Φs = ["0", "π/2", "-π/2", "π"]

ω, ϕ, z = let
	@load datadir(PREFIX, "latest.jld2") M
    M.Δ, M.ϕ, M.z
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

for (l, ϕi) in ProgressBar(zip(Φs, Φ))
    scatter!(p,
             z[1, ϕi .== ϕ] .* 1e6, z[2, ϕi .== ϕ] .* 1e6, ω[ϕi .== ϕ] /(2π*1e3),
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
