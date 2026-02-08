using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Logging

using JLD2
using Plots
default(show= false)
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
ps = []
titles= []
for (i, (l, ϕi)) in ProgressBar(enumerate(zip(Φs, Φ)))
    p_ = plot(;
        xlabel="Δ /2πkHz",
        yaxis=:log,
        ylims=(10^-2, 50),
        ylabel="z/μm",
        legend=:bottomleft,
        title="ϕ = $l"
    )

    scatter!(p_,
             ω[ϕi .== ϕ] /(2π*1e3),
             z[1, ϕi .== ϕ] .* 1e6,  
             label="z1"
    )
    scatter!(p_,
             ω[ϕi .== ϕ] /(2π*1e3),
             z[2, ϕi .== ϕ] .* 1e6,  
             label="z2"
    )

    savefig(p_, joinpath(
                        mkpath(plotsdir(PREFIX, "ind")), "latest-$(i).pdf"
    ))

    push!(ps, p_)
end

p = plot(ps...;
         layout=(4, 1),
         size=(400, 1200),
         subplot_titles = titles,
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
