#!/usr/bin/env julia
using DrWatson
HOME = ENV["HOME"]
DrWatson.quickactivate(
	joinpath(HOME, "projects", "bac", "CavityEquilibria"),
	"CavityEquilibria")

using Random

using JLD2
using Plots
if isinteractive()
    plotlyjs()
end

p = ARGS[1]
@load p d

p = plot(;
	 zlabel="Δ/ 2π kHz",
	 xlabel="z1/μm",
	 ylabel="z2/μm",
)


styles = [
          Dict(:color=>:blue),
          Dict(:color=>:orange),
          Dict(:color=>:red),
          Dict(:color=>:green),
]

for ((vϕ, M), style) in zip(d, styles)
    ω = M.Δ
    z = M.z
    scatter!(p,
		 z[1, :], # x
         z[2, :], #y
         ω,
         label="ϕ=$(round(vϕ[2]/π*10)/10)π"
    )
end

if isinteractive()
    gui(p)
end
