using DrWatson
DrWatson.quickactivate(
    "CavityEquilibria"
)

using Random

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


const PREFIX = "n2"
const params = DEFAULT_PARAMS
const N = 2
const κ = 18e4*2π
const ϕ = 0 

M = let
	@load datadir("P1", "latest.jld2") M
    M
end

z_sp = (copy(M.z))
Δs = unique(M.Δ)
function zsp(Δ)
    vcat(-1 .* z_sp, z_sp)
end


function main()
    z = []
    ω = []
    φ = []
    lc = ReentrantLock()
    
    convergence_failed = Threads.Atomic{Int}(0)
    outofbounds = Threads.Atomic{Int}(0)

    Threads.@threads for i in ProgressBar(1:length(Δs))
        Δ = Δs[i]
    	l = params.λ0 * (1- Δ/params.ω0)
    	vϕ = [0, ϕ]
    	vd = [-l/2, l/2]

        z_ = []
        ω_ = []
        φ_ = []

    	for zg in Iterators.product([zsp(Δ) for _ in 1:N]...)
    		try
    			r = find_roots(
    				collect(zg), vd, vϕ; 
    				Δ = Δ, κ = κ, params = params
    			)
    		
                if ((length(z_) == 0 || !any(v-> isapprox(v, r), z_)) && all(abs.(r) .< 1))
    				push!(z_, r)
    				push!(ω_, Δ)
                    push!(φ_, ϕ)
                elseif any(abs.(r) .> 1)
                    outofbounds[] += 1
    			end
    		catch err
    			if !(err isa ConvergenceError)
    				rethrow(err)
    			end
                convergence_failed[] += 1
    		end
    	end
        lock(lc) do
            append!(z, z_)
            append!(ω, ω_)
            append!(φ, φ_)
        end
    end

    @info "got $(outofbounds) solutions that were out of bounds"
    @info "got $(convergence_failed) solutions that did not converge"
    
    @show z
    z = Iterators.reduce(hcat, z)
    z = Float64.(z)

    ω = Float64.(ω)
    φ = Float64.(φ)
    
    M = (Δ=ω, z=z, ϕ=φ)
    let
        p = mkpath(datadir(PREFIX))
        @save joinpath(p, "latest.jld2") M
        @save joinpath(p, "$(now_nodots()).jld2") M
    end
    
    p = plot(;
      xaxis=:log,
      xlims=(10^-2, 50),
      xlabel="z1/μm",
      yaxis=:log,
      ylims=(10^-2, 50),
      ylabel="z2/μm",
      zlabel="Δ/ 2πkHz"
    )
    scatter!(p,
             z[1, :] .* 1e6, z[2, :] .* 1e6, ω /(2π*1e3)
    )
    
    if isinteractive()
        savefig(p, joinpath(mkpath(plotsdir(PREFIX)), "latest.html"))
        gui(p)
    end
    p, M
end

f = main()
