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


function local_minima(A)
    results = []
    for i in ProgressBar(2:size(A,1)-1)
        for j in 2:size(A,2)-1
            for k in 2:size(A,3)-1
                val = A[i,j,k]
                if val < A[i-1,j,k] && val < A[i+1,j,k] &&
                   val < A[i,j-1,k] && val < A[i,j+1,k] 
                    push!(results, (i,j,k))
                end
            end
        end
    end
    return results
end

z1, z2, Δs = let
    @load datadir("P1", "latest.jld2") M

    z1 = range(25e-6, 30e-6, 500)
    z2 = copy(z1)
    Δ = range(extrema(M.Δ)..., 100)

    vd(Δ; params) = begin
        l = 2*params.λ0 * (1-Δ/params.ω0)
        [-l, l]
    end
    vp = [0, ϕ]

    f(z1, z2, Δ) = norm(vL(
                           [z1, z2], 
                           vd(Δ; params=params),
                           vp;
                           Δ=Δ, κ=κ, params=params
              ))
    fs = f.(z1, z2', reshape(Δ, 1, 1, :))
    minimas = local_minima(fs)
    z = []
    ω = Float64[]
    for (i,j,k) in minimas
        push!(z, [z1[i], z2[j]])
        push!(ω, Δ[k])
    end
    z = Iterators.reduce(hcat, z)
    z = Float64.(z)

    z[1, :], z[2,:], ω
end


function main()
    z = []
    ω = []
    φ = []
    lc = ReentrantLock()
    
    convergence_failed = Threads.Atomic{Int}(0)
    outofbounds = Threads.Atomic{Int}(0)

    for (z1i, z2i, Δ) in ProgressBar(zip(z1, z2, Δs))
        zg = [z1i, z2i]
    	l = params.λ0 * (1- Δ/params.ω0)
    	vϕ = [0, ϕ]
    	vd = [-l/2, l/2]

        z_ = []
        ω_ = []
        φ_ = []
        try
    		r = find_roots(
    			zg, vd, vϕ; 
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
        if length(z_) == 0
            continue
        end
        lock(lc) do
            append!(z, z_)
            append!(ω, ω_)
            append!(φ, φ_)
        end
    end

    @info "got $(outofbounds) solutions that were out of bounds"
    @info "got $(convergence_failed) solutions that did not converge"
    
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
 #     xlims=(10^-2, 50),
      xlabel="z1/μm",
      yaxis=:log,
#      ylims=(10^-2, 50),
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
