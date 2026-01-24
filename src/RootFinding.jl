module RootFinding

using NLsolve

using ..Parameters
using ..Laser

export vL, find_roots, find_roots_log10, ConvergenceError

function δc(vz, vd; Δ, params::SystemParams)
    ħ = params.ħ
    α = params.α

    Δ -α/ħ * sum(Ec(0, yi, zi; Δ=Δ, params=params)^2 for(yi, zi) in zip(vd, vz))
end

function α_c_eq(vz, vd, vϕ; Δ, κ, params::SystemParams)
    ħ = params.ħ

    α = params.α

    δc_ = δc(vz, vd; Δ=Δ, params=params)
    α/ħ * sum(Ec(0, yi, zi; Δ=Δ, params=params)*Et(0, 0, zi; params=params)*exp(im*ϕi) for(yi, zi, ϕi) in zip(vd, vz, vϕ))/(δc_ - im*κ/2)
end

function vL(vz, vd, vϕ; Δ, κ, params::SystemParams)
    c = params.c
    k0 = params.k0
    ω0 = params.ω0
    kc = let
        ωc = Δ + ω0
        ωc/c
    end
    Wc = params.Wc
    zR = params.zR

    fz(z) = (im*(k0 - zR/(z^2 + zR^2)) - z/(z^2 + zR^2))
    gz(y, z) = -2*z/Wc^2 
    α_c = α_c_eq(vz, vd, vϕ; Δ=Δ, κ=κ, params=params)
    cα_c = conj(α_c)

    [begin
          t1 = (
                cα_c * Ec(0, yi, zi; Δ=Δ, params=params) 
                + conj(Et(0, 0, zi; params=params))*exp(-im*ϕi)
               )
          t2 = (
                α_c * gz(yi, zi)*Ec(0,yi, zi; Δ=Δ, params=params)
                + fz(zi)*Et(0, 0, zi; params=params) * exp(im * ϕi)
          )
          real(t1*t2)
     end
    for (yi, zi, ϕi) in zip(vd, vz, vϕ)]
end


@kwdef struct ConvergenceError <: Exception
    msg::String
    result::NLsolve.SolverResults
end

Base.showerror(io::IO, e::ConvergenceError) = print(io, "ConvergenceError: ", e.msg)


"""
    find_roots(zguess, vd, vϕ; Δ, κ, params)

Find roots of the optical force equation using nonlinear solving.

# Arguments
- `zguess::AbstractVector`: Initial guess for z positions (in meters; with default params)
- `vd::AbstractVector`: spacing of the tweezers
- `vϕ::AbstractVector`: phase of the tweezers in reference to a global phase
- `Δ`: Detuning angular frequency 
- `κ`: damping parameter
- `params::SystemParams`: parameters of the system
- `scale`: scaling factor for root finding (zguess gets scaled by scale and the root finding function by 1/scale)

# Returns
- `Vector{Float64}`: Root positions in meters

# Throws
- `ConvergenceError`: If the nonlinear solver fails to converge

# Examples
```julia
zguess = [1e-6, 2e-6, 3e-6]
vd = [0.0, 0.1, 0.2]
vϕ = [0.0, π/4, π/2]
roots = find_roots(zguess, vd, vϕ; Δ=1e6, κ=0.5, params=DEFAULT_PARAMETER)
```
"""
function find_roots(zguess, vd, vϕ; Δ, κ, params::SystemParams, scale=1e6)
    @assert zguess isa AbstractVector
    @assert length(zguess) == length(vd)
    @assert length(zguess) == length(vϕ)

    vL! = (F, zg)-> F .= vL(zg ./scale, vd, vϕ; Δ=Δ, κ=κ, params=params)
    root = nlsolve(vL!,
                   zguess .* scale
                   ;method=:trust_region,
                    ftol=1e-12,
                    xtol=1e-12
           )

    if !converged(root)
        throw(ConvergenceError(":trust_region failed to converge on a root", root))
    end
    root.zero ./ scale
end


end # module RootFinding
