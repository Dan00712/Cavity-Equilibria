module Parameters

export SystemParams, DEFAULT_PARAMS

@kwdef struct SystemParams
    # using SI Units
    ħ = 6.626e-34/2π       # Js
    ϵ0 = 8.854e-12      # F/m
    c = 2.991e+8            # m/s


    # Tweezer
    λ0 = 1550e-9        # m = 1.66μm
    k0 = 2π/λ0
    ω0 = c*k0

    Pt = 0.5            # W
    Wt = 1e-6           # m

    Ax = 1
    Ay = 1

    #	zR = 1/2 * k0 * Wt^2
    zR = π * Wt^2/λ0


    E0 = sqrt(4*Pt / π / ϵ0 / c / Wt^2 / Ax^2 / Ay)

    # SiO₂
    ρ = 2200              # kg/m^3
    ϵ = 2.07
    R = 100e-9            # m
    V = 4/3 * π * R^3
    α = 3 * V * ϵ0 * (ϵ-1)/(ϵ+2)


    # Cavity
    Wc = 20e-6        # m
    Lc = 19.8e-3      # m

    # keep ω_0 = const
    Vc = π*Wc^2*Lc/4
    m = ρ * (4/3 * π * R^3)

    ζ = 1
end
const DEFAULT_PARAMS = SystemParams()

function Base.show(io::IO, p::SystemParams)
    T = typeof(p)
    fields = fieldnames(T)

    print(io, nameof(T), "(")
    for (i, f) in enumerate(fields)
        val = getfield(p, f)
        print(io, f, "=", val)
        i < length(fields) && print(io, ", ")
    end
    print(io, ")")
end

end # module Parameters 
