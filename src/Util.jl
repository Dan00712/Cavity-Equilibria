module Util

using Dates

using IntervalArithmetic

using ..Parameters

export roottomidpoint, now_nodots, prepare_freq

function now_nodots()
    now() |> string |> x->split(x, ".")[1] |> x -> replace(x, ":" => "-")
end

function findclosest(a, binit, c)
"""
	We have that a * b = n * c
	binit does not satisfy this,
	but we want the closest be to it
"""
	if a == 0
		return binit
	end
	n = round(a*binit/c)
	b = n*c/a
	b
end

function prepare_freq(Δ, d; params::SystemParams)
	if d == 0
		return Δ
	end
	ωc = Δ + params.ω0
	kc = let
		ωc/params.c
	end
	kc = findclosest(d/2, kc, π)
	kc *params.c
end

function roottomidpoint(root)
    I = root.interval

    vz = mid.(I)
    err = radius.(I)

    vz, err
end

end # module Util
