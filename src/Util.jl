module Util

using IntervalArithmetic

export roottomidpoint

function roottomidpoint(root)
    I = root.interval

    vz = mid.(I)
    err = radius.(I)

    vz, err
end

end # module Util
