module Util

using Dates

export now_nodots

function now_nodots()
    now() |> string |> x->split(x, ".")[1] |> x -> replace(x, ":" => "-")
end

end # module Util
