# Gumbel Kernel
function gumbel(β, a, b)
    z = (β - b) / a
    return (1.0 / a) * exp(-z - exp(-z))
end