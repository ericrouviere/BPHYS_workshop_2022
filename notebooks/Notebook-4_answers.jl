

function buildSpringTable(q, k_min, k_max)
    log_diff = log(10, k_max) - log(10, k_min)
    K =  10 .^ (rand(q, q) .* log_diff .+ log(10,k_min))
    K = Float64.(Symmetric(K, :L))
    return K
end


------------------------------------------------------