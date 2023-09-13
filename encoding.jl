using Interpolations
using Statistics

struct REALVALUE_ENCODER
    D
    l
    signs
    shifts
    bipolar
    normalization
end

function realvalue_encoder(l, D; bipolar=false, normalization=true)
    random_signs = rand((-1,1),(Int64(ceil(1/l)+1),D))
    random_shifts = l*rand(D)

    x = range(0,1,length=200)
    ϕx = encode(REALVALUE_ENCODER(D, l, random_signs, random_shifts, bipolar, "to be determined"), x, normalized=false)
    nx = find_normalization(ϕx*ϕx' / D;n_iter=10)
    
    return REALVALUE_ENCODER(D, l, random_signs, random_shifts, bipolar, LinearInterpolation(x,nx))
end

function encode(encoder::REALVALUE_ENCODER, x; normalized=true)
    res = zeros(length(x),encoder.D)
    for i in 1:encoder.D
        res[:,i] .= @. (1 .- cos.(2π/encoder.l*(x.+encoder.shifts[i]))) .* encoder.signs[Int64.(ceil.((x.+encoder.shifts[i])/encoder.l)), i] /2 #.* boxsign[box]
    end
    if encoder.bipolar==true
        res .= sign.(res)
    end
    if normalized == true
        return res ./ encoder.normalization.(x)
    else
        return res
    end
end

function find_normalization(K;n_iter=10)
    n = sqrt.(mean(K, dims=1))[:]
    onetilde = mean(K ./ n, dims=1)[:] ./ n
    for i in 1:n_iter
        n[:] = n .* sqrt.(onetilde)
        onetilde[:] = mean(K ./ n, dims=1)[:] ./ n
    end
    return n
end