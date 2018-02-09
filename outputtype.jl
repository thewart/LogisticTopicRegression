type TLMMsample
    z::Vector{Int64}
    τ_μ::Float64
    σ2_η::Vector{Float64}
    η::Matrix{Float64}

    μ::Vector{Float64}
    β::Matrix{Float64}

    u::Matrix{Float64}
    τ_u::Vector{Float64}

    TLMMsample() = new()
end

function init_params!(samp::TLMMsample,K,n,p)
    samp.η = rand(K,n)*2;
    samp.τ_μ = exp(randn())*2;
    samp.σ2_η = rand(K)*2;
    samp.τ_u = exp.(randn(K));

    samp.μ = Vector{Float64}(K);
    samp.u = Matrix{Float64}(n,K);
    samp.β = Matrix{Float64}(p,K);
    return samp
end

function init_params(K,n,p)
    return init_params!(TLMMsample(),K,n,p);
end

type TLMMfit{T<:PostPredSS}
    θ::Vector{TLMMsample}
    tss::Matrix{VectorPosterior{T}}
    prior::VectorPosterior{T}
    hyperparameter::Dict{Symbol,Float64}
end
