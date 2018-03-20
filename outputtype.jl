mutable struct TLMMsample
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

function init_params!(samp::TLMMsample,K::Int,n::Int,p::Int;ranef=false)
    samp.η = rand(K,n)*2;
    samp.τ_μ = exp(randn())*2;
    samp.σ2_η = rand(K)*2;

    samp.μ = Vector{Float64}(K);
    samp.β = Matrix{Float64}(p,K);
    if ranef
        samp.u = Matrix{Float64}(n,K);
        samp.τ_u = exp.(randn(K));
    end
    return samp
end

function init_params(K::Int,n::Int,p::Int)
    return init_params!(TLMMsample(),K,n,p);
end

struct TLMMfit{U<:PostPredSS}
    θ::Vector{TLMMsample}
    tss::Matrix{Vector{U}}
    prior::Vector{U}
    hyperparameter::Dict{Symbol,Float64}
end

function gf(value::Vector{TLMMsample},name::Symbol)
    nd = ndims(getfield(value[1],name));
    return cat(nd+1,getfield.(value,name)...)
end
