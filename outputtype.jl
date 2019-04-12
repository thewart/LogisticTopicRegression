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

function init_params!(samp::TLMMsample,K::Int,n::Int,p::Int;
    ranef=false,collapsed=false)

    samp.η = rand(K,n)*2;
    samp.σ2_η = rand(K)*2;

    samp.μ = Vector{Float64}(K);
    samp.β = Matrix{Float64}(p,K);
    if ranef
        samp.u = Matrix{Float64}(n,K);
        samp.τ_u = exp.(randn(K));
    end
    if collapsed
        samp.τ_μ = exp(randn())*2;
    end

    return samp
end

function init_params(K::Int,n::Int,p::Int)
    return init_params!(TLMMsample(),K,n,p);
end

function mean(samps::Vector{TLMMsample})
    meanfit = TLMMsample();
    tmp = gf(samps,:τ_μ);
    meanfit.τ_μ = mean(tmp);
    tmp = gf(samps,:σ2_η);
    meanfit.σ2_η = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));
    tmp = gf(samps,:η);
    meanfit.η = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));
    tmp = gf(samps,:μ);
    meanfit.μ = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));
    tmp = gf(samps,:β);
    meanfit.β = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));

    if isdefined(samps[1],:u)
        tmp = gf(samps,:u);
        meanfit.u = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));
        tmp = gf(samps,:τ_u);
        meanfit.τ_u = squeeze(mean(tmp,ndims(tmp)),ndims(tmp));
    end

    return meanfit
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
