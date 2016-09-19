#log posterior predictive density
function lppd{T<:Real}(y::AbstractVector{T},pp::VectorPosterior)
  out = 0.0;
  for i in 1:length(pp)
    out += lppd(y[pp.span[i]],pp[i]);
  end
  return out
end

function lppd{T<:Real}(y::Array{T,2},pp::VectorPosterior)
  out = 0.0;
  for j in 1:size(y)[2]
    out += lppd(pp,y[:,j]);
  end
  return out
end

#log posterior predictive, integrating across indicator variables
function lppd{T<:Real}(y::AbstractVector{T},ppv::Vector{VectorPosterior},π::Vector{Float64})
  K = length(ppv);
  lp = Vector{Float64}(K);
  for k in 1:K lp[k] = log(π[k]) + lppd(y,ppv[k]); end
  return logsumexp(lp)
end

function lppd{T<:Real}(y::Array{T,2},ppv::Vector{VectorPosterior},π::Vector{Float64})
  n = size(y)[2];
  lp = 0;
  for i in 1:n lp += lppd(y[:,i],ppv,π); end
  return lp
end

#predictive density for new observations in same groups
function lppd{T<:Real}(y::Vector{Array{T,2}},fit::Dict{Symbol,Union{AbstractArray,Dict{Symbol,Float64}}})

  n = length(y);
  m = size(η)[3];
  lp = Array{Float64,2}(n,m);
  for i in 1:n
    for j in 1:m
      π = softmax(fit[:η][:,i,j]);
      lp[i,j] = lppd(y[i],fit[:topic][:,j],π);
    end
  end
  return lp

end

#predictive density for new groups
function lppd{T<:Real}(y::Array{T,2},X::Vector{Float64},topic::Vector{VectorPosterior},
  β::Array{Float64,2},μ::Vector{Float64},σ2::Vector{Float64},nsim::Int64=1)

  π = Vector{Float64}(length(topic));
  ηhat = β'*X;
  broadcast!(+,ηhat,ηhat,μ)
  lp = 0;
  p = length(X);

  for i in 1:nsim
    for j in 1:length(π) π[j] = ηhat[j] + sqrt(σ2[j])*randn(); end
    softmax!(π);
    lp += lppd(y,topic,π);
  end

  return lp/nsim
end

function lppd{T<:Real}(y::Vector{Array{T,2}},X::Array{Float64,2},
  fit::Dict{Symbol,Union{AbstractArray,Dict{Symbol,Float64}}},nsim::Int64=1)

  n = length(y);
  m = size(fit[:β])[3];
  lp = Array{Float64,2}(n,m);
  for i in 1:n
    for j in 1:m
      lp[i,j] = lppd(y[i],X[:,i],fit[:topic][:,j],fit[:β][:,:,j],fit[:μ][:,j],fit[:σ2][:,j],nsim);
    end
  end
  return lp
end
