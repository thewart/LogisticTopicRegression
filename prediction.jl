#log posterior predictive density
function lppd{T<:Real,U<:PostPredSS}(y::AbstractVector{T},pp::Vector{U})
  out = 0.0;
  for i in 1:length(pp)
    out += lppd(y[i],pp[i]);
  end
  return out
end

function lppd{T<:Real,U<:PostPredSS}(y::Array{T,2},pp::Vector{U})
  out = 0.0;
  for j in 1:size(y)[2]
    out += lppd(y[:,j],pp);
  end
  return out
end

#log posterior predictive, integrating across indicator variables
function lppd{T<:Real,U<:PostPredSS}(y::AbstractVector{T},ppv::Vector{Vector{U}},π::Vector{Float64})
  K = length(ppv);
  lp = Vector{Float64}(K);
  for k in 1:K lp[k] = log(π[k]) + lppd(y,ppv[k]); end
  return logsumexp(lp)
end

function lppd{T<:Real,U<:PostPredSS}(y::Array{T,2},ppv::Vector{Vector{U}},π::Vector{Float64})
  n = size(y)[2];
  lp = 0;
  for i in 1:n lp += lppd(y[:,i],ppv,π); end
  return lp
end

#predictive density for new observations in same groups
function lppd{U<:PostPredSS}(y,docrng,fit::TLMMfit{U})

  n = length(docrng);
  m = length(fit.θ);
  lp = Array{Float64,2}(n,m);
  for i in 1:m
    for j in 1:n
        π = softmax(fit.θ[i][:,j]);
        lp[i,j] = lppd(y[:,docrng[i]],fit.tss[:,i],π);
    end
  end
  return lp

end

#predictive density for new groups
function lppd{T<:Real,U<:PostPredSS}(y::Array{T,2},X::Vector{Float64},topic::Vector{Vector{U}},
  β::Array{Float64,2},μ::Vector{Float64},σ2::Vector{Float64},nsim::Int64=1)

  π = Vector{Float64}(lngth(topic));
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

function lppd{T<:Real}(y::Vector{Array{T,2}},X::Array{Float64,2},fit,nsim::Int64=1)

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
