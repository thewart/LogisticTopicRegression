abstract PostPredSS
import Base.length
import Base.getindex
import Base.rand

type VectorPosterior{T<:PostPredSS}
  PP::Vector{T}
  span::Vector{Union{Int64,UnitRange{Int64}}}
end

function VectorPosterior{T<:PostPredSS}(pp::Vector{T})
  span = Vector{Union{Int64,UnitRange{Int64}}}(length(pp));
  j = 0;
  for i in 1:length(pp)
    dpp = dim(pp[i]);
    dpp > 1 ? span[i] = ((j+1):(j+dpp)) : (span[i] = j+1);
    j += dim(pp[i]);
  end
  VectorPosterior(pp,span)
end
function VectorPosterior(x...)
  l = cumsum(collect(x[2:2:length(x)]));
  pp = Vector{PostPredSS}(maximum(l));
  for i in 1:maximum(l)
    pp[i] = deepcopy(x[findfirst(l .>= i)*2 - 1]);
  end
  VectorPosterior(pp)
end

getindex(VP::VectorPosterior,i) = VP.PP[i];
length(VP::VectorPosterior) = length(VP.PP);

#log posterior predictive density
function lppd{T<:Real}(pp::VectorPosterior,y::AbstractVector{T})
  out = 0.0;
  for i in 1:length(pp)
    out += lppd(pp[i],y[pp.span[i]]);
  end
  return out
end

function lppd{T<:Real}(pp::VectorPosterior,y::Array{T,2})
  out = 0.0;
  for j in 1:size(y)[2]
    out += lppd(pp,y[:,j]);
  end
  return out
end

#lppd{T<:Real}(pp::VectorPosterior,y::T) =
#  for i in 1:length(pp) lppd(pp[i],y); end

#log posterior predictive, integrating across indicator variables
function lppd{T<:Real}(ppv::Vector{VectorPosterior},π::Vector{Float64},y::AbstractVector{T})
  lp = Float64[log(π[i]) + lppd(ppv[i],y) for i in 1:length(ppv)];
  return logsumexp(lp)
end

lppd{T<:Real}(ppv::Vector{VectorPosterior},π::Vector{Float64},y::Array{T,2}) =
  Float64[lppd(ppv,π,y[:,i]) for i in 1:size(y)[2]]

function lppd{T<:Real}(fit::Dict{Symbol,AbstractArray},y::Vector{Array{T,2}},pointwise::Bool=true)
  n = length(y);
  m = length(fit[:τ]);
  nd = map(i -> size(y[i])[2],1:n);
  lp = Array{Array{Float64,2}}(n);
  for i in 1:n lp[i] = Array{Float64}(m,nd[i]); end

  for i in 1:n
    for j in 1:m
      π = softmax(fit[:η][:,i,j]);
      lp[i][j,:] = lppd(fit[:topic][j],π,y[i]);
    end

    if !pointwise lp[i] = sum(lp[i],1); end
  end

  if !pointwise lp = vcat(lp...); end
  return lp
end


#add one new observation
function addsample!{T<:Real}(pp::VectorPosterior,ynew::AbstractVector{T})
  for i in 1:length(pp) addsample!(pp[i],ynew[pp.span[i]]); end
end

addsample!{T<:Real}(pp::VectorPosterior,ynew::Array{T,2}) =
  for i in 1:length(pp) addsample!(pp,ynew[:,j]); end

addsample!{T<:Real}(pp::VectorPosterior,ynew::T) =
  for i in 1:length(pp) addsample!(pp[i],ynew); end

#remove one observation
function pullsample!{T<:Real}(pp::VectorPosterior,yold::AbstractVector{T})
  for i in 1:length(pp) pullsample!(pp[i],yold[pp.span[i]]); end
end

pullsample!{T<:Real}(pp::VectorPosterior,yold::Array{T,2}) =
  for i in 1:length(pp) pullsample!(pp,yold[:,j]); end

pullsample!{T<:Real}(pp::VectorPosterior,yold::T) =
  for i in 1:length(pp) pullsample!(pp[i],yold); end


#get topic parameter distribution
function topicpd(pp::VectorPosterior)
  topic = Vector{Sampleable}(length(pp));
  for i in 1:length(pp) topic[i] = topicpd(pp[i]); end
  return topic
end

#utilities for distribution vectors
function rand{T<:Sampleable}(dv::Vector{T},n::Int64=1)
  l = map(length,dv);
  p = sum(l);
  X = Array{Float64}(p,n);
  for i in 1:n
    k = 0;
    for j in 1:p
      rand!(dv[j],sub(X,(k+1):(k+l[j]),i));
      k = k + l[j];
    end
  end
  return X
end

include("/home/seth/code/LogisticTopicRegression/gammapoisson.jl")
include("/home/seth/code/LogisticTopicRegression/dirichletmultinomial.jl")
include("/home/seth/code/polyagamma/polygamma.jl")
