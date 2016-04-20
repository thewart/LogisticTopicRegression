abstract PostPredSS
import Base.length
import Base.getindex

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

getindex(VP::VectorPosterior,i) = VP.PP[i];
length(VP::VectorPosterior) = length(VP.PP);

#log posterior predictive density
function lppd{T<:AbstractVector{Real}}(pp::VectorPosterior,y::T)
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

lppd{T<:Real}(pp::VectorPosterior,y::T) =
  for i in 1:length(pp) lppd(pp[i],y); end

#add one new observation
function addsample!{T<:AbstractVector{Real}}(pp::VectorPosterior,ynew::T)
  for i in 1:length(pp) addsample!(pp[i],ynew[pp.span[i]]); end
end

addsample!{T<:Real}(pp::VectorPosterior,ynew::Array{T,2}) =
  for i in 1:length(pp) addsample!(pp,ynew[:,j]); end

addsample!{T<:Real}(pp::VectorPosterior,ynew::T) =
  for i in 1:length(pp) addsample!(pp[i],ynew); end

#remove one observation
function pullsample!{T<:AbstractVector{Real}}(pp::VectorPosterior,yold::T)
  for i in 1:length(pp) pullsample!(pp[i],yold[pp.span[i]]); end
end

pullsample!{T<:Real}(pp::VectorPosterior,yold::Array{T,2}) =
  for i in 1:length(pp) pullsample!(pp,yold[:,j]); end

pullsample!{T<:Real}(pp::VectorPosterior,yold::T) =
  for i in 1:length(pp) pullsample!(pp[i],[yold]); end


include("/home/seth/code/LogisticTopicRegression/gammapoisson.jl")
include("/home/seth/code/polyagamma/polygamma.jl")
