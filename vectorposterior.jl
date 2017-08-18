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
function VectorPosterior(x::Union{PostPredSS,Int64}...)
  l = cumsum(collect(x[2:2:length(x)]));
  pp = Vector{PostPredSS}(maximum(l));
  for i in 1:maximum(l)
    pp[i] = deepcopy(x[findfirst(l .>= i)*2 - 1]);
  end
  VectorPosterior(pp)
end
VectorPosterior(pp::PostPredSS) = VectorPosterior([pp]);

getindex(VP::VectorPosterior,i) = VP.PP[i];
length(VP::VectorPosterior) = length(VP.PP);
length(PP::PostPredSS) = 1;
function vcat{T<:PostPredSS}(VP::Union{VectorPosterior,T}...)
  l = map(length,VP);
  cat = Vector{PostPredSS}(sum(l));
  j = 0;
  for i in 1:length(VP)
    cat[(j+1):(j+l[i])] = (typeof(VP[i]) <: VectorPosterior) ? VP[i].PP : VP[i];
    j += l[i];
  end
  return VectorPosterior(cat)
end


#add one new observation
function addsample!{T<:Real}(pp::VectorPosterior,ynew::AbstractVector{T})
  for i in 1:length(pp) addsample!(pp[i],ynew[pp.span[i]]); end
end

addsample!{T<:Real}(pp::VectorPosterior,ynew::Array{T,2}) =
  for i in 1:length(pp) addsample!(pp,ynew[:,j]); end

#addsample!{T<:Real}(pp::VectorPosterior,ynew::T) =
#  for i in 1:length(pp) addsample!(pp[i],ynew); end

#remove one observation
function pullsample!{T<:Real}(pp::VectorPosterior,yold::AbstractVector{T})
  for i in 1:length(pp) pullsample!(pp[i],yold[pp.span[i]]); end
end

pullsample!{T<:Real}(pp::VectorPosterior,yold::Array{T,2}) =
  for i in 1:length(pp) pullsample!(pp,yold[:,j]); end

#pullsample!{T<:Real}(pp::VectorPosterior,yold::T) =
#  for i in 1:length(pp) pullsample!(pp[i],yold); end


#get topic parameter distribution
function topicpd(pp::VectorPosterior)
  topic = Vector{Sampleable}(length(pp));
  for i in 1:length(pp) topic[i] = topicpd(pp[i]); end
  return topic
end

function topicppd(pp::VectorPosterior)
  topic = Vector{Sampleable}(length(pp));
  for i in 1:length(pp) topic[i] = topicppd(pp[i]); end
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

function mean{T<:Sampleable}(dv::Vector{T})
  l = map(length,dv);
  p = sum(l);
  X = Vector{Float64}(p);
  k = 0;
  for j in 1:p
    X[(k+1):(k+l[j])] = mean(dv[j]);
    k = k + l[j];
  end
  return X
end
