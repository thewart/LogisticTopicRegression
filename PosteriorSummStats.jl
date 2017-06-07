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

function addsample!{T<:Real}(pp::VectorPosterior,ynew::Vector{T})
  for i in 1:length(pp) addsample!(pp[i],ynew[pp.span[i]]); end
end

function pullsample!{T<:Real}(pp::VectorPosterior,yold::Vector{T})
  for i in 1:length(pp) pullsample!(pp[i],yold[pp.span[i]]); end
end

#this is a hack! should convert to taking abstract arrays
lppd{T<:Real}(pp::VectorPosterior,y::T) = lppd(pp,[y]);
addsample!{T<:Real}(pp::VectorPosterior,ynew::T) = addsample!(pp,[ynew]);
pullsample!{T<:Real}(pp::VectorPosterior,yold::T) = pullsample!(pp,[yold]);
