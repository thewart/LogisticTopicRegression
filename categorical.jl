type CategoricalPosterior <: PostPredSS
  K::Int64
  n::Vector{Float64}

  CategoricalPosterior(y::Vector{Float64}) = new(length(y),y)
end
CategoricalPosterior(y::Vector{Int},α::Vector{Float64}) = CategoricalPosterior(y.+α);
CategoricalPosterior(K::Int64,α::Float64=1.0) = CategoricalPosterior(fill(α,K));

dim(pp::CategoricalPosterior) = 1;

function addsample!(pp::CategoricalPosterior,ynew::Int64) pp.n[ynew] += 1.0; end

function pullsample!(pp::CategoricalPosterior,yold::Int64) pp.n[yold] -= 1.0; end

function lppd(y::Int64,pp::CategoricalPosterior)
  N = sum(pp.n);
  return lgamma(N) - lgamma(N+1) + lgamma(pp.n[y]+1.0) - lgamma(pp.n[y])
end
lppd{T<:Real}(y::Vector{T},pp::CategoricalPosterior) = lppd(pp,round(Int64,y));

function topicpd(pp::CategoricalPosterior)
  return Dirichlet(pp.n)
end
function topicppd(pp::CategoricalPosterior)
  return Categorical(pp.n./sum(pp.n))
end
