type DirMultPosterior <: PostPredSS
  K::Int64
  n::Vector{Float64}

  DirMultPosterior(y::Vector{Float64}) = new(length(y),y)
end
DirMultPosterior(y::Vector{Int},α::Vector{Float64}) = DirMultPosterior(y.+α);
DirMultPosterior(K::Int64,α::Float64=1.0) = DirMultPosterior(fill(α,K));

dim(pp::DirMultPosterior) = pp.K;

function addsample!(pp::DirMultPosterior,ynew::Vector{Int})
  for k in 1:pp.K
    pp.n[k] += ynew[k];
  end
end
addsample!{T<:Real}(pp::DirMultPosterior,ynew::Vector{T}) =
    addsample!(pp,round(Int64,ynew));

function pullsample!(pp::DirMultPosterior,yold::Vector{Int})
  for k in 1:pp.K
    pp.n[k] -= yold[k];
  end
end
pullsample!{T<:Real}(pp::DirMultPosterior,yold::Vector{T}) =
    pullsample!(pp,round(Int64,yold));

function lppd(y::Vector{Int},pp::DirMultPosterior)
  Nprior = sum(pp.n);
  Nnew = sum(y);
  d = lgamma(Nprior) - lgamma(Nprior+Nnew) + lfact(Nnew);
  for k in 1:pp.K
    d += lgamma(pp.n[k] + y[k]) - lgamma(pp.n[k]) - lfact(y[k]);
  end
  return d
end
lppd{T<:Real}(y::Vector{T},pp::DirMultPosterior) = lppd(pp,round(Int64,y));

function topicpd(pp::DirMultPosterior)
  return Dirichlet(pp.n)
end

randtopic(pp::DirMultPosterior,n::Int64=1) = Multinomial(n,rand(topicpd(pp)));
