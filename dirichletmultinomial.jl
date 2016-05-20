type DirMultPosterior <: PostPredSS
  K::Int64
  n::Vector{Float64}

  DirMultPosterior(y::Vector{Float64}) = new(length(y),y)
end

DirMultPosterior(y::Vector{Int},α::Vector{Float64}) = DirMultPosterior(y.+α);
DirMultPosterior(K::Int64,α::Float64) = DirMultPosterior(fill(α,K));

dim(pp::DirMultPosterior) = pp.K;

function addsample!(pp::DirMultPosterior,y::Vector{Int})
  for k in 1:pp.K
    pp.n[k] += y[k];
  end
end

function pullsample!(pp::DirMultPosterior,y::Vector{Int})
  for k in 1:pp.K
    pp.n[k] -= y[k];
  end
end

function lppd{T<:Real}(pp::DirMultPosterior,y::Vector{T})
  Nprior = sum(pp.n);
  Nnew = sum(y);
  d = lgamma(Nprior) - lgamma(Nprior+Nnew) + lfact(Nnew);
  for k in 1:K
    d += lgamma(pp.n[k] + y[k]) - lgamma(pp.n[k]) - lfact(y[k]);
  end
  return d
end
