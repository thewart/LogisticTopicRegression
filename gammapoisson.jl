type PoissonPosterior <: PostPredSS
  n::Float64
  ys::Float64
end
PoissonPosterior(n::Real,ys::Real) = PoissonPosterior(Float64(n),Float64(ys))
PoissonPosterior(y::Vector{Int},a::Float64,b::Float64) = PoissonPosterior(length(y)+b, sum(y)+a)

dim(pp::PoissonPosterior) = 1;

function addsample!(pp::PoissonPosterior,ynew::Int)
  pp.n += 1;
  pp.ys += ynew;
end
addsample!{T<:Real}(pp::PoissonPosterior,ynew::T) = addsample!(pp,round(Int,ynew));

function pullsample!(pp::PoissonPosterior,yold::Int)
  pp.n -= 1;
  pp.ys -= yold;
end
pullsample!{T<:Real}(pp::PoissonPosterior,yold::T) = pullsample!(pp,round(Int,yold));

function lppd(y::Int,pp::PoissonPosterior)
  p = inv(pp.n+1);
  r = pp.ys;
  return lgamma(r+y) - lgamma(r) + r*log(1-p) + y*log(p) - lgamma(y+1);
end
lppd{T<:Real}(y::T,pp::PoissonPosterior) = lppd(round(Int,y),pp);

function topicpd(pp::PoissonPosterior)
  return Gamma(pp.ys,inv(pp.n))
end

function topicppd(pp::PoissonPosterior)
  return NegativeBinomial(pp.ys,1-inv(pp.n+1))
end
