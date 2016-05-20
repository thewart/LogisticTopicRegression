type PoissonPosterior <: PostPredSS
  n::Float64
  ys::Float64
end
PoissonPosterior(n::Real,ys::Real) = PoissonPosterior(float64(n),float64(ys))
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

function lppd(pp::PoissonPosterior,y::Int)
  p,r = ppdparam(pp);
  return lgamma(r+y) - lgamma(r) + r*log(1-p) + y*log(p) - lgamma(y+1);
end
lppd{T<:Real}(pp::PoissonPosterior,y::T) = lppd(pp,round(Int,y));

function ppdparam(pp::PoissonPosterior)
  p = inv(pp.n+1);
  r = pp.ys;
  return p,r
end

function topicpd(pp::PoissonPosterior)
  return Gamma(pp.ys,inv(pp.n))
end
