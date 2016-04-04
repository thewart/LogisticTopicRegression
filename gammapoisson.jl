#Gamma-poisson
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
addsample!{T<:Real}(pp::PoissonPosterior,ynew::T) = addsample!(pp,int(ynew));

function pullsample!(pp::PoissonPosterior,yold::Int)
  pp.n -= 1;
  pp.ys -= yold;
end
pullsample!{T<:Real}(pp::PoissonPosterior,yold::T) = addsample!(pp,int(yold));

function lppd(pp::PoissonPosterior,y::Int)
  p = inv(pp.n+1);
  r = pp.ys;
  return lgamma(r+y) - lgamma(r) + r*log(1-p) + y*log(p) # - lgamma(y+1) ;
end

