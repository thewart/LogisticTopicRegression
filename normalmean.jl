mutable struct NormalMeanPosterior <: PostPredSS
  sumy::Float64
  n::Float64
  σ2::Float64
  NormalMeanPosterior(y::Float64=0.0,n::Float64=0.001,σ2::Float64=1.0) = new(y,n,σ2)
end
NormalMeanPosterior(y::Vector{Float64},n0=0.001::Float64,σ2::Float64=1.0) =
  NormalMeanPosterior(sum(y),n0+length(y),σ2);

dim(pp::NormalMeanPosterior) = 1;

function addsample!(pp::NormalMeanPosterior,ynew::Float64)
  pp.n += 1;
  pp.sumy += ynew;
end

function pullsample!(pp::NormalMeanPosterior,yold::Float64)
  pp.n -= 1;
  pp.sumy -= yold;
end

function lppd(y::Float64,pp::NormalMeanPosterior)
    σ2 = pp.σ2/pp.n + pp.σ2;
    μ = pp.sumy/pp.n;
  return -0.5*log(2π*σ2) - 0.5*(y-μ)^2/σ2
end

function topicpd(pp::NormalMeanPosterior)
  σ2 = pp.σ2/pp.n;
  μ = pp.sumy/pp.n;
  return Normal(μ,sqrt(σ2))
end

function topicppd(pp::NormalMeanPosterior)
    σ2 = pp.σ2/pp.n + pp.σ2;
    μ = pp.sumy/pp.n;
    return Normal(μ,sqrt(σ2))
end

randtopic(pp::NormalMeanPosterior) = Normal(rand(topicpd(pp)),sqrt(pp.σ2));
