type NormalMeanPosterior <: PostPredSS
  sumy::Float64
  ss::Float64
  invσ2::Float64
  NormalMeanPosterior(y::Float64,ss::Float64,invσ2::Float64) = new(y*invσ2,ss,invσ2)
end
NormalMeanPosterior(y::Vector{Float64},ss0=0.001,invσ2::Float64=1.0) =
  NormalMeanPosterior(sum(y)*invσ2,ss0+length(y)*invσ2,invσ2);

dim(pp::NormalMeanPosterior) = 1;

function addsample!(pp::NormalMeanPosterior,ynew::Float64)
  pp.ss += pp.invσ2;
  pp.sumy += ynew*pp.invσ2;
end

function pullsample!(pp::NormalMeanPosterior,yold::Float64)
  pp.ss -= pp.invσ2;
  pp.sumy -= yold*pp.invσ2;
end

function lppd(y::Float64,pp::NormalMeanPosterior)
  σ2 = inv(pp.ss) + pp.invσ2;
  μ = pp.sumy*inv(pp.ss)*pp.invσ2;
  return -0.5*log(2π*σ2) - 0.5*(y-μ)^2/σ2
end

function topicpd(pp::NormalMeanPosterior)
  σ2 = inv(pp.ss);
  μ = pp.sumy*inv(pp.ss)*pp.invσ2;
  return Normal(μ,sqrt(σ2))
end

function topicppd(pp::NormalMeanPosterior)
  σ2 = inv(pp.ss) + pp.invσ2;
  μ = pp.sumy*inv(pp.ss)*pp.invσ2;
  return Normal(μ,sqrt(σ2))
end
