module LogTopReg

abstract PostPredSS
using Distributions
using StatsBase
using StatsFuns
using PolyaGamma
import Base.length
import Base.getindex
import Base.rand
import Base.vcat
import Base.mean

export
  #types
  PostPredSS,
  VectorPosterior,
  PoissonPosterior,
  DirMultPosterior,

  #functions
  refβ,
  lppd,
  addsample!,
  pullsample!,
  topicpd,
  softmax,
  vcat,
  getindex,
  rand,
  length,

  #the whole point
  topiclmm

include("vectorposterior.jl")
include("gammapoisson.jl");
include("dirichletmultinomial.jl");
include("topicllm.jl");

end
