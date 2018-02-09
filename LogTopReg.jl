module LogTopReg

using Distributions
using StatsBase
using StatsFuns
using PolyaGamma
using RCall
import Base.length
import Base.getindex
import Base.rand
import Base.vcat
import Base.mean

abstract type PostPredSS end

export
  #types
  TLMMsample,
  TLMMfit,
  PostPredSS,
  VectorPosterior,
  PoissonPosterior,
  DirMultPosterior,
  CategoricalPosterior,
  NormalMeanPosterior,

  #functions
  simplemix,
  hyperparameter,
  refβ,
  writefit,
  lppd,
  addsample!,
  pullsample!,
  topicpd,
  topicppd,
  softmax,
  vcat,
  getindex,
  rand,
  length,
  sim_dat,
  mlinit_dirichlet,
  randtopic,
  init_params,
  init_params!,

  #the whole point
  topiclmm

include("vectorposterior.jl")
include("outputtype.jl")
include("gammapoisson.jl")
include("dirichletmultinomial.jl")
include("categorical.jl")
include("normalmean.jl")
include("topicllm.jl")
include("topicllm_uncolapsed.jl")
include("prediction.jl")
include("llmhelpers.jl")
include("llmsamplers.jl")

end
