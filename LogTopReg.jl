module LogTopReg

using Distributions
using StatsBase
using StatsFuns
using PolyaGamma
import Base.rand
import Base.mean

abstract type PostPredSS end

export
  #types
  TLMMsample,
  TLMMfit,
  PostPredSS,
  PoissonPosterior,
  DirMultPosterior,
  CategoricalPosterior,
  NormalMeanPosterior,

  #functions
  gf,
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
  rand,
  sim_dat,
  randtopic,
  init_params,
  init_params!,
  init_z,

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
