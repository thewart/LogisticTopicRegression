function makeprec(XtX::Array{Float64,2},ZtZ::Array{Float64,2},
  τ_β::Float64,τ_u::Float64,τ_μ::Float64)

  n = size(XtX)[1];
  V = Array{Float64,2}(n,n);
  for i in 1:n
    for j in 1:n
      V[j,i] = XtX[j,i].*τ_β + ZtZ[j,i].*τ_u + (i==j ? 1+τ_μ : τ_μ);
    end
  end

  return inv(V)
end

function makeprec(XtX::Array{Float64,2},τ_β::Float64,τ_μ::Float64)

  n = size(XtX)[1];
  V = Array{Float64,2}(n,n);
  for i in 1:n
    for j in 1:n
      V[j,i] = XtX[j,i].*τ_β + (i==j ? 1+τ_μ : τ_μ);
    end
  end

  return inv(V)
end

function hyperparameter(;ν0_σ2η=0.01,σ0_σ2η = 0.01,
  τ0_τ = 0.01,ν0_τ = 0.01,τ0_u = 0.01,ν0_u=0.01,τ_β=5.0)
  return Dict(:ν0_σ2η => ν0_σ2η, :σ0_σ2η => σ0_σ2η, :τ0_τ => τ0_τ,
  :ν0_τ => ν0_τ, :τ0_u => τ0_u, :ν0_u => ν0_u, :τ_β => τ_β)
end

function init_topic(pss0,K,z,y)
  topic = Vector{typeof(pss0)}(K);
  map!(k -> deepcopy(pss0),topic,1:K);
  for i in 1:length(z) addsample!(topic[z[i]],y[:,i]); end
  return topic
end

function init_z(η,K,docrng)
  n = length(docrng);
  nobs = maximum(last.(docrng));
  z = Vector{Int64}(nobs);
  for i in 1:n z[docrng[i]] = sample(1:K,Weights(softmax(η[:,i])),length(docrng[i])); end
  return z
end

function refβ(β::Array{Float64,2},refk::Int64)
  return β .- β[:,refk:refk]
end
refβ(β::Array{Float64,3},refk::Int64) = β .- β[:,refk:refk,:]
refβ(β::Array{Float64},μ::Array{Float64,1}) = refβ(β,findmax(μ)[2]);
refβ(β::Array{Float64},μ::Array{Float64,2}) = refβ(β,findmax(mean(μ,2))[2]);

function countz(z,docrng,K)
    n = length(docrng);
    nk = fill(0,(K,n));
    for i in 1:n
        for k in 1:K
            nk[k,i] = countnz(z[docrng[i]].==k)
        end
    end
    return nk
end

function writefit(fit::Dict{Symbol,Union{AbstractArray,Dict{Symbol,Float64}}},path::String)
  if !isdir(path) mkpath(path); end
  n = length(fit[:z]);
  K,nsave = size(fit[:topic]);
  p = size(fit[:β])[1];
    b = length(fit[:topic][1]);
  writecsv(string(path,"dims.csv"),[n,K,p,b,nsave]);
  writecsv(string(path,"sigma.csv"),hcat(fit[:τ],fit[:σ2]',fit[:τ_u]'));
  writecsv(string(path,"beta.csv"),fit[:β][:]);
  writecsv(string(path,"eta.csv"),fit[:η][:]);
  writecsv(string(path,"u.csv"),fit[:u][:]);
  #writecsv(string(path,"loglik.csv"),fit[:loglik]');
  writecsv(string(path,"mu.csv"),fit[:μ]');
  writecsv(string(path,"z.csv"),vcat(fit[:z]...));
  writecsv(string(path,"nd.csv"),map(x -> size(x)[1],fit[:z]));
  writecsv(string(path,"topicmean.csv"),
            mapreduce(y -> mean(topicppd(y)),vcat,fit[:topic]));
  writecsv(string(path,"topicvar.csv"),
            mapreduce(x -> map(var,topicppd(x)),vcat,fit[:topic]));

  θ = Vector{Float64}(0);
  for i in 1:(length(fit[:topic]))
      append!(θ,vcat(map(x -> params(x)[1],topicppd(fit[:topic][i]))...));
  end
  writecsv(string(path,"topicparams.csv"),θ);
  writecsv(string(path,"paramlengths.csv"),
  map(x -> length(params(x)[1]),topicppd(fit[:topic][1])));
end

function sim_dat(n::Int64=100,nd::Int64=50,K::Int64=5,d::Int64=10,l::Int64=3,p::Int64=1,
  σ_μ::Float64=0.5,σ_η::Float64=0.25,σ_Θ::Float64=1.0,σ_β::Float64=0.0);

  sim = Dict{Symbol,Union{AbstractArray,Dict{Symbol,Float64}}}();
  nd = repeat([nd],inner=[n]);

  y = Vector{Array{Int64,2}}(n);
  μ = randn(K)*σ_μ;

  β = randn(K,p)*σ_β;
  X = randn(p,n);

  θ = rand(Normal(0,σ_Θ),l,d,K);
  θ = mapslices(softmax,θ,1)
  η = Array{Float64}(K,n);
  nk = Array{Int64}(K,n);
  for i in 1:n
    η[:,i] = randn(K).*σ_η .+ μ .+ β*X[:,i];
    nk[:,i] = rand(Multinomial(nd[i],softmax(η[:,i])));
    y[i] = Array{Int64,2}(d,0)
    for k in 1:K
      if nk[k,i] > 0
        y[i] = hcat(y[i],hcat(map(j -> mapslices(x -> rand(Categorical(x)),θ[:,:,k],1)[:],1:nk[k,i])...));
      end
    end
  end

  sim[:y] = y;
  sim[:Θ] = θ;
  sim[:η] = η;
  sim[:X] = X;
  sim[:β] = β;
  sim[:μ] = μ;

  return sim
end
