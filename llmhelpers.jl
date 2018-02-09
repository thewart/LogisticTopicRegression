function makecov(XtX::Array{Float64,2},ZtZ::Array{Float64,2},
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

function makecov(XtX::Array{Float64,2},τ_β::Float64,τ_μ::Float64)

  n = size(XtX)[1];
  V = Array{Float64,2}(n,n);
  for i in 1:n
    for j in 1:n
      V[j,i] = XtX[j,i].*τ_β + (i==j ? 1+τ_μ : τ_μ);
    end
  end

  return inv(V)
end

# function init_post(n::Int64,nd::Vector{Int64},K::Int64,p::Int64,nsave::Int64,isu::Bool)
#   post = Dict{Symbol,Union{AbstractArray,Dict{Symbol,Float64}}}();
#   post[:z] = Vector{Array{Int,2}}(n);
#   for i in 1:n post[:z][i] = Array{Int}(nd[i],nsave); end
#   post[:μ] = Array{Float64}(K,nsave);
#   post[:σ2] = Array{Float64}(K,nsave);
#   post[:τ] = Vector{Float64}(nsave);
#   post[:topic] = Array{VectorPosterior,2}(K,nsave);
#   post[:η] = Array{Float64}(K,n,nsave);
#   post[:β] = Array{Float64}(p,K,nsave);
#   return post
#   if isu
#     post[:u] = Array{Float64}(n,K,nsave);
#     post[:τ_u] = Array{Float64,2}(K,nsave);
#   end
# end

function hyperparameter(;ν0_σ2η=1.0,σ0_σ2η = 1.0,
  τ0_τ = 0.25,ν0_τ = 1.0,τ0_u = 0.25,ν0_u=1.0,τ_β=1.0)
  return Dict(:ν0_σ2η => ν0_σ2η, :σ0_σ2η => σ0_σ2η, :τ0_τ => τ0_τ,
  :ν0_τ => ν0_τ, :τ0_u => τ0_u, :ν0_u => ν0_u, :τ_β => τ_β)
end

function init_topic(pss0,K,z,y)
  n = length(y);
  nd = size.(y,2);

  topic = Vector{typeof(pss0)}(K);
  map!(k -> deepcopy(pss0),topic,1:K);
  for i in 1:n
    for j in 1:nd[i] addsample!(topic[z[i][j]],y[i][:,j]); end
  end
end

function init_z(η,K,docrng)
  n = length(grps);
  z = Vector{Int64}(n);
  for i in 1:n z[docrng[i]] = sample(1:K,Weights(softmax(η[:,i])),length(docrng[i])); end
  return z
end

function refβ(β::Array{Float64,2},refk::Int64)
  return β .- β[:,refk:refk]
end
refβ(β::Array{Float64,3},refk::Int64) = β .- β[:,refk:refk,:]
refβ(β::Array{Float64},μ::Array{Float64,1}) = refβ(β,findmax(μ)[2]);
refβ(β::Array{Float64},μ::Array{Float64,2}) = refβ(β,findmax(mean(μ,2))[2]);

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

function mlinit_dirichlet(y::Vector{Array{Int64,2}},K::Int64,d::Int64,L::Vector{Int64},nruns::Int64=10,noise::Float64=0.5,maxiter::Int64=100);

  Y = hcat(y...)';
  @rput K d L Y noise maxiter
  R"library(rstan)"
  R"topetho <- stan_model('~/code/topetho/categorical_topetho.stan')";
  R"dat <- list(Y=Y,n=nrow(Y),K=K,B=d,Bs=L,alpha_p=1,alpha_t=1)";

  ll = Vector{Float64}(nruns);
  r = Vector{Array{Float64,2}}(nruns);
  @time for i in 1:nruns
    R"init <- list(pi=gtools::rdirichlet(1,alpha = rep(1,dat$K)) %>% as.vector(),
    theta_raw=apply(Y,2,function(x) table(x) %>% prop.table()) %>% unlist() %>% matrix(nrow=dat$K,ncol=sum(dat$Bs),byrow = T))";
    R"init$theta_raw <- init$theta_raw * pmax(1-rnorm(length(init$theta_raw),sd=noise),0.01)";
    R"optout <- optimizing(topetho,dat,verbose=F,init=init,iter=400)";
    ll[i] = rcopy(R"optout$value");
    r[i] = rcopy(R"optout$par[str_detect(names(optout$par),'^r\\[')] %>% matrix(nrow=dat$K,ncol=dat$n,byrow = T)")
  end

  maxr = r[findmax(ll)[2]];

  zflat = mapslices(x -> rand(Categorical(x)),maxr,1);
  z = Vector{Vector{Int64}}(length(y));
  nd = map(x -> size(x)[2],y);
  guh = 1;
  for i in 1:length(y)
    z[i] = zflat[guh:(guh+nd[i]-1)];
    guh += nd[i];
  end

  return z
end
