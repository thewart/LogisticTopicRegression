using LogTopReg, StatsFuns, StatsBase

σ_μ = 0.5;
σ = 0.25;
σ_θ = 1.0;
σ_β = 0.0;

p=1;
n = 200;
nd = 100;
K = [5,10,20,40,80];
d = 30;
l = 3;

pss0 = VectorPosterior(CategoricalPosterior(l),d);
ns = 500;
warmup = 100;

Θmse = Vector{Vector{Float64}}(length(K));
pmse = Array{Float64,2}(n,length(K));
Θr2 = Vector{Vector{Float64}}(length(K));
pr2 = Array{Float64,2}(n,length(K));
topiccor = Vector{Vector{Float64}}(length(K));

for i in 1:length(K)
  sim = sim_dat(n,nd,K[i],d,l,p,σ_μ,σ,σ_θ,σ_β);
  z = mlinit_dirichlet(sim[:y],K[i],d,fill(l,d),10);
  fit = topiclmm(sim[:y],sim[:X],pss0,K[i],hyperparameter(τ_β=1e-10),zinit=z,iter=ns);

  θhat = zeros(size(sim[:Θ]));
  for k = 1:K[i]
    for t = (warmup+1):ns
      θhat[:,:,k] += hcat(map(mean,topicpd(fit[:topic][k,t]))...);
    end
  end
  θhat = θhat./ns;

  topicord = Vector{Int64}(K[i]);
  topiccor[i] = Vector{Float64}(K[i]);
  for k in 1:K[i]
    avail = setdiff(1:K[i],topicord[1:(k-1)]);
    topiccor[i][k],topicord[k] = findmax(mapslices(x -> cor(vec(x),vec(θhat[:,:,k])),sim[:Θ][:,:,avail],(1,2)));
    topicord[k] = avail[topicord[k]];
  end

  Θmse[i] = mapslices(sumabs2,θhat-sim[:Θ][:,:,topicord],(1,2))[:]./(l*d);
  Θr2[i] = 1-mapslices(sumabs2,θhat-sim[:Θ][:,:,topicord],(1,2))[:]./(mapslices(x->sumabs2(x-mean(x)),sim[:Θ][:,:,topicord],(1,2))[:]);

  ηhat = mapslices(mean,fit[:η][:,:,(warmup+1):ns],3)[:,:,1];
  phat = mapslices(softmax,ηhat,1);
  pi = mapslices(softmax,sim[:η][topicord,:],1);
  pmse[:,i] = mapslices(sumabs2,phat-pi,1)./n;
  pr2[:,i] = 1-mapslices(sumabs2,phat-pi,1)./mapslices(x->sumabs2(x-mean(x)),pi,1);

end
