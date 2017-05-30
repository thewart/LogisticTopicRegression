function simulation(n::Int64=100,nd::Int64=50,K::Int64=5,d::Int64=10,l::Int64=3,p::Int64=1,
  σ_μ::Float64=0.5,σ_η::Float64=0.25,σ_Θ::Float64=1.0,σ_β::Float64=0.0);

  nd = repeat([nd],inner=[n]);

  y = Vector{Array{Int64,2}}(n);
  μ = randn(K)*σ_μ;

  β = randn(K,p)*σ_β;
  X = randn(p,n);

  θ = rand(Normal(0,σ_θ),l,d,K);
  θ = mapslices(softmax,θ,1)
  η = Array{Float64}(K,n);
  nk = Array{Int64}(K,n);
  for i in 1:n
    η[:,i] = randn(K).*σ .+ μ .+ β*X[:,i];
    nk[:,i] = rand(Multinomial(nd[i],softmax(η[:,i])));
    y[i] = Array{Int64,2}(d,0)
    for k in 1:K
      if nk[k,i] > 0
        y[i] = hcat(y[i],hcat(map(j -> mapslices(x -> rand(Categorical(x)),θ[:,:,k],1)[:],1:nk[k,i])...));
      end
    end
  end

  Y = hcat(y...)';
  @rput K d l Y
  R"library(rstan)"
  R"topetho <- stan_model('~/code/topetho/categorical_topetho.stan')";
  R"dat <- list(Y=Y,n=nrow(Y),K=K,B=d,Bs=rep(l,d),alpha_p=1,alpha_t=1)";

  runs = 10;
  ll = Vector{Float64}(runs);
  r = Vector{Array{Float64,2}}(runs);
  @time for i in 1:runs
    R"init <- list(pi=gtools::rdirichlet(1,alpha = rep(1,dat$K)) %>% as.vector(),theta_raw=apply(Y,2,function(x) table(x) %>% prop.table()))";
    R"init$theta_raw <- init$theta_raw * pmax(1-rnorm(length(init$theta_raw),sd=0.5),0.01)";
    R"optout <- optimizing(topetho,dat,verbose=F,iter=100)";
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

  pss0 = VectorPosterior(CategoricalPosterior(l),d);
  @time fit = topiclmm(y,X,pss0,K,hyperparameter(τ_β=1e-10),zinit=z,iter=ns);

  return fit, η, Θ
end
