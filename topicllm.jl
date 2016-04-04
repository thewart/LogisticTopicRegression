function topiclmm{T}(y::Array{T},pss0::VectorPosterior,K::Int,iter::Int=1000)

  ## initialize
  pss = Array{VectorPosterior}(K);
  map!(k -> deepcopy(pss0),pss,1:K);
  n = length(y);
  z = rand(1:K,n);
  nk = map(k -> countnz(z.==k),1:K);
  for i in 1:n addsample!(pss[z[i]],y[i]); end

  σ2_η = 1.0;
  μ_η = 0.0;
  η = zeros(K);
  λ = Array{Float64}(K);

  zout = Array{Int}(n,iter);
  ηout = Array{Float64}(K,iter);
  logpost = Array{Float64}(K);
  π = Array{Float64}(K);

  for t in 1:iter
    ##  sample topic memberships
    softmax!(π,η);
    for i in 1:n
      pullsample!(pss[z[i]],y[i]);
      nk[z[i]] -= 1;

      for k in 1:K
        logpost[k] = log(π[k]);
        logpost[k] += lppd(pss[k],y[i]);
      end
      logpostnorm = logpost - logsumexp(logpost);
      z[i] = findfirst(rand(Multinomial(1,exp(logpostnorm))));

      addsample!(pss[z[i]],y[i]);
      nk[z[i]] += 1;
    end

    ## sample η
    for k in 1:K
      c = logsumexp(η[setdiff(1:K,k)]);
      ρ = η[k] - c;
      λ = rpolyagamma(ρ,n);

      v = inv(λ + 1/σ2_η);
      η[k] = randn()*sqrt(v) + v*(μ_η/σ2_η + nk[k]-n/2 + λ*c);
    end

    zout[:,t] = z;
    ηout[:,t] = η;
  end

  return zout,ηout,pss
end
