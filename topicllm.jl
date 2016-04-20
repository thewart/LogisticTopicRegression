function topiclmm{T<:Vector{Array{Real,2}}}(y::T,pss0::VectorPosterior,K::Int,iter::Int=1000)

  ## initialize
  pss = Array{VectorPosterior}(K);
  map!(k -> deepcopy(pss0),pss,1:K);
  n = length(y);
  nd = map(i -> size(y[i])[2],1:n);
  z = map(d -> rand(1:K,size(d)),y);
  nk = Array{Int64}(K,n);
  for i in 1:n nk[:,i] = map(k -> countnz(z[i].==k),1:K); end
  for i in 1:n
    for j in 1:nd[i] addsample!(pss[z[i][j]],y[i][:,j]); end
  end

  ν0_σ2η = 1.0;
  σ0_σ2η = 1.0;
  τ0_τ = 1.0;
  ν0_τ = 0.1;

  σ2_η = fill(1.0,K);
  τ_μ = 0.1;
  η = randn(K,n);
  μ_η = Array{Float64}[K];
  λ = Array{Float64}(n);

  zout = Array{Int}(n,iter);
  ηout = Array{Float64}(K,iter);
  logpost = Array{Float64}(K);
  π = Array{Float64}(K);

  for t in 1:iter

    ##  sample topic memberships
    for i in 1:n
      softmax!(π,η[:,i]);
      for j in 1:nd[i]
        zj = z[i][j];
        pullsample!(pss[zj],y[i][:,j]);
        nk[zj,i] -= 1;

        for k in 1:K
          logpost[k] = log(π[k]);
          logpost[k] += lppd(pss[k],y[i][:,j]);
        end
        logpostnorm = logpost - logsumexp(logpost);
        z[i][j] = findfirst(rand(Multinomial(1,exp(logpostnorm))));

        addsample!(pss[zj],y[i][:,j]);
        nk[zj,i] -= 1;
      end
    end

    iΣ_ηk = myinvcov(n,τ_μ);
    for k in 1:K
      ## sample η and λ
      iV = iΣ_ηk./σ2_η[k];
      for i in 1:n
        c = logsumexp(η[setdiff(1:K,k),i]);
        ρ = η[k,i] - c;
        λ[i] = rpolyagamma(ρ,nd[i]);
        w[i] = nk[k,i] + nd[i]/2 + c*λ[i];
        iV[i,i] += λ[i];
      end
      L = inv(chol(iV));
      η[k,:] = L*L'*w + L*randn(n);

      ## sample variance
      α = 0.5(n+ν0_σ2η);
      β = 0.5(σ0_σ2η*ν0_σ2η + η[k,:]*iΣ_ηk*η[k,:]');
      σ2_η[k] = rand(InverseGamma(α,β),1);

      ## sample mean
      σ2hat = inv(1/(τ_μ*σ2_η[k]) + n/σ2_η[k]);
      μhat = σ2hat*sum(η[k,:])/σ2_η[k];
      μ_η[k] = randn()*sqrt(σ2hat) + μhat;
    end
    ## sample variance of means
    α = 0.5(K+ν0_τ);
    β = 0.5(τ0_τ*ν0_τ + sum(μ_η.^2./σ2_η));
    τ_μ = rand(InverseGamma(α,β));

    ## save samples
    #zout[:,t] = z;
    pssout[t] = deepcopy(pss);
    ηout[:,:,t] = η;

  end

  return zout,ηout,pss
end

function myinvcov(d::Int64,τ::Float64)
  M = Array{Float64}(d,d);
  a = inv(1+τ) + τ^2*(d-1)/(d*τ^2 + (d+1)*τ +1);
  b = -τ/(1+d*τ);
  for i in 1:d
    for j in 1:i
      i==j ? M[i,i] = a : M[i,j] = M[j,i] = b;
    end
  end
  return M
end
