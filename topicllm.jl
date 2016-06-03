function topiclmm{T<:Real}(y::Vector{Array{T,2}},X::Array{Float64,2},pss0::VectorPosterior,K::Int,
                           iter::Int=1000,thin::Int=1)

  ## initialize
  Base.Test.@test maximum(pss0.span) == size(y[1])[1];
  n = length(y);
  Base.Test.@test size(X)[2] == n;

  topic = Vector{VectorPosterior}(K);
  map!(k -> deepcopy(pss0),topic,1:K);
  nd = map(i -> size(y[i])[2],1:n);
  z = map(d -> rand(1:K,size(d)[2]),y);
  nk = Array{Int64}(K,n);
  for i in 1:n nk[:,i] = map(k -> countnz(z[i].==k),1:K); end
  for i in 1:n
    for j in 1:nd[i] addsample!(topic[z[i][j]],y[i][:,j]); end
  end
  XtX = X'X;

  ν0_σ2η = 1.0;
  σ0_σ2η = 1.0;
  τ0_τ = 0.25;
  ν0_τ = 1.0;
  τ_β = 0.5;

  σ2_η = fill(1.0,K);
  τ_μ = 0.1;
  η = randn(K,n);
  μ_η = Array{Float64}(K);
  w = Array{Float64}(n);

  saveiter = thin:thin:iter;
  nsave = length(saveiter);
  iter = maximum(saveiter);

  post = Dict{Symbol,AbstractArray}();
  post[:z] = Vector{Array{Int,2}}(n);
  for i in 1:n post[:z][i] = Array{Int}(nd[i],nsave); end
  post[:μ] = Array{Float64}(K,nsave);
  post[:σ2] = Array{Float64}(K,nsave);
  post[:τ] = Vector{Float64}(nsave);
  post[:topic] = Vector{typeof(topic)}(nsave);
  post[:η] = Array{Float64}(K,n,nsave);

  logpost = Array{Float64}(K);
  π = Array{Float64}(K);

  for t in 1:iter

    ##  sample topic memberships
    for i in 1:n
      softmax!(π,η[:,i]);
      for j in 1:nd[i]
        zj = z[i][j];
        pullsample!(topic[zj],y[i][:,j]);
        nk[zj,i] -= 1;

        for k in 1:K
          logpost[k] = log(π[k]);
          logpost[k] += lppd(topic[k],y[i][:,j]);
        end
        logpostnorm = logpost - logsumexp(logpost);
        z[i][j] = findfirst(rand(Multinomial(1,exp(logpostnorm))));

        addsample!(topic[z[i][j]],y[i][:,j]);
        nk[z[i][j],i] += 1;
      end
    end

    iΣ = makecov(XtX,τ_μ,τ_β);
    for k in 1:K
      ## sample η and λ
      iΣ_k = iΣ./σ2_η[k];
      for i in 1:n
        c = logsumexp(η[setdiff(1:K,k),i]);
        ρ = η[k,i] - c;
        λ = rpolyagamma(ρ,nd[i]);
        w[i] = nk[k,i] - nd[i]/2 + c*λ;
        iΣ_k[i,i] += λ;
      end
      L = inv(chol(iΣ_k));
      η[k,:] = L*L'*w + L*randn(n);

      ## sample variance
      α = 0.5(n+ν0_σ2η);
      #β = 0.5(σ0_σ2η*ν0_σ2η + (η[k,:]*iΣ_ηk*η[k,:]')[1]);
      β = 0.5(σ0_σ2η*ν0_σ2η + (η[k,:]*iΣ_k*η[k,:]')[1]);
      σ2_η[k] = rand(InverseGamma(α,β));

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
    if t ∈ saveiter
      j = findin(saveiter,t)[1];
      for i in 1:n post[:z][i][:,j] = z[i]; end
      post[:topic][j] = deepcopy(topic);
      post[:η][:,:,j] = η;
      post[:μ][:,j] = μ_η;
      post[:σ2][:,j] = σ2_η;
      post[:τ][j] = τ_μ;
    end

  end

  return post
end

function makecov(XtX::Array{Float64,2},τ::Float64,τ_β::Float64)
  n = size(XtX)[1];
  V = Array{Float64,2}(n,n);
  for i in 1:n
    for j in 1:n
      V[i,j] = XtX[i,j].*τ_β + (i==j ? 1+τ : τ);
    end
  end

  L = inv(chol(V));
  return L*L'
end
