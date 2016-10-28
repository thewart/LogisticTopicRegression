function simplemix{T<:Real}(y::Array{T,2},pss0::VectorPosterior,K::Int,
                           hy::Dict{Symbol,Float64}=hyperparameter();zinit::Vector{Int64}=Vector{Int64}(),
                           iter::Int=1000,thin::Int=1,report_loglik::Bool=false)

  ## initialize
  Base.Test.@test maximum(pss0.span[length(pss0)]) == size(y)[1];
  n = size(y)[2];

  #K0 = rand(round(Int64,K/2):K);
  topic = Vector{VectorPosterior}(K);
  map!(k -> deepcopy(pss0),topic,1:K);

  if isempty(zinit)
    wv = weights(rand(Dirichlet(K,rand(Uniform(1/K,1)))));
    z = sample(1:K,wv,size(y)[2]);
  else
    z = zinit;
  end
  nk = map(k -> countnz(z.==k),1:K);

  for i in 1:n addsample!(topic[z[i]],y[:,i]); end

  σ2_η = 1.0;
  μ_η = 0.0;
  η = zeros(K);  w = Array{Float64}(n);

  saveiter = thin:thin:iter;
  nsave = length(saveiter);
  iter = maximum(saveiter);

  post = Dict{Symbol,Union{AbstractArray,Dict{Symbol,Float64}}}();
  post[:z] = Array{Int,2}(n,nsave);
  #post[:μ] = Array{Float64}(K,nsave);
  #post[:σ2] = Array{Float64}(K,nsave);
  #post[:τ] = Vector{Float64}(nsave);
  post[:topic] = Array{VectorPosterior,2}(K,nsave);
  #post[:η] = Array{Float64}(K,n,nsave);
  post[:η] = Array{Float64}(K,nsave);
  #post[:β] = Array{Float64}(p,K,nsave);

  if report_loglik post[:loglik] = Array{Float64,2}(n,nsave); end
  post[:hyperparameter] = hy;

  logpost = Array{Float64}(K);
  π = Array{Float64}(K);

  #iΣ = makecov(XtX,τ_μ,τ_β);

  for t in 1:iter

    softmax!(π,η);
    ##  sample topic memberships
    for i in 1:n
      pullsample!(topic[z[i]],y[:,i]);
      nk[z[i]] -= 1;

      for k in 1:K
        logpost[k] = log(π[k]);
        logpost[k] += lppd(y[:,i],topic[k]);
      end
      logpostnorm = logpost - logsumexp(logpost);
      z[i] = findfirst(rand(Multinomial(1,exp(logpostnorm))));

      addsample!(topic[z[i]],y[:,i]);
      nk[z[i]] += 1;
    end

      ## sample η and λ
    for k in 1:K
        c = logsumexp(η[setdiff(1:K,k)]);
        ρ = η[k] - c;
        λ = rpolyagamma(ρ,n);
        v = inv(λ + 1/σ2_η);
        η[k] = randn()*sqrt(v) + v*(μ_η/σ2_η + nk[k]-n/2 + λ*c);
    end

      ## sample variance
      #a = 0.5(n+ν0_σ2η);
      #b = 0.5(σ0_σ2η*ν0_σ2η + dot(ηk,iΣ*ηk));
      #σ2_η[k] = rand(InverseGamma(a,b));

      ## sample mean
      #σ2hat = σ2_η[k]/(1/τ_μ + suminvaddIτXtX);
      #μhat = σ2hat/σ2_η[k]*sum(invaddIτXtX*η[k,:]);
      #σ2hat = inv(1/(τ_μ*σ2_η[k]) + n/σ2_η[k]);
      #μhat = σ2hat*sum(η[k,:])/σ2_η[k];
      #μ_η[k] = randn()*sqrt(σ2hat) + μhat;

    ## sample variance of means
    #a = 0.5(K+ν0_τ);
    #b = 0.5(τ0_τ*ν0_τ + sum(μ_η.^2./σ2_η));
    #τ_μ = rand(InverseGamma(a,b));

    #iΣ = makecov(XtX,τ_μ,τ_β);

    ## save samples
    if t ∈ saveiter
      j = findin(saveiter,t)[1];
      post[:z][:,j] = z;
      #for k in 1:K
    #    post[:β][:,k,j] = ΣβX*(η[k,:] .- μ_η[k]) + sqrt(σ2_η[k]).*Lβ*randn(p);
        #  post[:lpθ][j] += logpdf(MvNormalCanon(iΣ./σ2_η[k]),η[k,:])[1] +
        #                  logpdf(σ2prior,σ2_η[k]);
      #end
      #post[:lpθ][j] += logpdf(τprior,τ_μ);
      post[:topic][:,j] = deepcopy(topic);
      post[:η][:,j] = η;
      #post[:μ][:,j] = μ_η;
      #post[:σ2][:,j] = σ2_η;
      #post[:τ][j] = τ_μ;
    end

  end

  return post
end
