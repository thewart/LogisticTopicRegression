function topiclmm{T<:Real}(y::Vector{Array{T,2}},X::Array{Float64,2},pss0::VectorPosterior,K::Int,
                           hy::Dict{Symbol,Float64}=hyperparameter();
                           iter::Int=1000,thin::Int=1,report_loglik::Bool=false)

  ## initialize
  Base.Test.@test maximum(pss0.span[length(pss0)]) == size(y[1])[1];
  n = length(y);
  Base.Test.@test size(X)[2] == n;
  p = size(X)[1];

  ν0_σ2η = hy[:ν0_σ2η];
  σ0_σ2η = hy[:σ0_σ2η];
  τ0_τ = hy[:τ0_τ];
  ν0_τ = hy[:ν0_τ];
  τ_β = hy[:τ_β];

  #K0 = rand(round(Int64,K/2):K);
  K0 = K;
  topic = Vector{VectorPosterior}(K);
  map!(k -> deepcopy(pss0),topic,1:K);
  nd = map(i -> size(y[i])[2],1:n);
  z = map(d -> rand(1:K0,size(d)[2]),y);
  nk = Array{Int64}(K,n);
  for i in 1:n nk[:,i] = map(k -> countnz(z[i].==k),1:K); end
  for i in 1:n
    for j in 1:nd[i] addsample!(topic[z[i][j]],y[i][:,j]); end
  end
  XtX = X'X;
  Lβ = inv( chol(X*X' + diagm(fill(τ_β,p)) ));
  Σβ = Lβ*Lβ';

  σ2prior = InverseGamma(0.5*ν0_σ2η,0.5*σ0_σ2η*ν0_σ2η);
  τprior = InverseGamma(0.5*ν0_τ,0.5*τ0_τ*ν0_τ);

  σ2_η = fill(1.0,K);
  τ_μ = 0.1;
  η = randn(K,n);
  μ_η = Array{Float64}(K);
  w = Array{Float64}(n);

  saveiter = thin:thin:iter;
  nsave = length(saveiter);
  iter = maximum(saveiter);

  post = Dict{Symbol,Union{AbstractArray,Dict{Symbol,Float64}}}();
  post[:z] = Vector{Array{Int,2}}(n);
  for i in 1:n post[:z][i] = Array{Int}(nd[i],nsave); end
  post[:μ] = Array{Float64}(K,nsave);
  post[:σ2] = Array{Float64}(K,nsave);
  post[:τ] = Vector{Float64}(nsave);
  post[:topic] = Array{VectorPosterior,2}(K,nsave);
  post[:η] = Array{Float64}(K,n,nsave);
  post[:β] = Array{Float64}(p,K,nsave);
  if report_loglik post[:loglik] = Array{Float64,2}(n,nsave); end
  post[:hyperparameter] = hy;

  ηcpd = Vector{MultivariateNormal}(K);
  σ2cpd = Vector{InverseGamma}(K);
  τcpd = Vector{InverseGamma}(K);

  logpost = Array{Float64}(K);
  π = Array{Float64}(K);

  iΣ = makecov(XtX,τ_μ,τ_β);

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
          logpost[k] += lppd(y[i][:,j],topic[k]);
        end
        logpostnorm = logpost - logsumexp(logpost);
        z[i][j] = rand(Categorical(exp(logpostnorm)));

        addsample!(topic[z[i][j]],y[i][:,j]);
        nk[z[i][j],i] += 1;
      end
    end

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

      L = inv(chol(Hermitian(iΣ_k)));
      ηk = L*L'*w + L*randn(n);
      η[k,:] = ηk;

      ## sample variance
      a = 0.5(n+ν0_σ2η);
      b = 0.5(σ0_σ2η*ν0_σ2η + dot(ηk,iΣ*ηk));
      σ2_η[k] = rand(InverseGamma(a,b));

      ## sample mean
      σ2hat = inv(1/(τ_μ*σ2_η[k]) + n/σ2_η[k]);
      μhat = σ2hat*sum(η[k,:])/σ2_η[k];
      μ_η[k] = randn()*sqrt(σ2hat) + μhat;

    end
    ## sample variance of means
    a = 0.5(K+ν0_τ);
    b = 0.5(τ0_τ*ν0_τ + sum(μ_η.^2./σ2_η));
    τ_μ = rand(InverseGamma(a,b));

    iΣ = makecov(XtX,τ_μ,τ_β);

    ## save samples
    if t ∈ saveiter
      j = findin(saveiter,t)[1];
      for i in 1:n
        post[:z][i][:,j] = z[i];
        if report_loglik, post[:loglik][i,j] =
          sum(lppd(y[i],topic,softmax(η[:,i]))); end
      end
      for k in 1:K
        post[:β][:,k,j] = Σβ*X*(η[k,:] .- μ_η[k]) + sqrt(σ2_η[k]).*Lβ*randn(p);
        #  post[:lpθ][j] += logpdf(MvNormalCanon(iΣ./σ2_η[k]),η[k,:])[1] +
        #                  logpdf(σ2prior,σ2_η[k]);
      end
      #post[:lpθ][j] += logpdf(τprior,τ_μ);
      post[:topic][:,j] = deepcopy(topic);
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
      V[j,i] = XtX[j,i].*τ_β + (i==j ? 1+τ : τ);
    end
  end

  L = inv(chol(V));
  return L*L'
end

function hyperparameter(;ν0_σ2η=1.0,σ0_σ2η = 1.0,
  τ0_τ = 0.25,ν0_τ = 1.0,τ_β = 1.0)
  return Dict(:ν0_σ2η => ν0_σ2η, :σ0_σ2η => σ0_σ2η, :τ0_τ => τ0_τ,
  :ν0_τ => ν0_τ, :τ_β => τ_β)
end

function refβ(β::Array{Float64,2},refk::Int64)
  return β .- β[:,refk]
end
refβ(β::Array{Float64,3},refk::Int64) = β .- β[:,refk,:]
refβ(β::Array{Float64},μ::Array{Float64,1}) = refβ(β,findmax(μ)[2]);
refβ(β::Array{Float64},μ::Array{Float64,2}) = refβ(β,findmax(mean(μ,2))[2]);

function writefit(fit::Dict{Symbol,AbstractArray},path::String)
  if !isdir(path) mkpath(path); end
  n = length(fit[:z]);
  K,nsave = size(fit[:topic]);
  p = size(fit[:β])[1];
    b = length(fit[:topic][1]);
  writecsv(string(path,"dims.csv"),[n,K,p,b,nsave]);
  writecsv(string(path,"sigma.csv"),hcat(fit[:τ],fit[:σ2]'));
  writecsv(string(path,"beta.csv"),fit[:β][:]);
  writecsv(string(path,"eta.csv"),fit[:η][:]);
  writecsv(string(path,"loglik.csv"),fit[:loglik]');
  writecsv(string(path,"mu.csv"),fit[:μ]');
  writecsv(string(path,"z.csv"),vcat(fit[:z]...));
  writecsv(string(path,"nd.csv"),map(x -> size(x)[1],fit[:z]));
  writecsv(string(path,"topicmean.csv"),
            mapreduce(y -> mean(topicppd(y)),vcat,fit[:topic]));
end
