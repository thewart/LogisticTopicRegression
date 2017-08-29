## random genetic (or other) effects
function topiclmm{T<:Real}(y::Vector{Array{T,2}},X::Array{Float64,2},Z::Array{Float64,2},pss0::VectorPosterior,K::Int,
                           hy::Dict{Symbol,Float64}=hyperparameter();
                           zinit::Vector{Vector{Int64}}=Vector{Vector{Int64}}(),
                           θinit::Dict{Symbol,Array{Float64}}=Dict{Symbol,Array{Float64}}(),
                           iter::Int=1000,thin::Int=1)

  ## initialize
  Base.Test.@test maximum(pss0.span[length(pss0)]) == size(y[1])[1];
  n = length(y);
  Base.Test.@test size(X)[2] == n;
  p = size(X)[1];

  ν0_σ2η = hy[:ν0_σ2η];
  σ0_σ2η = hy[:σ0_σ2η];
  τ0_τ = hy[:τ0_τ];
  ν0_τ = hy[:ν0_τ];
  τ0_u = hy[:τ0_u];
  ν0_u = hy[:ν0_u];
  τ_β = hy[:τ_β];

  θinit = merge(init_params(K,n),θinit);
  σ2_η = θinit[:σ2_η];
  τ_μ = θinit[:τ_μ][1];
  η = θinit[:η];
  τ_u = θinit[:τ_u];

  #K0 = rand(round(Int64,K/2):K);
  topic = Vector{VectorPosterior}(K);
  map!(k -> deepcopy(pss0),topic,1:K);
  nd = size.(y,2);

  if isempty(zinit)
    nk, z = init_topic!(topic,η,y);
  else
    nk = init_topic!(topic,z,y);
  end
  ZtZ = Z'Z;
  XtX = X'X;

  Lβ = inv( chol( X*X' + I*inv(τ_β)));
  ΣβX = Lβ*Lβ'*X;

  μ_η = Array{Float64}(K);
  w = Array{Float64}(n);
  u = Array{Float64}(n,K);
  β = Array{Float64,2}(p,K);

  saveiter = thin:thin:iter;
  nsave = length(saveiter);
  iter = maximum(saveiter);

  post = init_post(n,nd,K,p,nsave,true);
  post[:hyperparameter] = hy;

  for t in 1:iter

    ##  sample topic memberships
    for i in 1:n
      sample_z!(z[i],topic,view(nk,:,i),softmax(η[:,i]),y[i],K);
    end

    for k in 1:K
      ## sample η and λ
      iΣ = makecov(XtX,ZtZ,τ_β,τ_u[k],τ_μ);
      iΣ_k = iΣ./σ2_η[k];

      nk = sample_η(η[k,:],η[vcat(1:(k-1),(k+1):end),:],iΣ_k,nd,nk[k,:]);
      η[k,:] = ηk;

      ## sample variance
      σ2_η[k] = sample_variance(ηk,iΣ,ν0_σ2η,σ0_σ2η);

      ## sample mean
      Vμ = inv(τ_u[k]*ZtZ + τ_β*XtX + I);
      σ2hat = σ2_η[k]/(1/τ_μ + sum(Vμ));
      μhat = σ2hat/σ2_η[k]*sum(Vμ*η[k,:]);
      μ_η[k] = randn()*sqrt(σ2hat) + μhat;

      resid = (η[k,:] .- μ_η[k]);

      innerV = inv(XtX*τ_β + I);
      Lu = inv( chol(Hermitian(Z*innerV*Z' + I*inv(τ_u[k]))));
      u[:,k] = Lu*Lu'*Z*innerV*resid + sqrt(σ2_η[k]).*Lu*randn(n);

      τ_u[k] = sample_variance(u[:,k],σ2_η[k],ν0_u,τ0_u);

      resid .-= Z'u[:,k];
      β[:,k] = ΣβX*resid + sqrt(σ2_η[k]).*Lβ*randn(p);

    end
    ## sample variance of means
    τ_u = sample_variance(μ_η,σ2_η,ν0_τ,τ0_τ);

    ## save samples
    if t ∈ saveiter
      j = findin(saveiter,t)[1];
      for i in 1:n
        post[:z][i][:,j] = z[i];
      end
      post[:β][:,:,j] = β;
      post[:u][:,:,j] = u;
      post[:topic][:,j] = deepcopy(topic);
      post[:η][:,:,j] = η;
      post[:μ][:,j] = μ_η;
      post[:σ2][:,j] = σ2_η;
      post[:τ][j] = τ_μ;
      post[:τ_u][:,j] = τ_u;
    end
  end
  return post
end

## vanilla
function topiclmm{T<:Real}(y::Vector{Array{T,2}},X::Array{Float64,2},pss0::VectorPosterior,K::Int,
                           hy::Dict{Symbol,Float64}=hyperparameter();
                           zinit::Vector{Vector{Int64}}=Vector{Vector{Int64}}(),
                           θinit::Dict{Symbol,Array{Float64}}=Dict{Symbol,Array{Float64}}(),
                           iter::Int=1000,thin::Int=1)

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

  θinit = merge(init_params(K,n),θinit);
  σ2_η = θinit[:σ2_η];
  τ_μ = θinit[:τ_μ][1];
  η = θinit[:η];

  topic = Vector{VectorPosterior}(K);
  map!(k -> deepcopy(pss0),topic,1:K);
  nd = size.(y,2);

  if isempty(zinit)
    nk, z = init_topic!(topic,η,y);
  else
    nk = init_topic!(topic,z,y);
  end

  XtX = X'X;
  Lβ = inv( chol( X*X' + I*inv(τ_β)));
  ΣβX = Lβ*Lβ'*X;
  invaddIτXtX = inv(I+τ_β*XtX);
  suminvaddIτXtX = sum(invaddIτXtX);

  μ_η = Array{Float64}(K);
  w = Array{Float64}(n);
  β = Array{Float64,2}(p,K);

  saveiter = thin:thin:iter;
  nsave = length(saveiter);
  iter = maximum(saveiter);

  post = init_post(n,nd,K,p,nsave,false);
  post[:hyperparameter] = hy;

  iΣ = makecov(XtX,τ_μ,τ_β);
  for t in 1:iter

    ##  sample topic memberships
    for i in 1:n
      sample_z!(z[i],topic,view(nk,:,i),softmax(η[:,i]),y[i],K);
    end

    for k in 1:K
    ## sample η and λ
      iΣ_k = iΣ./σ2_η[k];

      ηk = sample_η(η[k,:],η[vcat(1:(k-1),(k+1):end),:],iΣ_k,nd,nk[k,:]);
      η[k,:] = ηk;

      σ2_η[k] = sample_variance(ηk,iΣ,ν0_σ2η,σ0_σ2η);

      σ2hat = σ2_η[k]/(1/τ_μ + suminvaddIτXtX);
      μhat = σ2hat/σ2_η[k]*sum(invaddIτXtX*η[k,:]);
      μ_η[k] = randn()*sqrt(σ2hat) + μhat;

      β[:,k] = ΣβX*(η[k,:] .- μ_η[k]) + sqrt(σ2_η[k]).*Lβ*randn(p);
    end

    τ_u = sample_variance(μ_η./σ2_η,1.0,ν0_τ,τ0_τ);
    ## save samples
    if t ∈ saveiter
      j = findin(saveiter,t)[1];
      for i in 1:n
        post[:z][i][:,j] = z[i];
      end
      post[:β][:,:,j] = β;
      post[:topic][:,j] = deepcopy(topic);
      post[:η][:,:,j] = η;
      post[:μ][:,j] = μ_η;
      post[:σ2][:,j] = σ2_η;
      post[:τ][j] = τ_μ;
    end
  end
  return post
end
