## random genetic (or other) effects
function topiclmm{U<:PostPredSS}(y,X,Z,pss0::VectorPosterior{U},K::Int64;
                           hy::Dict{Symbol,Float64}=hyperparameter(),
                           init::TLMMsample=init_params(K,length(y),size(X,1)),
                           iter::Int=1000,thin::Int=1)

  ## initialize
  Base.Test.@test maximum(pss0.span[length(pss0)]) == size(y[1])[1];
  n = length(y);
  Base.Test.@test size(X)[2] == n;
  p = size(X,1);
  nd = size.(y,2);

  if !isdefined(init,:z)
      init.z = init_z(init.η,K,nd)
  end
  topic = init_topic(pss0,K,init.z,y);
  nk = hcat(countz.(init.z,K)...);

  ZtZ = Z'Z;
  XtX = X'X;
  Lβ = inv( chol( X*X' + I*inv(hy[:τ_β])));
  ΣβX = Lβ*Lβ'*X;

  saveiter = thin:thin:iter;
  nsave = length(saveiter);
  iter = maximum(saveiter);

  samples = Vector{TLMMsample}(nsave);
  tss = Matrix{VectorPosterior{U}}(K,nsave);
  s = init;

  for t in 1:iter

    ##  sample topic memberships
    for i in 1:n
      sample_z!(s.z[i],topic,view(nk,:,i),softmax(s.η[:,i]),y[i],K);
    end

    for k in 1:K
      ## sample η and λ
      iΣ = makecov(XtX,ZtZ,hy[:τ_β],s.τ_u[k],s.τ_μ);
      iΣ_k = iΣ./s.σ2_η[k];

      nk = sample_η(s.η[k,:],s.η[vcat(1:(k-1),(k+1):end),:],iΣ_k,nd,nk[k,:]);
      s.η[k,:] = ηk;

      ## sample variance
      s.σ2_η[k] = sample_variance(ηk,iΣ,hy[:ν0_σ2η],hy[:σ0_σ2η]);

      ## sample mean
      Vμ = inv(s.τ_u[k]*ZtZ + hy[:τ_β]*XtX + I);
      σ2hat = s.σ2_η[k]/(1/s.τ_μ + sum(Vμ));
      μhat = σ2hat/s.σ2_η[k]*sum(Vμ*ηk);
      s.μ[k] = randn()*sqrt(σ2hat) + μhat;

      resid = (s.η[k,:] .- s.μ[k]);

      innerV = inv(XtX*hy[:τ_β] + I);
      Lu = inv( chol(Hermitian(Z*innerV*Z' + I*inv(s.τ_u[k]))));
      s.u[:,k] = Lu*Lu'*Z*innerV*resid + sqrt(s.σ2_η[k]).*Lu*randn(n);

      s.τ_u[k] = sample_variance(s.u[:,k],s.σ2_η[k],hy[:ν0_u],hy[:τ0_u]);

      resid .-= Z's.u[:,k];
      s.β[:,k] = ΣβX*resid + sqrt(s.σ2_η[k]).*Lβ*randn(p);

    end
    ## sample variance of means
    s.τ_μ = sample_variance(s.μ,s.σ2_η,hy[:ν0_τ],hy[:τ0_τ]);

    ## save samples
    if t ∈ saveiter
      j = findin(saveiter,t)[1];
      samples[j] = deepcopy(s);
      tss[:,j] = deepcopy(topic);
    end
  end

  return TLMMfit(samples,tss,pss0,hy);
end

## vanilla
function topiclmm{U<:PostPredSS}(y,X,pss0::VectorPosterior{U},K::Int64;
                           hy::Dict{Symbol,Float64}=hyperparameter(),
                           init::TLMMsample=init_params(K),
                           iter::Int=1000,thin::Int=1)

  ## initialize
  Base.Test.@test maximum(pss0.span[length(pss0)]) == size(y[1])[1];
  n = length(y);
  Base.Test.@test size(X)[2] == n;
  p = size(X)[1];
  nd = size.(y,2);

  if !isdefined(init,:z)
      init.z = init_z(init.η,K,nd)
  end
  topic = init_topic(pss0,K,init.z,y);
  nk = hcat(countz.(init.z,K)...);

  XtX = X'X;
  Lβ = inv( chol( X*X' + I*inv(hy[:τ_β])));
  ΣβX = Lβ*Lβ'*X;
  invaddIτXtX = inv(I+hy[:τ_β]*XtX);
  suminvaddIτXtX = sum(invaddIτXtX);

  saveiter = thin:thin:iter;
  nsave = length(saveiter);
  iter = maximum(saveiter);

  samples = Vector{TLMMsample}(nsave);
  tss = Matrix{VectorPosterior{U}}(K,nsave);
  s = init;

  for t in 1:iter

      iΣ = makecov(XtX,hy[:τ_β],τ_μ);
    ##  sample topic memberships
    for i in 1:n
      sample_z!(s.z[i],topic,view(nk,:,i),softmax(s.η[:,i]),y[i],K);
    end

    for k in 1:K
    ## sample η and λ
      iΣ_k = iΣ./s.σ2_η[k];

      ηk = sample_η(s.η[k,:],s.η[vcat(1:(k-1),(k+1):end),:],iΣ_k,nd,nk[k,:]);
      s.η[k,:] = ηk;

      s.σ2_η[k] = sample_variance(ηk,iΣ,hy[:ν0_σ2η],hy[:σ0_σ2η]);

      σ2hat = s.σ2_η[k]/(1/s.τ_μ + suminvaddIτXtX);
      μhat = σ2hat/s.σ2_η[k]*sum(invaddIτXtX*ηk);
      s.μ[k] = randn()*sqrt(σ2hat) + μhat;

      s.β[:,k] = ΣβX*(s.η[k,:] .- s.μ[k]) + sqrt(s.σ2_η[k]).*Lβ*randn(p);
    end
    s.τ_μ = sample_variance(s.μ./s.σ2_η,1.0,hy[:ν0_τ],hy[:τ0_τ]);

    ## save samples
    if t ∈ saveiter
      j = findin(saveiter,t)[1];
      samples[j] = deepcopy(s);
      tss[:,j] = deepcopy(topic);
    end
  end

  return TLMMfit(samples,tss,pss0,hy);
end
