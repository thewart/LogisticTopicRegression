function topiclmm{U<:PostPredSS,D<:Sampleable}(y,X,pss0::Vector{U},
                           topicparam::Array{D,2},K,docrng,
                           hy::Dict{Symbol,Float64}=hyperparameter();
                           init::TLMMsample=init_params(K,length(docrng),size(X,1)),
                           iter=1000,thin=1)

  ## initialize
  #Base.Test.@test maximum(pss0.span[length(pss0)]) == size(y[1])[1];
  n = length(docrng);
  nw = size(y,2);
  #Base.Test.@test size(X)[2] == n;
  p = size(X,1);
  nd = length.(docrng);

  if !isdefined(init,:z)
      init.z = Vector{Int64}(nw);
  end
  nk = Matrix{Int64}(K,n);
  topic = Vector{Vector{U}}(K);

  XtX = X'X;
  Lβ = inv( chol( X*X' + I*inv(hy[:τ_β])));
  ΣβX = Lβ*Lβ'*X;
  invaddIτXtX = inv(I+hy[:τ_β]*XtX);
  suminvaddIτXtX = sum(invaddIτXtX);

  saveiter = thin:thin:iter;
  nsave = length(saveiter);
  iter = maximum(saveiter);

  topicparamout = Array{D,3}(length(pss0),K,nsave);
  samples = Vector{TLMMsample}(nsave);
  tss = Matrix{Vector{U}}(K,nsave);
  s = init;

  for t in 1:iter
    iΣ = makecov(XtX,hy[:τ_β],s.τ_μ);

    ## reset topic
    map!(k -> deepcopy(pss0),topic,1:K);

    ##  sample topic memberships
    for i in 1:n
        logπ = log.(softmax(s.η[:,i]));
        s.z[docrng[i]] = sample_z(topicparam,logπ,y[:,docrng[i]]);
        nk[:,i] = counts(s.z[docrng[i]],1:K);
    end
    ##reconstruct topic
    for w in 1:nw addsample!(topic[s.z[w]],y[:,w]); end


    for k in 1:K
        ##sample new topic parameters -- will not work for DirMultPosterior!
        topicparam[:,k] = randtopic(topic[k]);

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
          topicparamout[:,:,j] = topicparam;
    end

  end
  return TLMMfit(samples,tss,pss0,hy), topicparamout
end

function sample_z{D<:Sampleable,V<:AbstractVector}(topic::Array{D,2},logπ,y::V)
    K = length(logπ);
    logpost = Vector{Float64}(K);
    for k in 1:K
        logpost[k] = logπ[k];
        for j in 1:length(y) logpost[k] += logpdf(topic[j,k],y[j]); end
    end
    logpostnorm = logpost - logsumexp(logpost);
    return rand(Categorical(exp.(logpostnorm)));
end

function sample_z{D<:Sampleable,M<:AbstractMatrix}(topic::Array{D,2},logπ,y::M)
    return [sample_z(topic,logπ,y[:,i]) for i in 1:size(y,2)];
end
