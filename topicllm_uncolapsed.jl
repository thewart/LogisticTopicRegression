function topiclmm_uncollapsed{U<:PostPredSS,D<:Sampleable}(y,X,docrng,pss0::Vector{U},
    topicparam::Vector{Vector{D}},K;
    hy::Dict{Symbol,Float64}=hyperparameter(),
    init::TLMMsample=init_params(K,length(docrng),size(X,1)),
    iter=1000,thin=1)

    ## initialize
    #Base.Test.@test maximum(pss0.span[length(pss0)]) == size(y[1])[1];
    n = length(docrng);
    nw = size(y,2);
    #Base.Test.@test size(X)[2] == n;
    p = size(X,1);
    nd = length.(docrng);

    XtX = X'X;
    Lβ = inv( chol( X*X' + I*inv(hy[:τ_β])));
    ΣβX = Lβ*Lβ'*X;
    invaddIτXtX = inv(I+hy[:τ_β]*XtX);
    suminvaddIτXtX = sum(invaddIτXtX);

    saveiter = thin:thin:iter;
    nsave = length(saveiter);
    iter = maximum(saveiter);

    topicparamout = Array{Vector{D},2}(K,nsave);
    samples = Vector{TLMMsample}(nsave);
    tss = Matrix{Vector{U}}(K,nsave);
    s = init;

    for t in 1:iter
        iΣ = makeprec(XtX,hy[:τ_β],s.τ_μ);

        ##  sample topic memberships
        for i in 1:n
            s.z[docrng[i]] = sample_z(y[:,docrng[i]],topicparam,softmax(s.η[:,i]));
        end
        nk = countz(s.z,docrng,K);

        ##reconstruct topic
        topic = [deepcopy(pss0) for k in 1:K];
        for w in 1:nw addsample!(topic[s.z[w]],y[:,w]); end
        for k in 1:K
            ##sample new topic parameters -- will not work for DirMultPosterior!
            topicparam[k] = randtopic.(topic[k]);

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
            topicparamout[:,j] = topicparam;
        end

    end
    return TLMMfit(samples,tss,pss0,hy), topicparamout
end

function topiclmm_uncollapsed{U<:PostPredSS,D<:Sampleable}(y,docrng,pss0::Vector{U},
    topicparam::Vector{Vector{D}},K;
    hy::Dict{Symbol,Float64}=hyperparameter(),
    init::TLMMsample=init_params(K,length(docrng),size(X,1)),
    iter=1000,thin=1)

    ## initialize
    #Base.Test.@test maximum(pss0.span[length(pss0)]) == size(y[1])[1];
    n = length(docrng);
    nw = size(y,2);
    #Base.Test.@test size(X)[2] == n;
    nd = length.(docrng);

    saveiter = thin:thin:iter;
    nsave = length(saveiter);
    iter = maximum(saveiter);

    topicparamout = Array{Vector{D},2}(K,nsave);
    samples = Vector{TLMMsample}(nsave);
    tss = Matrix{Vector{U}}(K,nsave);
    s = init;

    nk = countz(s.z,docrng,K);
    for t in 1:iter

        ##  sample topic memberships
        if K>1
            for i in 1:n
                s.z[docrng[i]] = sample_z(y[:,docrng[i]],topicparam,softmax(s.η[:,i]));
            end
            nk = countz(s.z,docrng,K);
        end

        ##reconstruct topic
        topic = [deepcopy(pss0) for k in 1:K];
        for w in 1:nw addsample!(topic[s.z[w]],y[:,w]); end
        for k in 1:K
            ##sample new topic parameters -- will not work for DirMultPosterior!
            topicparam[k] = randtopic.(topic[k]);

            if K>1
                ## sample η and λ
                ηk = sample_η(s.η[k,:],s.μ[k],s.σ2_η[k],s.η[vcat(1:(k-1),(k+1):end),:],nd,nk[k,:]);
                s.η[k,:] = ηk;

                σ2hat = inv(inv(s.τ_μ) + n/s.σ2_η[k]);
                μhat = σ2hat/s.σ2_η[k]*sum(ηk);
                s.μ[k] = randn()*sqrt(σ2hat) + μhat;

                s.σ2_η[k] = sample_variance(ηk-s.μ[k],1.0,hy[:ν0_σ2η],hy[:σ0_σ2η]);
            end
        end

        ## save samples
        if t ∈ saveiter
            j = findin(saveiter,t)[1];
            samples[j] = deepcopy(s);
            tss[:,j] = deepcopy(topic);
            topicparamout[:,j] = topicparam;
        end

    end
    return TLMMfit(samples,tss,pss0,hy), topicparamout
end
