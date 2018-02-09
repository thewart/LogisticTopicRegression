function sample_z!(z,topic,nk,π,y,K)

  logpost = Vector{Float64}(K);
  for j in 1:length(z)
    zj = z[j];
    pullsample!(topic[zj],y[:,j]);
    nk[zj] -= 1;

    for k in 1:K
      logpost[k] = log(π[k]) + lppd(y[:,j],topic[k]);
    end
    logpostnorm = logpost - logsumexp(logpost);
    z[j] = rand(Categorical(exp.(logpostnorm)));

    addsample!(topic[z[j]],y[:,j]);
    nk[z[j]] += 1;
  end
end

function sample_η(η::Vector{Float64},ηp::Array{Float64,2},
  iΣ::Array{Float64,2},nd::Vector{Int64},nk::Vector{Int64})
  n = length(nd);
  w = Vector{Float64}(n);
  for i in 1:n
    c = logsumexp(ηp[:,i]);
    ρ = η[i] - c;
    λ = rpolyagamma(ρ,nd[i]);
    w[i] = nk[i] - nd[i]/2 + c*λ;
    iΣ[i,i] += λ;
  end
  return iΣ \ w + chol(Hermitian(iΣ)) \ randn(n)
end

function sample_variance(y::Vector{Float64},V::Array{Float64,2},ν0::Float64,σ0::Float64)
  a = 0.5(length(y)+ν0);
  b = 0.5(σ0*ν0 + dot(y,V*y));
  return rand(InverseGamma(a,b));
end

function sample_variance(y::Vector{Float64},σ2::Float64,ν0::Float64,σ0::Float64)
  a = 0.5(length(y)+ν0);
  b = 0.5(σ0*ν0 + dot(y,y)/σ2);
  return rand(InverseGamma(a,b));
end
