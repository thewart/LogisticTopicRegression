function makecov(XtX::Array{Float64,2},ZtZ::Array{Float64,2},
  τ_β::Float64,τ_u::Float64,τ_μ::Float64)

  n = size(XtX)[1];
  V = Array{Float64,2}(n,n);
  for i in 1:n
    for j in 1:n
      V[j,i] = XtX[j,i].*τ_β + ZtZ[j,i].*τ_u + (i==j ? 1+τ_μ : τ_μ);
    end
  end

  return inv(V)
end

function makecov(XtX::Array{Float64,2},τ_β::Float64,τ_μ::Float64)

  n = size(XtX)[1];
  V = Array{Float64,2}(n,n);
  for i in 1:n
    for j in 1:n
      V[j,i] = XtX[j,i].*τ_β + (i==j ? 1+τ_μ : τ_μ);
    end
  end

  return inv(V)
end


function init_post(n::Int64,nd::Vector{Int64},K::Int64,p::Int64,nsave::Int64,isu::Bool)
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
  return post
  if isu
    post[:u] = Array{Float64}(n,K,nsave);
    post[:τ_u] = Array{Float64,2}(K,nsave);
  end
end

function hyperparameter(;ν0_σ2η=1.0,σ0_σ2η = 1.0,
  τ0_τ = 0.25,ν0_τ = 1.0,τ0_u = 0.25,ν0_u=1.0,τ_β=1.0)
  return Dict(:ν0_σ2η => ν0_σ2η, :σ0_σ2η => σ0_σ2η, :τ0_τ => τ0_τ,
  :ν0_τ => ν0_τ, :τ0_u => τ0_u, :ν0_u => ν0_u, :τ_β => τ_β)
end

function refβ(β::Array{Float64,2},refk::Int64)
  return β .- β[:,refk:refk]
end
refβ(β::Array{Float64,3},refk::Int64) = β .- β[:,refk:refk,:]
refβ(β::Array{Float64},μ::Array{Float64,1}) = refβ(β,findmax(μ)[2]);
refβ(β::Array{Float64},μ::Array{Float64,2}) = refβ(β,findmax(mean(μ,2))[2]);

function writefit(fit::Dict{Symbol,Union{AbstractArray,Dict{Symbol,Float64}}},path::String)
  if !isdir(path) mkpath(path); end
  n = length(fit[:z]);
  K,nsave = size(fit[:topic]);
  p = size(fit[:β])[1];
    b = length(fit[:topic][1]);
  writecsv(string(path,"dims.csv"),[n,K,p,b,nsave]);
  writecsv(string(path,"sigma.csv"),hcat(fit[:τ],fit[:σ2]',fit[:τ_u]'));
  writecsv(string(path,"beta.csv"),fit[:β][:]);
  writecsv(string(path,"eta.csv"),fit[:η][:]);
  writecsv(string(path,"u.csv"),fit[:u][:]);
  #writecsv(string(path,"loglik.csv"),fit[:loglik]');
  writecsv(string(path,"mu.csv"),fit[:μ]');
  writecsv(string(path,"z.csv"),vcat(fit[:z]...));
  writecsv(string(path,"nd.csv"),map(x -> size(x)[1],fit[:z]));
  writecsv(string(path,"topicmean.csv"),
            mapreduce(y -> mean(topicppd(y)),vcat,fit[:topic]));
  writecsv(string(path,"topicvar.csv"),
            mapreduce(x -> map(var,topicppd(x)),vcat,fit[:topic]));

  θ = Vector{Float64}(0);
  for i in 1:(length(fit[:topic]))
      append!(θ,vcat(map(x -> params(x)[1],topicppd(fit[:topic][i]))...));
  end
  writecsv(string(path,"topicparams.csv"),θ);
  writecsv(string(path,"paramlengths.csv"),
  map(x -> length(params(x)[1]),topicppd(fit[:topic][1])));
end
