
#add one new observation
function addsample!{T<:Real,U<:PostPredSS}(pp::Vector{U},ynew::AbstractVector{T})
    for i in 1:length(pp) addsample!(pp[i],ynew[i]); end
end

addsample!{T<:Real,U<:PostPredSS}(pp::Vector{U},ynew::Array{T,2}) =
for i in 1:length(pp) addsample!(pp,ynew[:,j]); end

#addsample!{T<:Real}(pp::VectorPosterior,ynew::T) =
#  for i in 1:length(pp) addsample!(pp[i],ynew); end

#remove one observation
function pullsample!{T<:Real,U<:PostPredSS}(pp::Vector{U},yold::AbstractVector{T})
    for i in 1:length(pp) pullsample!(pp[i],yold[i]); end
end

pullsample!{T<:Real,U<:PostPredSS}(pp::Vector{U},yold::Array{T,2}) =
for i in 1:length(pp) pullsample!(pp,yold[:,j]); end


#utilities for distribution vectors
function rand{T<:Sampleable}(dv::Vector{T},n::Int64=1)
    p = length(dv);
    X = Array{Float64}(p,n);
    for i in 1:n
        for j in 1:p
            X[j,i] = rand(dv[j]);
        end
    end
    return X
end

function mean{T<:Sampleable}(dv::Vector{T})
    l = map(length,dv);
    p = sum(l);
    X = Vector{Float64}(p);
    k = 0;
    for j in 1:p
        X[(k+1):(k+l[j])] = mean(dv[j]);
        k = k + l[j];
    end
    return X
end

function topicpd{U<:PostPredSS}(pp::Vector{U})
    return topicpd.(pp)
end
