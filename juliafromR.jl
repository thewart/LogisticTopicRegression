function getfocaldata()
  R"source('/home/seth/code/LogisticTopicRegression/R2julia.R')";
  X = @> rcopy("X") transpose();

  nc = Int64(rcopy("length(Y)-ncovcols"));
  k = @>> rcopy("Y[,lapply(.SD,
    function(x) length(unique(x))),.SD=-(1:ncovcols)] %>% as.matrix") convert(Array{Int64}) vec();
  pss0 = VectorPosterior(CategoricalPosterior(k[1]),1);
  for i in 2:nc pss0 = vcat(pss0,CategoricalPosterior(k[i])); end

  #pss0 = VectorPosterior(PoissonPosterior(0.1,0.1),np) #,
  #  DirMultPosterior(ns[1]),1);
  #for i in 2:length(ns) pss0 = vcat(pss0,DirMultPosterior(ns[i])); end

  n = size(X)[2];
  Y = Vector{Array{Int64,2}}(n);
  for i in 1:n
    @rput i
    Y[i] = @> rcopy(
    #"Y[obsgroup[[i]],names(Y) %in% c(ptetho,unlist(stetho)),with=F]") Array{Int64}() transpose();
    "Y[obsgroup[[i]],-(1:ncovcols),with=F]") Array{Int64}() transpose();
  end

  return Y,X,pss0
end
