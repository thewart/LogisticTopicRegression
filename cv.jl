using RCall, LogTopReg, Distributions
import Lazy.@>, Lazy.@>>

R"source('/home/seth/code/LogisticTopicRegression/R2julia.R')";
R"assfold <- read.csv('~/analysis/logtopreg/crossvalidation/folds_5k.csv',header = F)[,1]"
R"ID <- Y[,unique(FocalID)]";
nfold = rcopy("max(assfold)");
n = rcopy("Y[,length(ID)]");
r = readcsv("/home/seth/analysis/logtopreg/crossvalidation/optimout_5k.csv")
nc = Int64(rcopy("length(Y)-ncovcols"));
k = @>> rcopy("Y[,lapply(.SD,
  function(x) length(unique(x))),.SD=-(1:ncovcols)] %>% as.matrix") convert(Array{Int64}) vec();
pss0 = VectorPosterior(CategoricalPosterior(k[1]),1);
for i in 2:nc pss0 = vcat(pss0,CategoricalPosterior(k[i])); end


Y1 = Vector{Array{Int64,2}}(n);
Y2 = Vector{Array{Int64,2}}(n);
f = 1;
@rput f
R"Y1 <- Y[assfold!=f]"
R"Y2 <- Y[assfold==f]"
for i in 1:n
  @rput i
  Y1[i] = @> rcopy(R"Y1[FocalID==ID[i],-(1:ncovcols),with=F]") Array{Int64}() transpose();
  Y2[i] = @> rcopy(R"Y2[FocalID==ID[i],-(1:ncovcols),with=F]") Array{Int64}() transpose();
end

zflat = @> mapslices(x -> rand(Categorical(x)),r[r[:,11].==f,1:10],2) vec;
z = Vector{Vector{Int64}}(length(Y1));
nd = map(x -> size(x)[2],Y1);
guh = 1;
for i in 1:n
  z[i] = zflat[guh:(guh+nd[i]-1)];
  guh += nd[i];
end

fit = topiclmm(Y1,randn(1,n),pss0,10,hyperparameter(τ_β=1e-10),zinit=z,iter=100);
lppdout = lppd(Y2,fit);
