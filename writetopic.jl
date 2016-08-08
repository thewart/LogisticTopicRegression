function writefit(fit::Dict{Symbol,AbstractArray},path::ASCIIString)
  if !isdir(path) mkdir(path); end
  n = length(fit[:z]);
  K,nsave = size(fit[:topic]);
  p = size(fit[:β])[1];
  b = size(fit[:topic][1]);
  writecsv(string(path,"dims.csv"),[n,K,p,b,nsave]);
  writecsv(string(path,"sigma.csv"),hcat(fit[:τ]',fit[:σ2]'));
  writecsv(string(path,"beta.csv"),fit[:β][:]);
  writecsv(string(path,"eta.csv"),fit[:η][:]);
  writecsv(string(path,"loglik.csv"),fit[:η]');
  writecsv(string(path,"mu.csv"),fit[:μ]');
  writecsv(string(path,"topicmean.csv"),
            mapreduce(y -> mean(topicppd(y)),vcat,fit[:topic]));
end
