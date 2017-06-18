K = c(5,10,20,40,80) #from .txt file
r2dat <- data.table(r2=scan("~/analysis/logtopreg/simulation/topiccor_Kseq.csv"),
                    K=ordered(rep(K,K)),type="State content")

n <- 200 #from .txt file
r2dat <- rbind(r2dat,data.table(r2=scan("~/analysis/logtopreg/simulation/pr2_Kseq.csv"),
                   K=rep(K,each=n) %>% ordered,type="Individual phenotypes"))
ggplot(r2dat,aes(y=r2,x=K)) + geom_boxplot() + facet_grid( ~ type,scale="free") + ylab("Correlation") + theme_light()
