d <- read.csv('../mKate2counts_TotHour.csv', row.names=1, as.is=TRUE)

w <- d[,c("TotHour",colnames(d)[grep("^R",colnames(d))])]
colnames(w)[1] <- "time"

tall <- data.frame( time=rep(w$time), 
                    counts=as.integer(as.matrix(w[,-1])), 
                    well=rep(colnames(w)[-1], each=length(w$time)))
 
s <- apply(d[,grep("^R",colnames(d))], 1, sum)
s <- data.frame(time=d$TotHour,counts=as.integer(s))

plot(counts ~ time, data=tall, type="n")
for(well in unique(tall$well)) lines(counts ~ time, data=tall[tall$well==well,])

m <- lm(counts ~ well * time, data=tall)

