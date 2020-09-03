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


d2 <- read.table('../data/cleaned_dynamic_data.txt', as.is=TRUE, header=TRUE)

library(ggplot2)

ggplot(d2, aes(tdata, y = log2(value), color = variable)) + 
    geom_line(aes(y = log2(DMSO_mean), col = "DMSO_mean")) + 
    geom_line(aes(y = log2(X10uM_mean), col = "X10uM_mean")) +
    geom_line(aes(y = log2(cleaned_2.5uM_mean), col = "cleaned_2.5uM_mean")) +
    geom_line(aes(y = log2(cleaned_625nM_mean), col = "cleaned_625nM_mean")) +
    geom_line(aes(y = log2(cleaned_156.25nM_mean), col = "cleaned_156.25nM_mean")) +
    geom_line(aes(y = log2(cleaned_39nM_mean), col = "cleaned_39nM_mean")) +
    geom_line(aes(y = log2(cleaned_9.76nM_mean), col = "cleaned_9.76nM_mean")) +
    geom_line(aes(y = log2(cleaned_2.44nM_mean), col = "cleaned_2.44nM_mean"))

