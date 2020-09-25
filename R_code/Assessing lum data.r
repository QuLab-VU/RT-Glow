library(diprate)
w <- read.csv("../data/H1048 Luminescence Raw/H1048_SNS-314_RawLum.csv", 
                as.is=TRUE, row.names=1)

library(ggplot2)

colnames(w)[1] <- "time"

tall <- data.frame( time=rep(w$time), 
                    counts=as.integer(as.matrix(w[,-1])), 
                    dc=rep(gsub("^X","",colnames(w)[-1]), each=length(w$time)))
tall$rep <- sapply(tall$dc, function(x) {
    z <- as.character(x); as.integer(substr(z,nchar(z),nchar(z)))})

tall$dc <- as.character(tall$dc)
tall$du <- sapply(tall$dc, function(x) ifelse(grepl("u",x),'µM','nM'))

tall$dc <- as.numeric(gsub("DM","0",substr(tall$dc,1,nchar(tall$dc)-3)))

tall$drug1.conc <- tall$dc * sapply(tall$du, function(x) switch(x, µM=1e-6, nM=1e-9, 1))
tall$uc <- paste0('time_',match(tall$time, sort(unique(tall$time))),'_',tall$drug1.conc)

z <- tall
z$dc <- z$du <- NULL
colnames(z)[colnames(z) == 'counts'] <- 'cell.count'

z <- z[order(z$drug1.conc,z$time),]
rownames(z) <- NULL
z$uid <- paste(z$drug1.conc,z$rep,sep="_")

plotGC(z$time, z$cell.count, z$drug1.conc, z$uid)

a <- sapply(unique(tall$uc), function(cond) sum(tall[tall$uc==cond,"counts"]))

s <- tall[!duplicated(tall$uc),]
s$cell.count <- a
s$dc <- s$rep <- s$du <- s$uc <- s$counts <- NULL

do.call(plotGC, getGCargs(s, dat.col=c("time","cell.count","drug1.conc")))


