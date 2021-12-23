
# d <- data.frame(lum=c(11,1,2),direct=c(10,2,2), row.names=c("negative Emax", "positive Emax", "no effect"))

getEmaxType <- function(x) {
    if(x < 0) { out <- "negative rate" } else if (x >= 1) {
                out <- "no effect" } else {
                out <- "positive rate"}
    out <- factor(out, levels=c("negative rate","positive rate","no effect"))
    return(out)
}

d <- read.csv("~/Desktop/lum emax.csv")
d$lum_type <- sapply(d$lum.Emax, getEmaxType)
d$direct_type <- sapply(d$direct.Emax, getEmaxType)


caret::confusionMatrix(d$lum_type, reference=d$direct_type)

cor(d$lum.Emax,d$direct.Emax, method="spearman")



r <- read.csv("~/Desktop/IC50comparison.csv")
s <- lm(loglum ~ logdir, data=r)
summary(s)

