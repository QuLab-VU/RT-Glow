# Wandishin dose-response analysis
# 2022-01-04

library(diprate)

# New function  to obtain drc fits from Clayton's data structure
dipDRC.fromrates <- function (dtf, print.dip = FALSE, norm = FALSE, 
    var = c("cell.line", "drug1", "drug1.conc", "dip_rate","DIPrNorm"),
    plotIt = TRUE, toFile = FALSE, fct = LL.4(), 
    uidName = "uid", 
    ...) 
{
    if (plotIt & toFile) 
        pdf("dipDRC_graph.pdf")
    concName <- var[grep("[Cc]onc", var)]
    Uconc <- sort(unique(dtf[, concName]), decreasing = TRUE)
    dip.rates <- dtf
    
    if (norm) {
        f <- formula(paste0("DIPrNorm ~ ", concName))
        out <- tryCatch({
            drm(f, data = dip.rates, fct = fct)
        }, error = function(cond) {
            return(dip.rates)
        })
    }
    else {
        f <- formula(paste0("dip_rate ~ ", concName))
        out <- tryCatch({
            drc::drm(f, data = dip.rates, fct = fct)
        }, error = function(cond) {
            return(dip.rates)
        })
    }
    if (plotIt & class(out) == "drc") {
        plot.dipDRC(out, ...)
        abline(h = 0, col = grey(0.5), lty = 2)
    }
    if (plotIt & class(out) != "drc") {
        temp <- out
        temp[temp$drug1.conc == 0, "drug1.conc"] <- min(temp[temp$drug1.conc != 
            0, "drug1.conc"])/10
        myargs <- list(...)
        myargs["type"] <- NULL
        myargs <- c(list(formula("dip_rate ~ log10(drug1.conc)"), 
            data = temp), myargs)
        do.call(plot, myargs)
        abline(h = 0, col = grey(0.5), lty = 2)
        legend("bottomleft", "No DRC fit", bty = "n")
    }
    if (plotIt & toFile) 
        dev.off()
    invisible(out)
}



# Load data
d <- read.csv("../data/KScomparisonDF.csv")

# remove column of unknown information (row numbers?)
d <- d[,-1]

# make column of unique conditions
d$ucond <- paste(d$cell.line,d$drug1,d$Count_Type, sep="_")

# save unique unique conditions
uc <- unique(d$ucond)

# try to get drc fits using LL.4 function
z <- lapply(uc[grepl("Direct",uc)], tryCatch({function(u) {
    a <- d[d$ucond==u,]
    a_ctrl <- a[a$drug1.conc==0,]
    b <- rbind(a[a$drug1.conc!=0,], a_ctrl[which(a_ctrl$dip_rate == sample(a_ctrl$dip_rate, 10)),]) 

    dipDRC.fromrates(dtf=b, main=unique(b$ucond), ylim=c(-0.04,0.04))}}, 
        error={function(cond){
            message(cond)
            return(NA)}
            }
))

# check control dip rates
dip_rate_ctrl <- sapply(uc, function(u) mean(d[d$ucond==u & d$drug1.conc==0,"dip_rate"]))

# try to get drc fits using LL.3u function
z2 <- lapply(uc[grepl("Direct",uc)], tryCatch({function(u) {
    a <- d[d$ucond==u,]
    a_ctrl <- a[a$drug1.conc==0,]
    b <- rbind(a[a$drug1.conc!=0,], a_ctrl[which(a_ctrl$dip_rate == sample(a_ctrl$dip_rate, 10)),]) 

    dipDRC.fromrates(dtf=b, main=unique(b$ucond), ylim=c(-1,1), fct=LL.3u(), norm=TRUE, type="all")}}, 
        error={function(cond){
            message(cond)
            return(NA)}
            }
))

names(z2) <- uc[grepl("Direct",uc)]

ci_LL3_fits <- lapply(z2[sapply(z2, class)=="drc"], confint, USE.NAMES = TRUE)
