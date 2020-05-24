# Model for static luminescence relationship to cell count

a <- read.csv('../H841_Static_Lum_5k.csv', as.is=TRUE)
a <- data.frame(t(a))

n_cells <- as.integer(gsub("^X","",rownames(a)))
a <- data.frame(lum=c(as.matrix(a)), n_cells=rep(n_cells), measurement=rep(1:4, each=9))
bkg <- mean(a[a$n_cells==0,'lum'])
a$lum <- a$lum - bkg
a <- a[a$n_cells != 0,]


plot(lum ~ n_cells, data=a)
m <- lm(lum ~ n_cells, data=a)
abline(m, col='orange')


# SSasymp(input, Asym, R0, lrc)
# input: a numeric vector of values at which to evaluate the model.
# Asym: a numeric parameter representing the horizontal asymptote on the right side (very large values of input).
# R0: a numeric parameter representing the response when input is zero.
# lrc: a numeric parameter representing the natural logarithm of the rate constant.

# function pulled from SSasymp help description in R
exp_decay <- function(x, ymax, R0, log_alpha) ymax+(R0-ymax)*exp(-exp(log_alpha)*x)

curve(exp_decay(x, coef(fit)['ymax'], coef(fit)['R0'], coef(fit)['log_alpha']), 
      from=0, to=25000, col='purple', lwd=2, add=TRUE) 

dev.new()
plot(log2(lum) ~ log2(n_cells), data=a)
m <- lm(log2(lum) ~ log2(n_cells), data=a)
abline(m, col='green')
m2 <- lm(log2(lum) ~ log2(n_cells), data=a[a$n_cells > 40,])
abline(m2, col='orange')

