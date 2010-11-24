library(mcprofile)
str(cta)
boxplot(foci ~ conc, cta, xlab="concentration", col="lightgreen")
## change class of cta$concentration into factor
cta$concf <- factor(cta$conc, levels=unique(cta$conc))

# glm fit assuming a Poisson distribution for foci counts
# parameter estimation on the log link
# estimating marginal means by removing the intercept
fm <- glm.nb(foci ~ concf-1, data=cta)

### Comparing each dose to the control by Dunnett-type comparisons
# Constructing contrast matrix
CM <- contrMat(table(cta$concf), type="Williams")

# calculating signed root deviance profiles
(dmcp <- mcpcalc(fm, CM))
(dmcpc <- mcpcalc(fm, CM, method="BFGS"))
(dmcpc <- mcpcalc(fm, CM, method="constrOptim"))
glht(fm, linfct=CM)


plot(dmcp)
plot(dmcpc)


hb <- hoa(dmcp)
hbm <- hoa(dmcpc)
plot(hbm, ylim=c(-5,5))
plot(hb, ylim=c(-5,5))


wb <- wald(dmcp)
wbm <- wald(dmcpc)
plot(wbm, ylim=c(-5,5))
plot(wb, ylim=c(-5,5))

confint(dmcp, adjust="single-step")
confint(dmcpc, adjust="single-step")
confint(glht(fm, linfct=CM))

# plot profiles
plot(dmcp)

# multiplicity adjusted p-values
(adjpv <- test(dmcpc,adjust="single-step",alternative="two.sided",margin=1.8))

plot(adjpv, alpha=0.05, order=FALSE)
plot(dmcp, adjpv)

# simultaneous confidence intervals
(ci <- confint(dmcp, adjust="single-step"))
# exponent of confidence limits --> ratio of group means
exp(ci)

plot(exp(ci))
abline(v=1, lty=2)

plot(dmcp, ci)



data("alzheimer", package = "coin")
y <- factor(alzheimer$disease == "Alzheimer", labels = c("other", "Alzheimer"))
gmod <- glm(y ~ smoking:gender-1, data = alzheimer,family = binomial())
summary(gmod)

K <- diag(8)
rownames(K) <- c("None:Female","<10:Female","10-20:Female",">20:Female",
                 "None:Male","<10:Male","10-20:Male",">20:Male")

gm <- mcpcalc(gmod, K)
ci <- confint(gm, adjust="single-step")
ci@confint <- data.frame(binomial()$linkinv(as.matrix(ci@confint)))
ci@estimate <- binomial()$linkinv(ci@estimate)

plot(ci)


gmod2 <- glm(y ~ smoking + gender-1, data = alzheimer,family = binomial())
K2 <- cbind(contrMat(table(alzheimer$smoking), type="Dunnett"), gender=0)
gm2 <- mcpcalc(gmod2, K2)
ci2 <- confint(gm2, adjust="single-step")

par(mar=c(5,8,2,2))
plot(exp(ci2))
box()
abline(v=1, lty=2)

gg <- glht(gmod2, linfct=mcp(smoking="Dunnett"))


data("babies", package = "cond")
fit <- glm(cbind(r1,r2) ~ lull-1 + as.factor(day), data=babies, family=binomial())
K <- rbind(c(-1,1,rep(0,17)))
bm <- mcpcalc(fit, K)
wbm <- wald(bm)
hbm <- hoa(bm)


fit <- glm(cbind(r1,r2) ~ lull +as.numeric(day), data=babies, family=binomial())
K <- rbind(c(-1,1,0))
bm <- mcpcalc(fit, K)


###
nt <- 5
nb <- 20
Ns <- 6
mt <- rbeta(nt, 1,1)
mb <- rbeta(nb, 1,1)
m <- rep(mb, nt) * rep(mt, each=nb)

trt <- as.factor(rep(LETTERS[1:nt], each=nb))
bl <- as.factor(rep(letters[1:nb], nt))
resp <- rbinom(length(m), size=Ns, m)
dat <- data.frame(s=resp, f=Ns-resp, trt, bl)

fit <- glm(cbind(f,s) ~ 0 + trt + bl, data=dat, family=binomial(link="logit"))

CMt <- contrMat(table(dat$trt), type="AVE")
CMb <- matrix(0, nrow=nrow(CMt), ncol=nb-1)
CM <- cbind(CMt, CMb)

mm <- mcpcalc(fit, CM)
ci <- confint(mm, adjust="single-step")
hmm <- hoa(mm)
cih <- confint(hmm, adjust="single-step")
wmm <- wald(mm)
ciw <- confint(wmm, adjust="single-step")


plot(hmm)

