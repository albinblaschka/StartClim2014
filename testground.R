



# Beispiel aus Paket-Publikation ----
data(Tasmania)
attach(Tasmania)
tasmvabund <- mvabund(copepods)
plot(tasmvabund ~ treatment, col= as.numeric(block))

# Eigene Tests und Entwürfe ----
library(mvabund)
library(reshape)
library(vegan)

sobotik <- read.csv2('https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/species.csv',
					 header=TRUE, sep=';',stringsAsFactors = FALSE)

sobotik.head <- read.csv2('https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/sites_wide.csv',
						 header=TRUE, sep=';')

sobotik$voucher <- sobotik$history <- sobotik$comment <- sobotik$layer <- NULL
names(sobotik) <- c('plot','species','cov')

sobotik$cov <- gsub("+","0.1", sobotik$cov, fixed = TRUE)
sobotik$cov <- gsub("r","0", sobotik$cov, fixed = TRUE)
sobotik$cov <- as.numeric(sobotik$cov)
vegtab <- cast(sobotik, plot~species, value='cov', fun.aggregate='sum')
names(vegtab) <- make.cepnames(names(vegtab))
rownames(vegtab) <- vegtab$plot
vegtab$plot <- NULL
vegtab[] <- lapply(vegtab,function(x) replace(x, is.na(x), 0))


# https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf

ord <- metaMDS(vegtab)
plot(ord, disp="sites", type="n")
ordihull(ord, sobotik.head$geology, col="blue")
points(ord, disp="sites", pch=as.numeric(sobotik.head$geology))
ord.fit <- envfit(ord ~ msm * geology + year_releve, data = sobotik.head, perm=999)
ord.fit

plot(ord, dis="site")
plot(ord.fit)
ordisurf(ord, sobotik.head$year_releve, add=TRUE)

ord1 <- capscale(vegtab ~ sobotik.head$msm + sobotik.head$geology)
ord1
plot(ord1, display="site")
anova(ord1, by="axis")

# https://cran.r-project.org/web/packages/vegan/vignettes/diversity-vegan.pdf

sor <- betadiver(vegtab, 'sor')
mod <- with(sobotik.head, betadisper(sor, msm))
mod
plot(mod)
boxplot(mod)
ordipointlabel(mod, display="sites", add=TRUE)

smo<- beals(vegtab)
j <- which(colnames(vegtab) == 'Vaccmyrt')
plot(beals(vegtab, species=j, include=FALSE), vegtab[,j],ylab="Occurrence", main="Vaccmyrt",
	 xlab="Probability of occurrence")

# http://cc.oulu.fi/~jarioksa/opetus/metodi/sessio3.pdf

library(labdsv)

importance(vegtab, sobotik.head$geology)
indi <- indval(vegtab, sobotik.head$msm)
summary(indi)
summary(indi, type = "long")

# mvabund ----
geology <- as.list(sobotik.head$geology)
loc <- as.factor(sobotik.head$location)
msl <- as.factor(sobotik.head$msm)

smvabund <- mvabund(vegtab)
plot(smvabund ~ msl, col= as.numeric(geology))
boxplot.mvabund(smvabund)

vv <- as.matrix(vegtab)
mm <- manyglm(vv ~ sobotik.head$msm * sobotik.head$geology, family="negative.binomial")
plot.manyglm(mm)
