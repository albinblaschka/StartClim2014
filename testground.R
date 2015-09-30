
# Taxonomie-Check? ----

#library(taxize)

#res <- gnr_resolve(names = qry$species, data_source_ids = '158')
#res

#correct <- tnrs(query = qry$species, source = "iPlant_TNRS")[ , -c(5:7)]
#correct


# Beispiel aus Paket-Publikation ----
data(Tasmania)
attach(Tasmania)
tasmvabund <- mvabund(copepods)
plot(tasmvabund ~ treatment, col= as.numeric(block))

# Eigene Tests und Entwürfe ----
library(mvabund)
library(reshape)
library(vegan)
library(RODBC) 

sobotik <- read.csv2('https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/species.csv',
					 header=TRUE, sep=';',stringsAsFactors = FALSE)

sobotik.head <- read.csv2('https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/sites_wide.csv',
						 header=TRUE, sep=';')

mdbConnect <- odbcConnectAccess2007('Vegetation_StartClim2014.accdb') 

startclim <- sqlQuery(mdbConnect, "SELECT plot, species, cov FROM tbl_startclim", errors = TRUE)
startclim.head <- sqlQuery(mdbConnect, "SELECT plot, date, year_releve, month_releve, day_releve, geology, exp, hhl, cover,
						                       grasses, forbs, legumes, mosses, rock, location, latitude, longitude, accuracy, msm, observer
						                 FROM Kopfdaten_Startclim", errors = TRUE)

sobotik$voucher <- sobotik$history <- sobotik$comment <- sobotik$layer <- NULL
names(sobotik) <- c('plot','species','cov')

tbljoined <- rbind(sobotik, startclim)
headjoined <- rbind(sobotik.head,startclim.head)

tbljoined <- tbljoined[tbljoined$plot !='john51',] 
tbljoined <- tbljoined[tbljoined$plot !='john52',] 
tbljoined <- tbljoined[tbljoined$plot !='john53',] 
tbljoined <- tbljoined[tbljoined$plot !='john54',] 

headjoined <- headjoined[headjoined$plot !='john51',] 
headjoined <- headjoined[headjoined$plot !='john52',] 
headjoined <- headjoined[headjoined$plot !='john53',] 
headjoined <- headjoined[headjoined$plot !='john54',] 
headjoined <- headjoined[headjoined$plot !='start51',] 

tbljoined$cov <- gsub("+","0.1", tbljoined$cov, fixed = TRUE)
tbljoined$cov <- gsub("r","0", tbljoined$cov, fixed = TRUE)
tbljoined$cov <- as.numeric(tbljoined$cov)

vegtab <- cast(tbljoined, plot~species, value='cov', fun.aggregate='sum')
names(vegtab) <- make.cepnames(names(vegtab))
rownames(vegtab) <- vegtab$plot
vegtab$plot <- NULL
vegtab[] <- lapply(vegtab,function(x) replace(x, is.na(x), 0))


# https://cran.r-project.org/web/packages/vegan/vignettes/intro-vegan.pdf

ord <- metaMDS(vegtab)
plot(ord, dis="sites", type="n")
ordihull(ord, headjoined$geology, col="blue")

points(ord, disp="sites", pch=as.numeric(headjoined$geology))
ord.fit <- envfit(ord ~ msm * geology , data = headjoined, perm=999)
ord.fit

plot(ord, dis="site")
plot(ord.fit)
ordisurf(ord, headjoined$year_releve, add=TRUE)

ord1 <- capscale(vegtab ~ headjoined$msm + headjoined$geology)
ord1
plot(ord1, display="site")
anova(ord1, by="axis")

# https://cran.r-project.org/web/packages/vegan/vignettes/diversity-vegan.pdf

sor <- betadiver(vegtab, 'sor')
mod <- with(headjoined, betadisper(sor, msm))
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

importance(vegtab, headjoined$geology)
indi <- indval(vegtab, headjoined$msm)
summary(indi)
summary(indi, type = "long")

# mvabund ----
geology <- as.list(headjoined$geology)
loc <- as.factor(headjoined$location)
msl <- as.factor(headjoined$msm)

smvabund <- mvabund(vegtab)
plot(smvabund ~ msl, col= as.numeric(geology))
boxplot.mvabund(smvabund)

vv <- as.matrix(vegtab)
mm <- manyglm(vv ~ headjoined$msm * headjoined$geology, family="negative.binomial")
plot.manyglm(mm)
