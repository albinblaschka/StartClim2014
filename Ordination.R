
setwd("~/Data/StartClim2014")
library(reshape)
library(vegan)
library(MASS)

species.sb <- read.csv2("https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/species.csv",
				  header = TRUE, sep = ";", dec=',')
header.sb <- read.csv2("https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/sites_wide.csv",
					header = TRUE, sep = ";", dec='.')

species.sc <- read.csv2("vegetation2015.csv",header = TRUE, sep = ";", dec=',')
header.sc <- read.csv2("kopf2015.csv",header = TRUE, sep = ";", dec=',')
header.sc$ID <- NULL

species.sb$layer <- species.sb$voucher <- species.sb$comment <- species.sb$history <- NULL

# Analyse Daten aus dem Almprojekt ----

# Transformation, einfach...

species.sb <- species.sb[species.sb$cov !='r',]

#species.sb$cov <- gsub('[r]', 0.1, species.sb$cov)
species.sb$cov <- gsub('[+]', 1, species.sb$cov)

# cov ist Text, umwandlung zu Zahl (Int)
species.sb$cov <- as.integer(species.sb$cov)

#Kreuztabelle
species.sb.wide <- cast(species.sb, plot~species, value='cov')
species.sb.wide[] <- lapply(species.sb.wide,function(x) replace(x, is.na(x), 0))

# Artkürzel...
names(species.sb.wide) <- make.cepnames(names(species.sb.wide))

dis.sb <- vegdist(species.sb.wide)
mds0.sb <- metaMDS(dis.sb)
stressplot(mds0.sb, dis.sb)
ordiplot(mds0.sb, type = "t", cex=0.6)
ordihull(mds0.sb, header.sb$msm, col="blue")

clus.sb <- hclust(dis.sb, 'single')
plot(clus.sb, cex=0.6)

# Analyse Daten StartClim 2015 ----

# Transformation, einfach...

species.sc <- species.sc[species.sc$cov !='r',]

#species.sc$cov <- gsub('[r]', 0.1, species.sc$cov)
species.sc$cov <- gsub('[+]', 1, species.sc$cov)

# cov ist Text, umwandlung zu Zahl (Int)
species.sc$cov <- as.integer(species.sc$cov)

#Kreuztabelle
species.sc.wide <- cast(species.sc, plot~species, value='cov')
species.sc.wide[] <- lapply(species.sc.wide,function(x) replace(x, is.na(x), 0))

# Artkürzel...
names(species.sc.wide) <- make.cepnames(names(species.sc.wide))

dis.sc <- vegdist(species.sc.wide)
mds0.sc <- metaMDS(dis.sc)
stressplot(mds0.sc, dis.sc)
ordiplot(mds0.sc, type = "t", cex=0.6)

clus.sc <- hclust(dis.sc, 'single')
plot(clus.sc, cex=0.6)

# Analyse gemeinsam ----
species <- rbind(species.sb, species.sc)
header <- rbind(header.sb, header.sc)

species <- species[species$plot !='john51',] 
species <- species[species$plot !='john52',] 
species <- species[species$plot !='john53',] 
species <- species[species$plot !='john54',] 

header <- header[header$plot !='john51',] 
header <- header[header$plot !='john52',] 
header <- header[header$plot !='john53',] 
header <- header[header$plot !='john54',] 
header <- header[header$plot !='start51',] 

species.wide <- cast(species, plot~species, value='cov')
species.wide[] <- lapply(species.wide,function(x) replace(x, is.na(x), 0))
names(species.wide) <- make.cepnames(names(species.wide))


rankindex(scale(header$msm), species.wide, c("euc","man","bray","jac","kul"))

dis <- vegdist(species.wide, method = 'jaccard')
mds0 <- metaMDS(dis)
mds0
stressplot(mds0, dis, main='Almprojekt + StartClim, ohne Bucheck')
ordiplot(mds0, type = "t", cex=0.6, main='Almprojekt + StartClim, ohne Bucheck')
#ordihull(mds0, header$geology, col="blue")
ordihull(mds0, header$msm, col="red")

ef <- envfit(mds0 ~ geology + msm + exp, header, permu = 999)
ef
plot(ef)


clus <- hclust(dis, 'single')

geoclus <- as.numeric(header$geology)
msmclus <- header$msm
names(geoclus) <- names(msmclus) <- header$plot

library(ape)

plot(as.phylo(clus), cex = 0.5, label.offset = 0.01, tip.color = msmclus/100, main='Almprojekt + StartClim, ohne Bucheck')
# plot(dendro, cex=0.6, main='Almprojekt + StartClim, ohne Bucheck')


