
# setwd("~/Data/StartClim2014")
setwd("Y:/Institut3/AGRAM/StartClim2014/Albin/analyse")
library(reshape)
library(vegan)
library(MASS)
library(ape)

species.sb <- read.csv2("https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/species.csv",
				  header = TRUE, sep = ";", dec=',')
header.sb <- read.csv2("https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/sites_wide.csv",
					header = TRUE, sep = ";", dec='.')

species.sc <- read.csv2("vegetation2015.csv",header = TRUE, sep = ";", dec=',')
header.sc <- read.csv2("kopf2015.csv",header = TRUE, sep = ";", dec=',')
header.sc$ID <- NULL

species.sb$layer <- species.sb$voucher <- species.sb$comment <- species.sb$history <- NULL
names(species.sb) <- c('plot', 'species', 'cov')


# Analyse Daten aus dem Almprojekt ----

# Transformation, einfach...

species.sb <- species.sb[species.sb$cov !='r',]

species.sb <- species.sb[species.sb$plot !='john51',] 
species.sb <- species.sb[species.sb$plot !='john52',] 
species.sb <- species.sb[species.sb$plot !='john53',] 
species.sb <- species.sb[species.sb$plot !='john54',] 

header.sb <- header.sb[header.sb$plot !='john51',] 
header.sb <- header.sb[header.sb$plot !='john52',] 
header.sb <- header.sb[header.sb$plot !='john53',] 
header.sb <- header.sb[header.sb$plot !='john54',] 
header.sc <- header.sc[header.sc$plot !='start51',] 

#species.sb$cov <- gsub('[r]', 0.1, species.sb$cov)
species.sb$cov <- gsub('[+]', 1, species.sb$cov)

# cov ist Text, umwandlung zu Zahl (Int)
species.sb$cov <- as.integer(species.sb$cov)

#Kreuztabelle
species.sb.wide <- cast(species.sb, plot~species, value='cov')
species.sb.wide[is.na(species.sb.wide)] <- 0
species.sb.wide$plot <- as.numeric(species.sb.wide$plot)
# Artkürzel...
names(species.sb.wide) <- make.cepnames(names(species.sb.wide))

dis.sb <- vegdist(species.sb.wide)
mds0.sb <- metaMDS(dis.sb)
stressplot(mds0.sb, dis.sb)

ef <- envfit(mds0.sb ~ msm, header.sb, permu = 999)
ef

ordiplot(mds0.sb, type = "t", cex=0.6)
#ordihull(mds0.sb, header.sb$geology, col="darkred")
ordisurf(mds0.sb, header.sb$msm, add = TRUE, col = "dodgerblue")
ordisurf(mds0.sc, header.sc$msm, add = TRUE, col = "green2")


clus.sb <- hclust(dis.sb, 'single')
geoclus.sb <- as.numeric(header.sb$geology)
msmclus.sb <- header.sb$msm
names(geoclus.sb) <- names(msmclus.sb) <- header.sb$plot

x <- scores(mds0.sb, display = "sites", choices = 1)
dendro.sb <- reorder(clus.sb, x)

plot(as.phylo(dendro.sb), cex = 0.5, label.offset = 0.01, tip.color = msmclus.sb/100,
	 main='Clusteranalyse Höhenprofil Johnsbach 1993-1996', sub = 'Quelle: Sobotik et al. 1998')
# plot(clus.sb, cex=0.6)

# Analyse Daten StartClim 2015 ----

# Transformation, einfach...

species.sc <- species.sc[species.sc$cov !='r',]

#species.sc$cov <- gsub('[r]', 0.1, species.sc$cov)
species.sc$cov <- gsub('[+]', 1, species.sc$cov)

# cov ist Text, umwandlung zu Zahl (Int)
species.sc$cov <- as.integer(species.sc$cov)

#Kreuztabelle
species.sc.wide <- cast(species.sc, plot~species, value='cov')
species.sc.wide[is.na(species.sc.wide)] <- 0
species.sc.wide$plot <- as.numeric(species.sc.wide$plot)

# Artkürzel...
names(species.sc.wide) <- make.cepnames(names(species.sc.wide))

dis.sc <- vegdist(species.sc.wide)
mds0.sc <- metaMDS(dis.sc)
stressplot(mds0.sc, dis.sc)
ordiplot(mds0.sc, type = "t", cex=0.6, main = 'StartClim 2015')
ordisurf(mds0.sc, header.sc$msm, add = TRUE, col = "dodgerblue")
ordisurf(mds0.sb, header.sb$msm, add = TRUE, col = "green")

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
species.wide[is.na(species.wide)] <- 0
species.wide$plot <- as.numeric(species.wide$plot)
names(species.wide) <- make.cepnames(names(species.wide))


rankindex(scale(header$year_releve), species.wide, c("euc","man","bray","jac","kul"))

dis <- vegdist(species.wide, method = 'jaccard')
mds0 <- metaMDS(dis)
mds0
stressplot(mds0, dis, main='Almprojekt + StartClim, ohne Bucheck')
ordiplot(mds0, type = "t", cex=0.6, main='Almprojekt + StartClim, ohne Bucheck')

ordisurf(mds0, header$msm, add = TRUE, col = "dodgerblue")
ordisurf(mds0.sb, header.sb$msm, add = TRUE, col = "green")

ef <- envfit(mds0 ~ year_releve + geology, header, permu = 999)
ef
#plot(ef)

clus <- hclust(dis, 'single')

geoclus <- as.numeric(header$geology)
msmclus <- header$msm
names(geoclus) <- names(msmclus) <- header$plot

plot(as.phylo(clus), cex = 0.5, label.offset = 0.01, tip.color = msmclus/100, main='Almprojekt + StartClim, ohne Bucheck')
# plot(dendro, cex=0.6, main='Almprojekt + StartClim, ohne Bucheck')


