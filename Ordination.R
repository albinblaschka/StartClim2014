

library(reshape)

library(vegan)
library(MASS)

dta <- read.csv2("https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/species.csv",
				  header = TRUE, sep = ";", dec=',')
header <- read.csv2("https://raw.githubusercontent.com/albinblaschka/vegsoup-data/master/sobotik1998/sites_wide.csv",
					header = TRUE, sep = ";", dec='.')

species <- data.frame(dta$plot,dta$abbr,dta$cov)
names(species) <- c('plot','species','cov')

# Transformation, einfach...
species$cov <- gsub('[r]', 0.1, species$cov)
species$cov <- gsub('[+]', 1, species$cov)

# cov ist Text, umwandlung zu Zahl (Int)
species$cov <- as.integer(species$cov)

#Kreuztabelle
species.wide <- cast(species, plot~species, value='cov')

# Artkürzel...
names(species.wide) <- make.cepnames(names(species.wide))
species.wide[] <- lapply(species.wide,function(x) replace(x, is.na(x), 0))

dis <- vegdist(species.wide)
clus <- hclust(dis, 'single')
plot(clus)


