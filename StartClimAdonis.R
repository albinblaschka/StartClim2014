## -----------------------------------------------------------------------------------------------------
##
## Vegetationsökologische Analyse im Rahmen des Projektes
## StartClim2014.X: Zur Bedeutung des Klimawandels für die Ernährung und Krankheiten alpiner Wildarten
## Teilprojekt Änderungen in der Vegetation durch Klimaänderungen
##
## Autoren: Mag. Thomas Guggenberger MSc., Dr.rer.nat. Albin Blaschka, HBLFA Raumberg-Gumpenstein
## Oktober 2015
##
## -----------------------------------------------------------------------------------------------------


## Initialisierung, Einlesen und Aufbereitung der Daten ----

setwd("//storage02/transfer/Institut3/AGRAM/StartClim2014/Albin/analyse")

library(RODBC)
library(reshape)
library(vegan)

mdbConnect <- odbcConnectAccess2007('../../Almprojekt_StartClim2014.accdb') 
releves <- sqlQuery(mdbConnect, "SELECT plot, species, cov FROM Vegetationsaufnahmen", errors = TRUE)
allSpecies <- sqlQuery(mdbConnect, "SELECT species FROM Gesamtartenliste", errors = TRUE)
header <- sqlQuery(mdbConnect, "SELECT site_id, year_releve, rdate, location, geology, expo,
                   msm, plot, cover, grasses, forbs, legumes,
                   mosses, rock, RFA, NDF, ADF, ADL, observer
                   FROM Kopfdaten_Analysen")
odbcClose(mdbConnect)

allSpecies$abbr <- make.cepnames(allSpecies$species)
releves <- merge(releves,allSpecies, by='species')
releves$species <- NULL

# Anpassen der Daten, Transformation... ----
releves$cov <- gsub('[r]', 0.1, releves$cov)
releves$cov <- gsub('[+]', 0.5, releves$cov)
releves$cov <- as.numeric(releves$cov)


releves.sb <- releves[grep('[john]',releves$plot),]
releves.sc <- releves[grep('[start]',releves$plot),]
header.sb <- header[header$observer == 'Sobotik et al.',]
header.sc <- header[header$observer == 'Blaschka',]

releves.wide <- cast(releves, plot~abbr, value='cov')
releves.wide[is.na(releves.wide)] <- 0
row.names(releves.wide) <- releves.wide$plot
releves.wide$plot <- NULL

releves.sb.wide <- cast(releves.sb, plot~abbr, value='cov')
releves.sb.wide[is.na(releves.sb.wide)] <- 0
row.names(releves.sb.wide) <- releves.sb.wide$plot
releves.sb.wide$plot <- NULL

releves.sc.wide <- cast(releves.sc, plot~abbr, value='cov')
releves.sc.wide[is.na(releves.sc.wide)] <- 0
row.names(releves.sc.wide) <- releves.sc.wide$plot
releves.sc.wide$plot <- NULL

# Plot 'start51' keine Daten, wird in Kopf eliminiert
#header <- header[!is.na(header$RFA),]

dis <- vegdist(releves.wide, method = 'jaccard')
dis.sb <- vegdist(releves.sb.wide, method = 'jaccard')
dis.sc <- vegdist(releves.sc.wide, method = 'jaccard')

x <- data.frame(labels(dis))
names(x) <- c('plot')
hh <- merge(x,header, by='plot')

x.sb <- data.frame(labels(dis.sb))
names(x.sb) <- c('plot')
hh.sb <- merge(x.sb,header.sb, by='plot')

x.sc <- data.frame(labels(dis.sc))
names(x.sc) <- c('plot')
hh.sc <- merge(x.sc,header.sc, by='plot')

# Analyse ----

### Almprojekt 1993 - 1996
##########################

ad.sb <- adonis(dis.sb ~ year_releve * geology * msm, hh.sb, method = 'jaccard', permutations = 999, strata=hh.sb$msm)
ad.sb
be.sb <- betadisper(dis.sb, hh.sb$msm)
be.sb
plot(be.sb, main = "Almprojekt 1993 - 1996, nach Seehöhe")
ordipointlabel(be.sb, display="centroids", add=TRUE)

### StartClim 2015, nach Seehöhe
################################

ad.sc <- adonis(dis.sc ~ msm  * geology, hh.sc, method = 'jaccard', permutations = 999, strata=hh.sc$msm)
ad.sc
be.sc <- betadisper(dis.sc, hh.sc$msm)
be.sc
plot(be.sc, main = "StartClim 2015, nach Seehöhe")
ordipointlabel(be.sc, display="sites", add=TRUE)

### Gesamter Datensatz, nach Seehöhe und Geologie
#################################################

ad <- adonis(dis ~ msm * geology, hh, method = 'jaccard', permutations = 999, strata=hh$msm)
ad
be <- betadisper(dis, hh$msm)
be
plot(be, main = "Gesamter Datensatz, nach Seehöhe")
ordipointlabel(be, display="centroids", add=TRUE)

### Gesamter Datensatz, nach Jahren, Geologie und Seehöhe
#########################################################

ad.year <- adonis(dis ~ year_releve * geology * msm, hh, method = 'jaccard', permutations = 5000, strata=hh$msm)
ad.year
be.year <- betadisper(dis, hh$year_releve)
be.year
plot(be.year, main = "Gesamter Datensatz, nach Jahren")
ordipointlabel(be.year, display="centroids", add=TRUE)

be.sb.year <- betadisper(dis.sb, hh.sb$year_releve)
be.sb.year
plot(be.sb.year, main = "Almprojekt 1993 - 1996, nach Jahren")
ordipointlabel(be.sb.year, display="centroids", add=TRUE)


