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
library(MASS)
library(ape)

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

releves.sb <- releves[grep('[john]',releves$plot),]
releves.sc <- releves[grep('[start]',releves$plot),]
header.sb <- header[header$observer == 'Sobotik et al.',]
header.sc <- header[header$observer == 'Blaschka',]

analysis.sb <- data.frame(header.sb$plot, header.sb$RFA, header.sb$NDF, header.sb$ADF, header.sb$ADL, header.sb$msm)
names(analysis.sb) <- c('plot', 'RFA', 'NDF', 'ADF', 'ADL', 'msm')

# Analyse Daten aus dem Almprojekt ----

# Transformation, einfach...

#releves.sb <- releves.sb[releves.sb$cov !='r',]

releves.sb$cov <- gsub('[r]', 0.1, releves.sb$cov)
releves.sb$cov <- gsub('[+]', 0.5, releves.sb$cov)

# cov ist Text, umwandlung zu Zahl (Int)
releves.sb$cov <- as.integer(releves.sb$cov)

#Kreuztabelle
releves.sb.wide <- cast(releves.sb, plot~abbr, value='cov')
releves.sb.wide[is.na(releves.sb.wide)] <- 0

# Artkürzel...

cap.sb <- capscale(releves.sb.wide ~ RFA + NDF + ADF + ADL, header.sb)
cap.sb
anova(cap.sb, by='term')
anova(cap.sb, by='margin')
anova(cap.sb, by='axis')

stems <- colSums(releves.sb.wide)

plot(cap.sb)
plot(cap.sb, dis="sites")
ordilabel(cap.sb, dis='sp', lab=names(releves.sb.wide), priority = stems)

ordisurf(cap.sb, header.sb$RFA, add = TRUE, col = "green")
ordisurf(cap.sb, header.sb$msm, add = TRUE, col = "gray")

analysis.sb$plot <- NULL
ef <- envfit(cap.sb, analysis.sb, permu = 999)
ef
plot(cap.sb, display="sites", main="Ordination Aufnahmen Almprojekt nach RFA, NDF, ADF und ADL",
	 						  sub="Grüne Linien: RFA; graue: Seehöhe")
plot(ef)
ordisurf(cap.sb, header.sb$RFA, add = TRUE, col = "green")
ordisurf(cap.sb, header.sb$msm, add = TRUE, col = "gray")

rankindex(scale(header.sb$msm), species.sb.wide, c("euc","man","bray","jac","kul"))
dis <- vegdist(species.sb.wide, method = 'jaccard')
adonis(dis ~ geology + msm, header.sb)
