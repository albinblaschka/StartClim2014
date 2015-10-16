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
plot(analysis.sb, panel=panel.smooth)


releves$cov <- gsub('[r]', 0.1, releves$cov)
releves$cov <- gsub('[+]', 0.5, releves$cov)

# cov ist Text, umwandlung zu Zahl 
releves$cov <- as.numeric(releves$cov)
releves.wide <- cast(releves, plot~abbr, value='cov')
releves.wide[is.na(releves.wide)] <- 0
row.names(releves.wide) <- releves.wide$plot
releves.wide$plot <- NULL

analysis <- data.frame(header$plot, header$RFA, header$NDF, header$ADF, header$ADL, header$msm)
names(analysis) <- c('plot', 'RFA', 'NDF', 'ADF', 'ADL', 'msm')


# Analyse Gesamt-Datensatz ----

# Plot 'start51' keine Daten, wird in Kopf eliminiert

header <- header[!is.na(header$RFA),]

rankindex(scale(analysis$RFA), releves.wide, c("euc","man","bray","jac","kul"))
dis <- vegdist(releves.wide, method = 'jaccard')
adonis(dis ~  geology, header)

cap <- capscale(releves.wide ~ geology + expo, header)
cap
ef <- envfit(cap, header$RFA * header$year_releve, strata = header$geology)
plot(cap, display="sites")
ordisurf(cap, header$RFA, add = TRUE, col = "green")
#ordisurf(cap, header$msm, add = TRUE, col = "gray")
ordihull(cap, header$geology, col="dodgerblue", label= TRUE )

plot(ef)



# Analyse Daten aus dem Almprojekt ----

# Transformation, einfach...

#releves.sb <- releves.sb[releves.sb$cov !='r',]

releves.sb$cov <- gsub('[r]', 0.1, releves.sb$cov)
releves.sb$cov <- gsub('[+]', 0.5, releves.sb$cov)

# cov ist Text, umwandlung zu Zahl 
releves.sb$cov <- as.numeric(releves.sb$cov)

#Kreuztabelle
releves.sb.wide <- cast(releves.sb, plot~abbr, value='cov')
releves.sb.wide[is.na(releves.sb.wide)] <- 0


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

rankindex(scale(header.sb$RFA), releves.sb.wide, c("euc","man","bray","jac","kul"))
dis <- vegdist(releves.sb.wide, method = 'jaccard')
adonis(dis ~ RFA * msm + geology, header.sb)
