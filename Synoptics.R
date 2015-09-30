

setwd("Y:/Institut3/AGRAM/StartClim2014/Albin/analyse")

library(RODBC) 
library(synoptic)

mdbConnect <- odbcConnectAccess2007('Vegetation_StartClim2014.accdb') 

rels <- sqlQuery(mdbConnect, "SELECT plot, species as taxon, layer, cov FROM tbl_sobotik", errors = TRUE)
specs <- sqlQuery(mdbConnect, "select species as taxon from Gesamtartenliste", errors = TRUE)
odbcClose(mdbConnect)

rels <- rels[rels$cov !='r',]

specs$abbr <- make.cepnames(specs$taxon)
rels <- merge(rels,specs, by='taxon')
rels$taxon <- NULL

rels$cov <- gsub("+","0.1", rels$cov, fixed = TRUE)
rels$cov <- as.numeric(rels$cov)

rels <- rels[,c('plot','abbr','layer','cov')]
specs <- specs[,c('abbr','taxon')]
head(rels)
head(specs)

specObj <- new('Species',rels)
siteObj <- new('Sites',rels)
taxObj <-  new('Taxonomy', specs)

st0 <- new('SpeciesTaxonomy')
st0@species <- specObj
st0@taxonomy <- taxObj

Vegsoup(st0,siteObj, taxObj, 'percentage', verbose = TRUE)


