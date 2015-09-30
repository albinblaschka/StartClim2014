

setwd("Y:/Institut3/AGRAM/StartClim2014/Albin/analyse")

library(RODBC) 
library(taxize)
library(synoptic)

mdbConnect <- odbcConnectAccess2007('Vegetation_StartClim2014.accdb') 

#res <- gnr_resolve(names = qry$species, data_source_ids = '158')
#res

#correct <- tnrs(query = qry$species, source = "iPlant_TNRS")[ , -c(5:7)]
#correct

rels <- sqlQuery(mdbConnect, "SELECT plot, species, layer, cov FROM tbl_sobotik", errors = TRUE)
specs <- sqlQuery(mdbConnect, "select species from Gesamtartenliste", errors = TRUE)
odbcClose(mdbConnect)

specs$abbr <- make.cepnames(specs$species)
rels <- merge(rels,specs, by='species')

head(rels)
head(specs)

specObj <- new('Species',rels)
siteObj <- new('Sites',rels)
taxObj <-  new("Taxonomy", specs)




