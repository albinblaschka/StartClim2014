

setwd("Y:/Institut3/AGRAM/StartClim2014/Albin/analyse")

library(RODBC) 
library(taxize)

mdbConnect <- odbcConnectAccess2007('Vegetation_StartClim2014.accdb') 
odbcGetInfo(mdbConnect)
# sqlTables(mdbConnect)

qry <- sqlQuery(mdbConnect, "SELECT DISTINCT species FROM tbl_sobotik", errors = TRUE)
head(qry)
odbcClose(mdbConnect)

res <- gnr_resolve(names = qry$species, data_source_ids = '158')
res

correct <- tnrs(query = qry$species, source = "iPlant_TNRS")[ , -c(5:7)]
correct