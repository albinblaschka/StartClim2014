

setwd("Y:/Institut3/AGRAM/StartClim2014/Albin/analyse")

library(RODBC) 


mdbConnect <- odbcConnectAccess2007('Vegetation_StartClim2014.accdb') 
odbcGetInfo(mdbConnect)

sqlTables(mdbConnect)
dta <- sqlFetch(mdbConnect, 'tbl_sobotik') 

odbcClose(mdbConnect)
