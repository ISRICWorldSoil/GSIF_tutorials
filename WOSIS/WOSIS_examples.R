## Data access to WOSIS / examples of how to query
## http://www.isric.org/projects/data-wosis-project
## Prepared by E. Ribeiro & T. Hengl (tom.hengl@isric.org)

## WOSIS Postgresql Connection and exercises
## NOTE: WOSIS only accepts WUR IP for connections

## 1. We require DBI (CommonDatabaseInterface) and PostgreSQL package
## DBI is a generic package that connect to multiple databases, DBI is a standard procedure on database connection
## Normall steps are: Connection, Cursor, Operation, Commit, Close Connection
install.packages(c("DBI","RPostgreSQL"))

## 2. PostGreSQL is a server database, meaning it works in a server and clients structured
## one server may have multiple databases and accept multiple connection in parallel
## First load the driver to be used, then we defined a connection string. The connection string indicates
## database name, login (user and password), location of server (host) and where is the database listening (port)
library(RPostgreSQL)
drv <- dbDriver("PostgreSQL")
con <- dbConnect(drv, dbname="spring", user="student", host="81.169.159.7", port="5432", password="student")

## This error indicates that we couldnt connect maybe there is something wrong with connection string
#Error in postgresqlNewConnection(drv, ...) : 
#RS-DBI driver: (could not connect student@81.169.159.7 on dbname "spring"

## 3. Question the following commant reports what sort of information ??
dbGetInfo(con)
dbListTables(con)

## 4. First query lets count the number of profiles on database, in WOSIS the profiles 
## are described on the "Profile"."Profile" table

rs <- dbSendQuery(con, 'SELECT COUNT(*) FROM "Profile"."Profile"') ## rs is  a result statment / cursor, we need to get the data from it
## 5. see the output, it looks like a table, the n=-1 indicates that we will get all records
## Question how many records we have in database
fetch(rs,n=-1)

## 6. What happen if we fetch again the records??? Why 0 columns and 0 rows ??? The cursor reached the end of query no more data
fetch(rs,n=-1)

## 7. Lets get 10 profiles
rs <- dbSendQuery(con, 'SELECT * FROM "Profile"."Profile"')
## Note:
#Error in postgresqlExecStatement(conn, statement, ...) : 
#  RS-DBI driver: (could not run statement: SSL SYSCALL error: EOF detected
#  )
## This error means that SQL query had problems in this case the connection may have been terminated. Try to run the 
#con <- dbConnect......
result <- fetch(rs, n=10)

## Question: why do we get this warning?
# Warning messages:
#  1: In postgresqlExecStatement(conn, statement, ...) :
#  RS-DBI driver warning: (unrecognized PostgreSQL field type geometry (id:35390690) in column 3)
#2: In postgresqlExecStatement(conn, statement, ...) :
#  RS-DBI driver warning: (unrecognized PostgreSQL field type uuid (id:2950) in column 18)

#Error in postgresqlExecStatement(conn, statement, ...) : 
#  RS-DBI driver: (connection with pending rows, close resultSet before continuing)
dbClearResult(rs)

## 8. Better to limit the number of profiles in the SQL it self this is more efficient
rs <- dbSendQuery(con, 'SELECT * FROM "Profile"."Profile" LIMIT 10')
result <- fetch(rs,n=-1)

## 9. The Profile.Profile table has a field for countryID,
## we can use the country ID to get all profiles in Iberia (Portugal+Spain)

## a. First search For countryID in Location.Countries
## We can run a describe table to see the fields inside:
rs <- dbSendQuery(con, 'SELECT * FROM "Location"."Country" WHERE "Name"=\'Spain\' OR "Name"=\'Portugal\'')
fetch(rs,n=-1)

## Why do we have to use '\'?
## What is the countryId for Portugal and Spain?

rs <- dbSendQuery(con,'SELECT * FROM "Profile"."Profile" WHERE "CountryId"=172 OR "CountryId"=196')
data <- fetch(rs,n=-1)
View(data)
## Around 156 observations BUT We would like to have a X Y column hoew to do it?
rs <- dbSendQuery(con,'SELECT ST_X("Location") AS x, ST_Y("Location") AS y FROM "Profile"."Profile" WHERE "CountryId"=172 OR "CountryId"=196') 
points <- fetch(rs,n=-1)
View(points)

## Get all data from WoSIS currently available:
wosis.pnts <- dbSendQuery(con, "SELECT * FROM web_services.latest_all;")
wosis.pnts.tbl <- fetch(wosis.pnts, n=-1)
summary(wosis.pnts.tbl$bd_ws)
wosis.pnts.tbl = plyr::rename(wosis.pnts.tbl, c("latest_ph_cacl2"="ph_cacl2", "latest_ph_h2o"="ph_h2o", "latest_ph_kcl"="ph_kcl", "latest_ph_naf"="ph_naf"))
str(wosis.pnts.tbl)

wosis.xy <- dbSendQuery(con, "SELECT profile_id, country_id, datasets, latitude, longitude, geom_accuracy FROM web_services.latest_profiles;")
wosis.xy.tbl <- fetch(wosis.xy, n=-1)
str(wosis.xy.tbl)

wosis_profiles <- list(wosis.xy.tbl, wosis.pnts.tbl)
names(wosis_profiles) = c("sites", "horizons")
save(wosis_profiles, file="wosis_profiles.rda", compress="xz")

## b. We can query the DB for some actual soil data like pH 
library(plotKML)
library(sp)
library(rgdal)
sql='SELECT "ProfileId",avg(value::float),ST_X("Location") as x, ST_Y("Location") as y FROM "Profile"."ProfileWithLayerAttributeGeo" WHERE "Name" LIKE \'%pH (H2O)%\' AND value!=\'\' AND value::float<14 GROUP BY "ProfileId","Location" ORDER BY "ProfileId"'
rs <- dbSendQuery(con,sql)
dataPH <- fetch(rs,n=-1)
names(dataPH) <- c("ProfileId","pH","x","y")
## We have NA values lets remove them
dataPH <- dataPH[complete.cases(dataPH),]
coordinates(dataPH) <- ~x+y
proj4string(dataPH)  <- CRS("+init=epsg:4326")
plotKML(dataPH["pH"], colour_scale=R_pal[["pH_pal"]])
## This can take time as you are plotting thousands of points!

## end of script;