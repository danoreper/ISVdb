# ISVdb
Code to generate the inbred strain variant database (v1.1), as well as to build the GUI. 
Builds imputed genotypes and diplotypes for CC strains, and also stores founder info.

# Requirements
* VCFTools>0.1.13
* R>=3.2.3
* Python>=2.7.12
* MySQL or MariaDB

# To run database construction:
1. Working directory must be ISVDB_LOCATION/src
2. Edit ../config/defaultvdb.yaml to specify an appropriate value for the user (currently root) and password (currently XXX), 
such that the user chosen has sufficient privileges to construct a database.
3. Enter the following at command line from the working directory:
ISVDB_LOCATION/src$ R CMD BATCH '--args ../config/defaultvdb.yaml' ./genomerep/variantdb/dbdriver.R

# Directory structure
TODO...


