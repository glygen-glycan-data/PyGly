
To use the build scripts, you will need an account on the respective
wikis: production (glycomotif), dev (glycomotifdev), and test
(glycomotiftest). In addition, define the connection parameters of each
wiki in your ~/.glycomotif.ini file:

    [glycomotif]
    protocol = 
    host = 
    port = 
    tsport = 
    username = 
    password = 
    
    [glycomotifdev]
    protocol = 
    host = 
    port = 
    tsport = 
    username = 
    password = 
    
    [glycomotiftest]
    protocol = 
    host = 
    port = 
    tsport = 
    username = 
    password =

By default, the scripts will attempt to change the dev wiki,
To explicitly set the environment to work with, the environment variable
SMWENV can be set to PROD, DEV, or TEST, or the first command-line
option of any script provided as --smwenv ENV with ENV one of PROD, DEV,
TEST. The scripts must be run from the EdwardsLab cluster as the update
URL for the triplestore is not available to the public. Note that these
connection parameters are not for public (anonymous) access to the wiki
and triplestore, in each case, these are:

    http://edwardslab.bcmb.georgetown.edu/<site>

where site is one of glycomotif, glycomotifdev, glycomotiftest. 

