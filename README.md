### Install from GitHub with
`devtools::install_github('ricebrian/MaizeNAM')`

### Get the complete SNP data
`require(MaizeNAM)`
`DownloadData()`

### Get genotypic data + phenotypic value for any given trait
`GetData = BLUP(trait="PlantHeight",family="all",env="all")`

### Check data
`head(GetData$Phen) # Deregressed phenotypic BLUPs`
`dim(GetData$Gen) # Panel of SNPs segregating in all familier`

### Analysis
`# The data above is already formatted to run genomic analysis on multiple packages (NAM,bWGR,mas,BGLR,rrBLUP,glmnet)`
