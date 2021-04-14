\name{BLUP}
\alias{BLUP}
\title{
Best Linear Unbias Predictor
}
\description{
Genetic values for a given trait computed by REML.
}
\usage{
BLUP(trait="PlantHeight",family="all",env="all",dereg=TRUE,use_bins=TRUE)
}
\arguments{
  \item{trait}{
Character. Trait of interest. The options are: GerminationCount, StandCount, DaystoSilk, DaysToTassel, GDDDaystoSilk, GDDDaystoTassel, GDDAnthesis.SilkingInterval, TasselLength, MainSpikeLength, TasselPrimaryBranches, TasselBranchLength, SpikeletsMainSpike, SpikeletsPrimaryBranch, SecondaryBranchNumber, TasselSterility, PlantHeight, EarHeight, NumberofLeaves, LeafLength, LeafWidth, UpperLeafAngle, MiddleLeafAngle, TilleringIndex, LeafSheathLength, Phylotaxy, NumberofTilleringPlants, SizeofTillers, NR, CS, GS, GLDH, AlaAT, AGPase, ShikDH, NIRStarch, NIRProtein, NIROil, NorthernLeafBlight, PercentageNorthernLeafBlight, IncubationPeriodNorthernLeafBlight, SouthernLeafBlight, CobDiameter, CobWeight, EarDiameter, EarLength, SeedSetLength, EarRowNumber, TotalKernelVolume, KernelFillPercentage, EarWeight, EarRankNumber, KernelWeight20, HerbicideSensitivity.
}
  \item{family}{
Numberic vector or "all". Which NAM families to use, ranging from 1 to 26.
}
  \item{env}{
Numberic vector or "all". Which environments to use. The environments are coded as follows: 1 = Aurora_2006, 2 = Aurora_2007, 3 = Aurora_2008, 4 = Clayton_2006, 5 = Clayton_2007, 6 = Columbia_2006, 7 = Columbia_2007, 8 = Homestead_2006, 9 = Homestead_2007, 10 = Ithaca_2007, 11 = Ithaca_2009, 12 = Ponce_2006, 13 = Ponce_2008, 14 = Urbana_2006, 15 = Urbana_2007.
}
  \item{dereg}{
Logical. Deregress BLUPs?
}
  \item{use_bins}{
Logical. If false, the genomic information is based on SNPs that are segregating across all NAM families. If true, the genomic information will be based on fixed-length bins that inform the percentage of B73. Bins were build with approximately 50 SNP each, based on the high-density genotypic information based on more than 200K markers.
}

}
\details{

The algorithm start from selecting the chosen families and environment that will be used for the best linear unbias predictor (BLUP). The BLUP values are calculated based on the following model: \eqn{Trait = Environment + Genotype}), where both environment and genotypes are random effects. If de-regression is used, the shrinkage of BLUP values is removed to remove downward biased of QTL effects, and to ensure all genetic values are in the same scale regardless of unbalanced number of observations.

We encourga the use of bins for genomic prediction studies, and the use of SNPs for association studies. Since markers were QCed to be segregating across all families, genomic prediction based on this SNP set may not be meaningful for cross-validations across-families.

}
\value{

This function returns a list with four objects. A numeric vector with the  BLUP solution of the phenotyes ("Phen"); the corresponding genotypes ("Gen"); a vector with the respective family ("Fam"); and a numeric vector with the number of SNPs per cromosome ("Chrom"). The output of this fuction has the exact input format for the NAM package (Xavier et al. 2015) to perform genome-wide association analysis. For the association analysis using the NAM package we recommend using the function NAM::gwas2 if use_bins=FALSE and NAM::gwas3 is use_bins=TRUE.

}
\references{
Xavier, A., Xu, S., Muir, W. M., & Rainey, K. M. (2015). NAM: association studies in multiple populations. Bioinformatics, 31(23), 3862-3864.
}
\author{
Alencar Xavier
}
\examples{
Test=BLUP(trait="PlantHeight",family=2:7,env=1:6)
}