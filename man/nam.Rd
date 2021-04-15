\name{Dataset}
\alias{geno}
\alias{bins}
\alias{pheno}
\alias{Z}
\alias{B}
\alias{NamData}
\alias{MetaData}
\title{
Dataset
}
\description{
Genotypes and phenotypes from the Maize Nested Association Mapping (NAM) dataset. See the section "details" for the description of data objects.

The Maize NAM population (panzea.org) is a nested association mapping panel designed to dissect the genetic architecture of complex traits. It comprises about 5000 recombinant inbred lines (RILs) derived from biparental crosses between a common parent, B73, and 25 selected stiff stalk, non-stiff stalk, popcorn, sweetcorn, tropical, and mixed lines, i.e., the founder lines. Each biparental population approximately contains 200 individuals. Additional genotype and phenotype data for 281 members from the Goodman-Buckler diversity panel (Flint-Garcia et al. 2004) can be download using LoadData("panel").   

Single nucleotide polymorphisms were collected using genotype by sequencing for 4984 RILs and 281 diverse inbred lines (Glaubitz et al. 2014). Chromosome and basepair coordinates were uplifted maize B73 reference genome AGPv4.

}
\usage{
data(pheno)
data(geno)
}
\details{

Datasets of the MaizeNAM project was downloaded on April 2nd 2020 from \url{panzea.org}. Phenotypic data was collected at four locations (Aurora, NY; Clayton, NC; Columbia, MO; and Urbana, IL) between 2006 - 2009 (see BLUP). Studies performed on the entire dataset with additional detail about the experimental design include Brown et al., 2011; Buckler et al., 2009; Cook et al., 2012; H. Y. Hung et al., 2012a; H. Y. Hung et al., 2012b; Kump et al., 2011; Peiffer et al., 2014; Poland et al., 2011; Tian et al., 2011. Please note that GBS data contained in this package was not avalible at the time these were studies were published. 

Genotypic matrix "B" contains bin information for individual lines and matrix "Z" contains 4085 SNPs segregating across RIL familes. Matrix "X" and "Y", downloaded using LoadData("panel"), contain 269365 SNPs. 

In each dataset, phenotypes are allocated into two objects, one with the lines ("data.line") and one with checks and parents ("data.checks"). Information on data objects include year, location, environment (combination of year and location), strain, family, set (set in each environment), spot (combination of set and environment), height (in centimeters), R8 (number of days to maturity), planting date (501 represents may 1), flowering (701 represents july 1), maturity (901 represents september 1), lodging (score from 1 to 5), yield (in Kg/ha), moisture, protein (percentage in the seed), oil (percentage in the seed), fiber (percentage in the seed), seed size (in grams of 100 seeds).


The dataset "NAMGeno" and "Diversitygeno" contains the genotypic information in the matrix form, with missing values imputed using the software LinkImputed (Money et al. 2015). Note: these datasets have not been filterd for minor allele frequency. Datasets "NAMmap" and "Diversitymap" contain the accompnaying SNP allele and positional information in data.frame format. 

Similar to the datasets previously described, phenotypes are allocated into two objects, lines ("data.line.in") and checks ("data.checks.in"). Information on these data objects include year, location, environment (combination of year and location), strain, family, set (set in each environment), spot (combination of set and environment), the spatial coordinates of the field plots (BLOCK, ROW and COLUMN), plant height (in centimeters), R1 (number of days to flowering), R8 (number of days to maturity), lodging (score from 1 to 5), yield (in bu/ac), leaf shape (ratio length:width), number of nodes in the main stem, number of pods in the main stem, number of pods per node, average canopy coverage, rate of canopy coverage, growing degree day to flowering (GDD_R1), growing degree day to maturity (GDD_R1), and length of reproductive period in terms of growing degree day (GDD_REP).

}
\author{
Brian Rice, Alencar Xavier
}
\references{

Brown, P. J., Upadyayula, N., Mahone, G. S., Tian, F., Bradbury, P. J., Myles, S., Holland, J. B., Flint-Garcia, S., McMullen, M. D., Buckler, E. S., & Rocheford, T. R. (2011). Distinct genetic architectures for male and female inflorescence traits of maize. PLoS Genetics, 7(11). https://doi.org/10.1371/journal.pgen.1002383

Buckler, E. S., Holland, J. B., Bradbury, P. J., Acharya, C. B., Brown, P. J., Browne, C., Ersoz, E., Flint-Garcia, S., Garcia, A., Glaubitz, J. C., Goodman, M. M., Harjes, C., Guill, K., Kroon, D. E., & Larsson, S. (2009). The Genetic Architecture of Maize Flowering Time. Science, 325(2009), 714–718. https://doi.org/10.1126/science.1174276

Cook, J. P., McMullen, M. D., Holland, J. B., Tian, F., Bradbury, P., Ross-Ibarra, J., Buckler, E. S., & Flint-Garcia, S. A. (2012). Genetic Architecture of Maize Kernel Composition in the Nested Association Mapping and Inbred Association Panels. Plant Physiology, 158(2), 824–834. https://doi.org/10.1104/pp.111.185033

Flint-Garcia, S.A., Thuillet, A.C., Yu, J., Pressoir, G., Romero, S.M., Mitchell, S. E., et al. (2005) Maize association population: a high-resolution platform for quantitative trait locus dissection. Plant J. 44: 1054–1064

Glaubitz JC, Casstevens TM, Lu F, Harriman J, Elshire RJ, Sun Q, et al. TASSEL-GBS: A high capacity genotyping by sequencing analysis pipeline. PLoS One. 2014;9(2)

Hung, H.-Y., Shannon, L. M., Tian, F., Bradbury, P. J., Chen, C., Flint-Garcia, S. A., McMullen, M. D., Ware, D., Buckler, E. S., Doebley, J. F., & Holland, J. B. (2012). ZmCCT and the genetic basis of day-length adaptation underlying the postdomestication spread of maize. Proceedings of the National Academy of Sciences, 109(28), E1913–E1921. https://doi.org/10.1073/pnas.1203189109

Hung, H. Y., Browne, C., Guill, K., Coles, N., Eller, M., Garcia, A., Lepak, N., Melia-Hancock, S., Oropeza-Rosas, M., Salvo, S., Upadyayula, N., Buckler, E. S., Flint-Garcia, S., McMullen, M. D., Rocheford, T. R., & Holland, J. B. (2012). The relationship between parental genetic or phenotypic divergence and progeny variation in the maize nested association mapping population. Heredity, 108(5), 490–499. https://doi.org/10.1038/hdy.2011.103

Kump, K. L., Bradbury, P. J., Wisser, R. J., Buckler, E. S., Belcher, A. R., Oropeza-Rosas, M. A., Zwonitzer, J. C., Kresovich, S., McMullen, M. D., Ware, D., Balint-Kurti, P. J., & Holland, J. B. (2011). Genome-wide association study of quantitative resistance to southern leaf blight in the maize nested association mapping population. Nature Genetics, 43(2), 163–168. https://doi.org/10.1038/ng.747

Peiffer, J. A., Romay, M. C., Gore, M. A., Flint-Garcia, S. A., Zhang, Z., Millard, M. J., Gardner, C. A. C., McMullen, M. D., Holland, J. B., Bradbury, P. J., & Buckler, E. S. (2014). The genetic architecture of maize height. Genetics, 196(4), 1337–1356. https://doi.org/10.1534/genetics.113.159152

Poland, J. A., Bradbury, P. J., Buckler, E. S., & Nelson, R. J. (2011). Genome-wide nested association mapping of quantitative resistance to northern leaf blight in maize. Proceedings of the National Academy of Sciences of the United States of America, 108(17), 6893–6898. https://doi.org/10.1073/pnas.1010894108

Tian, F., Bradbury, P. J., Brown, P. J., Hung, H., Sun, Q., Flint-Garcia, S., Rocheford, T. R., McMullen, M. D., Holland, J. B., & Buckler, E. S. (2011). Genome-wide association study of leaf architecture in the maize nested association mapping population. Nature Genetics, 43(2), 159–162. https://doi.org/10.1038/ng.746


}
\examples{
data(pheno)
data(geno)
}
