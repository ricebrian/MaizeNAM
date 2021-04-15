\name{More data}
\alias{DownloadData}
\alias{LoadData}
\title{
Full GBS markers, bins, and more
}
\description{
Download data not included in the package due to size.
}
\usage{
DownloadData( dataset = "all", redownload = FALSE )
LoadData( dataset = "panel" )
}

\arguments{
  \item{dataset}{ String: Which dataset to load? Options are "panel", "gbs" and "bins". For DownloadData, "all" is also an option.
 }
 \item{redownload}{ Logical: Should the data be downloaded even if it is already in the package?
 }
}
\details{
DownloadData downloads three datasets: (1) the full set of SNPs obtained using genotype-by-sequencing, with missing values pre-imputed; (2) Maize diversity panel; and (3) genotypic information from synthetic bins built with 50 SNPs each from the GBS data.  
}
\author{
Alencar Xavier
}
\examples{
\dontrun{
DownloadData( dataset = "bins" )
LoadData( dataset = "bins" )
}


}
