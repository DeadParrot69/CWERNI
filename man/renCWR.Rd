% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/renCWR.R
\name{renCWR}
\alias{renCWR}
\title{Perform CWR simulation, Calculate Reny entropy (alfa=1, alfa=2), re-simulate neutral communities using SADISA parameters estimated on CWR simulation,
simulation neutral communities, calculate their Reny entropy (alfa=1, alfa=2) to generate a distribution.}
\usage{
renCWR(tmax, b1, d1, k1, interval, b2, d2, m12, m21, orgsimmaxcount,
  totalresamp)
}
\arguments{
\item{tmax}{Arbitrary units of time the simulation should be run}

\item{b1}{Birthrate of the mutants}

\item{d1}{Deathrate of the mutants}

\item{k1}{Carrying capacity}

\item{interval}{After how many evenst should the population state be saved to the output matrix}

\item{b2}{Birthrate of the residents}

\item{d2}{Death rate of the residents}

\item{m12}{Mutation rate of mutants to residents}

\item{m21}{Mutation rate of residents to mutants}

\item{orgsimmaxcount}{how many CWR simulations should be performed}

\item{totalresamp}{how many neutral communities should be simulated for each single CWR simulation}
}
\value{
This function saves two .pdf plots for each CWR simulation performed, showing the Reny entropy of each CWR simulation and the distribution of Reny entorpy generated from the neutral re-simulations
}
\description{
Perform CWR simulation, Calculate Reny entropy (alfa=1, alfa=2), re-simulate neutral communities using SADISA parameters estimated on CWR simulation,
simulation neutral communities, calculate their Reny entropy (alfa=1, alfa=2) to generate a distribution.
}
\examples{
renCWR(10, 0.6, 0.1, 16000,200, 0.05,0.1, 0, 0.0005, 1,1)

}
\author{
Timo van Eldijk
}
