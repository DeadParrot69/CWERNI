% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/loglikCWR.R
\name{loglikCWR}
\alias{loglikCWR}
\title{Perform CWR simulation, estimate paramters using SADISA, re-simulate neutral communities to generate distribution of likelihoods}
\usage{
loglikCWR(tmax, b1, d1, k1, interval, b2, d2, m12, m21, orgsimmaxcount,
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
This function saves a .pdf plot for each CWR simulation performed, showing the distribution of likelihoods
}
\description{
Perform CWR simulation, estimate paramters using SADISA, re-simulate neutral communities to generate distribution of likelihoods
}
\examples{
loglikCWR(10, 0.6, 0.1, 16000,200, 0.05,0.1, 0, 0.0005, 1,1)




}
\author{
Timo van Eldijk
}
