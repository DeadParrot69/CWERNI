% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/multiCWRsimSPEC.R
\name{multiCWRsimSPEC}
\alias{multiCWRsimSPEC}
\title{Multiple simulations of the Community Wide Rescue model (CWR) saving specific timepoints and also saving the output matrix for every simulation.}
\usage{
multiCWRsimSPEC(tmax, b1, d1, b2, d2, m12, m21, k1, abun_original,
  interval, nsim, wantedtimes)
}
\arguments{
\item{tmax}{Arbitrary units of time the simulation should be run}

\item{b1}{Birthrate of the mutants}

\item{d1}{Deathrate of the mutants}

\item{b2}{Birthrate of the residents}

\item{d2}{Death rate of the residents}

\item{m12}{Mutation rate of mutants to residents}

\item{m21}{Mutation rate of residents to mutants}

\item{k1}{Carrying capacity}

\item{abun_original}{A vector of initial abundances for each species such as those created by generate_spat_abund}

\item{interval}{After how many evenst should the population state be saved to the output matrix}

\item{nsim}{The number of simulations to be performed}

\item{wantedtimes}{Specific timepoints for which the RAC is to be plotted}
}
\value{
this function plots the RAC of the simulations for the times specified in Wantedtimes as .tiff files in the current working directory also writes a .csv file outputting the poplation matrix for every simulation performed
#'@author Timo van Eldijk

this function plots all the results of the simulations in a .tiff file in the current working directory
}
\description{
Multiple simulations of the Community Wide Rescue model (CWR) saving specific timepoints and also saving the output matrix for every simulation.
}
\examples{
multiCWRsimSPEC(10,0.6 , 0.1, 0.05, 0.1, 0, 0.0005, 16000,
generate_spat_abund(theta = 200,Ivec = rep(40,1),Jvec = c(16000)),
200, 2, c(15,30,50,75,100))
}
