% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CTAssignment.R
\name{ctblock}
\alias{ctblock}
\title{Clinical Trial Blocking Randomization Function}
\usage{
ctblock(data = data.frame(), treatments = c("A", "B"),
  proportions = c(1, 1), treatmentnames = "Treatment")
}
\arguments{
\item{data}{A data frame containing patient information}

\item{treatments}{Name of all possible treatments to assign. Defaults to "A" and "B"}

\item{proportions}{A vector containing the relative probability with which eah treatment should be
chosen, when randomly assigned, or how many times each treatment should show in each block, when
blocking. Defaults to c(1,1)}

\item{treatmentnames}{Name of the data frame column where the assigned treatments should be written
Defaults to "Treatment"}
}
\description{
Assigns patients to treatments in clinical trials.
}
