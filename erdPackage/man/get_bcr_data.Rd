% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/BBS_functions.R
\name{get_bcr_data}
\alias{get_bcr_data}
\title{Aggregate eBird data over cells within BCRs}
\usage{
get_bcr_data(
  data,
  bcrs,
  .year,
  tgrid_min,
  tgrid_max,
  time_window = "full",
  min_lat,
  max_lat,
  min_lon,
  small_grid = 11,
  .cores = 4
)
}
\arguments{
\item{data}{filtered erd}

\item{bcrs}{sf representation of BCRs}

\item{.year}{single year for which index is desired}

\item{tgrid_min}{minimum week of focal period}

\item{tgrid_max}{maximum week of focal period}

\item{time_window}{one of "weekly" or "full"; controlls whether the bootstrap is stratified by week}

\item{min_lat}{minimum latitude}

\item{max_lat}{maximum latitude}

\item{min_lon}{minimum longitude}

\item{small_grid}{resolution of dggrid for small cells to use in bootstrapping}

\item{.cores}{number of cores for parallel processing; used only if time_window = "weekly"}
}
\value{
abundance data
}
\description{
Aggregate eBird data over cells within BCRs
}
