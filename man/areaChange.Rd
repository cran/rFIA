\name{areaChange}
\alias{areaChange}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Estimate land area change from the FIADB}
\description{
Produces estimates of annual net and component change in land area (acreage) from FIA data. Estimates can be produced for regions defined within the FIA Database (e.g. counties), at the plot level, or within user-defined areal units. Options to group estimates by species, size class, and other variables defined in the FIADB. If multiple reporting years (EVALIDs) are included in the data, estimates will be output as a time series. If multiple states are represented by the data, estimates will be output for the full region (all area combined), unless specified otherwise (e.g. grpBy = STATECD). }


\usage{
areaChange(db, grpBy = NULL, polys = NULL, returnSpatial = FALSE,
           byLandType = FALSE, landType = "forest", method = "TI",
           lambda = 0.5, treeDomain = NULL, areaDomain = NULL,
           totals = FALSE, variance = FALSE, byPlot = FALSE,
           condList = FALSE, chngType = 'net', nCores = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{db}{\code{FIA.Database} or \code{Remote.FIA.Database} object produced from \code{\link{readFIA}} or \code{\link{getFIA}}. If a \code{Remote.FIA.Database}, data will be read in and processed state-by-state to conserve RAM (see details for an example).}

\item{grpBy}{variables from PLOT, COND, or TREE tables to group estimates by (NOT quoted). Multiple grouping variables should be combined with \code{c()}, and grouping will occur heirarchically. For example, to produce seperate estimates for each ownership group within ecoregion subsections, specify \code{c(ECOSUBCD, OWNGRPCD)}.}

  \item{polys}{\code{sp} or \code{sf} Polygon/MultiPolgyon object; Areal units to bin data for estimation. Seperate estimates will be produces for region encompassed by each areal unit. FIA plot locations will be reprojected to match projection of \code{polys} object.}

  \item{returnSpatial}{logical; if TRUE, merge population estimates with \code{polys} and return as \code{sf} multipolygon object. When \code{byPlot = TRUE}, return plot-level estimates as \code{sf} spatial points.}

  \item{byLandType}{logical; if TRUE, return estimates grouped by individual land type classes ("timberland", "non-timberland forest", "non-forest", and "water").}

  \item{landType}{character, one of: "forest", "non-forest", "census water", "non-census water", "water", or "all"; Type of land that estimates will be produced for. Timberland is a subset of forestland (default) which has high site potential and non-reserve status (see details).}

  \item{method}{character; design-based estimator to use. One of:  "TI" (temporally indifferent, default), "annual" (annual), "SMA"" (simple moving average), "LMA" (linear moving average), or "EMA" (exponential moving average). See \href{https://www.fs.fed.us/nrs/pubs/jrnl/2020/nrs_2020_stanke_001.pdf}{Stanke et al 2020} for a complete description of these estimators.}

  \item{lambda}{numeric (0,1); if \code{method = 'EMA'}, the decay parameter used to define weighting scheme for annual panels. Low values place higher weight on more recent panels, and vice versa. Specify a vector of values to compute estimates using mulitple wieghting schemes, and use \code{plotFIA} with \code{grp} set to \code{lambda} to produce moving average ribbon plots. See \href{https://www.fs.fed.us/nrs/pubs/jrnl/2020/nrs_2020_stanke_001.pdf}{Stanke et al 2020} for examples.}

  \item{treeDomain}{logical predicates defined in terms of the variables in PLOT, TREE, and/or COND tables. Used to define the type of trees for which estimates will be produced (e.g. DBH greater than 20 inches: \code{DIA > 20}, Dominant/Co-dominant crowns only: \code{CCLCD \%in\% 2:3)}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only trees where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{areaDomain}{logical predicates defined in terms of the variables in PLOT and/or COND tables. Used to define the area for which estimates will be produced (e.g. within 1 mile of improved road: \code{RDDISTCD \%in\% 1:6}, Hard maple/basswood forest type: \code{FORTYPCD == 805)}. Multiple conditions are combined with \code{&} (and) or \code{|} (or). Only plots within areas where the condition evaluates to TRUE are used in producing estimates. Should NOT be quoted.}

  \item{totals}{logical; if TRUE, return total population estimates (e.g. total area) along with ratio estimates (e.g. mean trees per acre).}

    \item{variance}{logical; if TRUE, return estimated variance (\code{VAR}) and sample size (\code{N}). If FALSE, return 'sampling error' (\code{SE}) as returned by EVALIDator. Note: sampling error cannot be used to construct confidence intervals.}

  \item{byPlot}{logical; if TRUE, returns estimates for individual plot locations instead of population estimates.}

  \item{condList}{logical; if TRUE, returns condition-level summaries intended for subsequent use with \code{\link{customPSE}}.}

  \item{chngType}{character, one of "net" or "component"; if "net", produce estimates of net change in land area, and if "component", produce estimates of component change in land area (i.e., showing all shifts across classified attributes).}

  \item{nCores}{numeric; number of cores to use for parallel implementation. Check available cores using \code{\link{detectCores}}. Default = 1, serial processing.}
}


\details{
\strong{Estimation Details}

Estimation of forest variables follows the procedures documented in Bechtold and Patterson (2005) and \href{https://www.fs.fed.us/nrs/pubs/jrnl/2020/nrs_2020_stanke_001.pdf}{Stanke et al 2020}.

Estimates are returned in terms of \emph{net} annual change in land area by default, however users may choose to estimate components of land area change by setting \code{chngType = 'component'}. For example, imagine we are interested in estimating the change in forestland area across a region. During our sampling period, 4000 acres of forestland was \emph{diverted} to non-forest and an additional 6000 acres of non-forest \emph{reverted} to forestland. rFIA considers these shifts in land classifications change \emph{components}, and hence these point estimates would be returned when \code{chngType = 'component'}. However, we are often intereseted in \emph{net} change in land area, rather than individual components. Here net change in forestland area is +2000 acres (6000-4000) and represents the net result of diversion and reversion processes in the region over our study period.

Users may specify alternatives to the 'Temporally Indifferent' estimator using the \code{method} argument. Alternative design-based estimators include the annual estimator ("ANNUAL"; annual panels, or estimates from plots measured in the same year), simple moving average ("SMA"; combines annual panels with equal weight), linear moving average ("LMA"; combine annual panels with weights that decay \emph{linearly} with time since measurement), and exponential moving average ("EMA"; combine annual panels with weights that decay \emph{exponentially} with time since measurement). The "best" estimator depends entirely on user-objectives, see \href{https://www.fs.fed.us/nrs/pubs/jrnl/2020/nrs_2020_stanke_001.pdf}{Stanke et al 2020} for a complete description of these estimators and tradeoffs between precision and temporal specificity.

When \code{byPlot = FALSE} (i.e., population estimates are returned), the "YEAR" column in the resulting dataframe indicates the final year of the inventory cycle that estimates are produced for. For example, an estimate of current forest area (e.g., 2018) may draw on data collected from 2008-2018, and "YEAR" will be listed as 2018 (consistent with EVALIDator). However, when \code{byPlot = TRUE} (i.e., plot-level estimates returned), the "YEAR" column denotes the year that each plot was measured (MEASYEAR), which may differ slightly from its associated inventory year (INVYR).

Stratified random sampling techniques are most often employed to compute estimates in recent inventories, although double sampling and simple random sampling may be employed for early inventories. Estimates are adjusted for non-response bias by assuming attributes of non-response plot locations to be equal to the mean of other plots included within thier respective stratum or population.



\strong{Working with "Big Data"}

If FIA data are too large to hold in memory (e.g., R throws the "cannot allocate vector of size ..." errors), use larger-than-RAM options. See documentation of \code{link{readFIA}} for examples of how to set up a \code{Remote.FIA.Database}. As a reference, we have used rFIA's larger-than-RAM methods to estimate forest variables using the entire FIA Database (~50GB) on a standard desktop computer with 16GB of RAM. Check out \href{https://rfia.netlify.app/}{our website} for more details and examples.

Easy, efficient parallelization is implemented with the \code{\link{parallel}} package. Users must only specify the \code{nCores} argument with a value greater than 1 in order to implement parallel processing on their machines. Parallel implementation is achieved using a snow type cluster on any Windows OS, and with multicore forking on any Unix OS (Linux, Mac). Implementing parallel processing may substantially decrease free memory during processing, particularly on Windows OS. Thus, users should be cautious when running in parallel, and consider implementing serial processing for this task if computational resources are limited (\code{nCores = 1}).



\strong{Definition of forestland}

Forest land must be at least 10-percent stocked by trees of any size, including land that formerly had such tree cover and that will be naturally or artificially regenerated. Forest land includes transition zones, such as areas between heavily forested and nonforested lands that are at least 10-percent stocked with trees and forest areas adjacent to urban and builtup lands. The minimum area for classification of forest land is 1 acre and 120 feet wide measured stem-to-stem from the outer-most edge. Unimproved roads and trails, streams, and clearings in forest areas are classified as forest if less than 120 feet wide. Timber land is a subset of forest land that is producing or is capable of producing crops of industrial wood and not withdrawn from timber utilization by statute or administrative regulation. (Note: Areas qualifying as timberland are capable of producing at least 20 cubic feet per acre per year of industrial wood in natural stands. Currently inaccessible and inoperable areas are NOT included).
}


\value{
Dataframe or SF object (if \code{returnSpatial = TRUE}). If \code{byPlot = TRUE}, values are returned for each plot (proportion of plot in domain of interest; \code{PLOT_STATUS_CD = 1} when forest exists at the plot location). All variables with names ending in \code{SE}, represent the estimate of sampling error (\%) of the variable. When \code{variance = TRUE}, variables ending in \code{VAR} denote the variance of the variable and \code{N} is the total sample size (i.e., including non-zero plots).

\itemize{
      \item{\strong{YEAR}: reporting year associated with estimates}
      \item{\strong{PERC_CHNG}: estimate of annual percent change in land area within domain of interest (\% of previous)}
      \item{\strong{AREA_CHNG}: estimate of annual change in land area within domain of interest (acres)}
      \item{\strong{PREV_AREA}: estimate of total land area within domain of interest at first measurement (acres)}
      \item{\strong{nPlots}: number of non-zero plots used to compute area change estimates}
      }

Importantly, when \code{chngType = 'component'}, individual change components will be returned. If no grouping variables are specified in \code{grpBy}, results will be grouped by variables named \strong{STATUS1} and \strong{STATUS2}, indicating the land classification at first and second measurements, respectively. Otherwise, if \code{grpBy} is specified, change components will be estimated for all shifts in land area across classified attributes represented by the variables (first and second measurements again denoted by the suffix 1 and 2).
}

\references{
rFIA website: \url{https://rfia.netlify.app/}

FIA Database User Guide: \url{https://www.fia.fs.fed.us/library/database-documentation/}

Bechtold, W.A.; Patterson, P.L., eds. 2005. The Enhanced Forest Inventory and Analysis Program - National Sampling Design and Estimation Procedures. Gen. Tech. Rep. SRS - 80. Asheville, NC: U.S. Department of Agriculture, Forest Service, Southern Research Station. 85 p. \url{https://www.srs.fs.usda.gov/pubs/gtr/gtr_srs080/gtr_srs080.pdf}

Stanke, H., Finley, A. O., Weed, A. S., Walters, B. F., & Domke, G. M. (2020). rFIA: An R package for estimation of forest attributes with the US Forest Inventory and Analysis database. Environmental Modelling & Software, 127, 104664.
}


\author{
Hunter Stanke and Andrew Finley
}


\note{All sampling error estimates (SE) are returned as the "percent coefficient of variation" (standard deviation / mean * 100) for consistency with EVALIDator. IMPORTANT: sampling error cannot be used to construct confidence intervals. Please use \code{variance = TRUE} for that (i.e., return variance and sample size instead of sampling error).}



\seealso{\code{\link{area}}}

\examples{
## Load data from the rFIA package
data(fiaRI)
data(countiesRI)

## Most recents subset
fiaRI_mr <- clipFIA(fiaRI)

## Most recent estimates of change in forested area in RI
areaChange(db = fiaRI_mr)

\donttest{
## Same as above grouped by land class
areaChange(db = fiaRI_mr, byLandType = TRUE)

## Estimates for change in forest area where stems greater than 20 in DBH
## occur for all available inventories (time-series)
areaChange(db = fiaRI,
           landType = 'forest',
           treeDomain = DIA > 20)

## Same as above, but implemented in parallel (much quicker)
parallel::detectCores(logical = FALSE) # 4 cores available, we will take 2
areaChange(db = fiaRI,
           landType = 'forest',
           treeDomain = DIA > 20,
           nCores =2)

## Return estimates at the plot-level
areaChange(db = fiaRI,
           byPlot = TRUE)
  }
}
