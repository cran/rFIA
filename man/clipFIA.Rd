\name{clipFIA}
\alias{clipFIA}
\title{Spatial and temporal queries for FIADB
}

\description{
Performs space-time queries on Forest Inventory and Analysis Database (FIADB). Subset database to include only data associated with particular inventory years (i.e., most recent), and/or only data within a user-defined region. }

\usage{
clipFIA(db, mostRecent = TRUE, mask = NULL, matchEval = FALSE,
        evalid = NULL, designCD = NULL)
}

\arguments{
  \item{db}{FIA.Database object produced from \code{\link{readFIA}}.}

  \item{mostRecent}{logical; if TRUE, returns only data for most recent inventory.}

  \item{mask}{sp or sf Polygon/MultiPolgyon object; defines the boundaries of spatial intersection with FIA tables.}

  \item{matchEval}{logical; if TRUE, returns subset of data for which there are matching reporting years across states. Only useful if db contains mulitple state subsets of the FIA database.}

  \item{evalid}{character; unique value which identifies an inventory year and inventory type for a state. If you would like to subset data for an inventory year other than the most recent, use \code{\link{findEVALID}} to look locate this value (see Examples below).}

  \item{designCD}{character vector; plot designs to include. Default includes standard national plot design with other similar sampling designs. See FIA Database User Guide Appendix 1 for descriptions of plot designs (see References).}
}

\details{
Not required to run other \pkg{rFIA} functions, but may help conserve free memory and reduce processing time if user is interested in producing estimates for a specific inventory year or within a region not explicitly described in the database (w/in user defined polygons).

Spatial intersections do not adhere strictly to absolute plot locations, all plots which fall within an estimation unit (often a county) which intersects with a user defined region will be returned. The plots which fall slightly outside of the region do NOT bias estimates (removed from computations), but as FIA often employs stratified random sampling estimators, all plots within intersecting estimation units must be present to proudce unbiased variance estimates.
}

\value{
List object containing spatially intersected FIADB tables.
}

\author{
Hunter Stanke and Andrew Finley
}

\references{
FIA Database User Guide: \url{https://www.fia.fs.fed.us/library/database-documentation/}
}

\seealso{
\code{\link{findEVALID}}
}

\examples{
## Load data from rFIA package
data(fiaRI)

## Most recent inventory
clipFIA(fiaRI, mostRecent = TRUE)

## Only plots w/in estimation units w/in a user defined polygon
clipFIA(fiaRI, mask = countiesRI[1,])
}