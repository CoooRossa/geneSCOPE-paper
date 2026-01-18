\name{GetModuleQualityScores}
\alias{GetModuleQualityScores}
\title{Extract stored module-quality scores}
\description{
Reads module-quality results stored under \code{module_quality} in a
\code{scope_object} and returns member/module tables (including results stored
under \code{module_quality$hotspot_style}).
}
\usage{
GetModuleQualityScores(x, grid_name = NULL, stats_layer = NULL)
}
\arguments{
\item{x}{A \code{scope_object} or a list containing \code{module_quality}.}
\item{grid_name}{Optional grid filter.}
\item{stats_layer}{Optional stats layer filter.}
}
\value{
A list with \code{member_scores} and \code{module_scores}.
}
\examples{
\dontrun{
res <- GetModuleQualityScores(scope_obj, grid_name = "grid30")
}
}
