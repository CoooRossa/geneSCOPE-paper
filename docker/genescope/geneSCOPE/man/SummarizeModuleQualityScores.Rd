\name{SummarizeModuleQualityScores}
\alias{SummarizeModuleQualityScores}
\title{Summarize stored module-quality scores}
\description{
Produces module- and member-level summary tables from stored module-quality
results. The member table reports the top members per module.
}
\usage{
SummarizeModuleQualityScores(x, grid_name = NULL, stats_layer = NULL)
}
\arguments{
\item{x}{A \code{scope_object} containing stored module-quality scores.}
\item{grid_name}{Optional grid filter.}
\item{stats_layer}{Optional stats layer filter.}
}
\value{
A list with \code{module_table}, \code{member_table}, and \code{summary_text}
(plus \code{per_graph} and \code{per_module} for convenience).
}
\examples{
\dontrun{
summary <- SummarizeModuleQualityScores(scope_obj)
}
}
