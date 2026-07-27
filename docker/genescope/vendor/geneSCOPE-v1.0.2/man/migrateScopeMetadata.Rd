\name{migrateScopeMetadata}
\alias{migrateScopeMetadata}
\title{Migrate legacy object-level metadata out of gene metadata}
\usage{
migrateScopeMetadata(
  scope_obj,
  platform = NULL,
  dataset = NULL,
  drop_redundant_columns = TRUE,
  drop_legacy_pseudo_rows = TRUE
)
}
\arguments{
\item{scope_obj}{A \code{scope_object} to migrate.}

\item{platform}{Optional explicit platform override, such as
\code{"Xenium"}. This is recommended when a legacy object contains the
non-informative label \code{"geneSCOPE"}; the true platform is never guessed.}

\item{dataset}{Optional explicit dataset label override.}

\item{drop_redundant_columns}{Remove constant \code{platform}/\code{dataset}
columns from gene metadata. Non-constant columns are retained.}

\item{drop_legacy_pseudo_rows}{Remove reserved
\code{__scope_platform__} and \code{__scope_metadata__} rows.}
}
\value{
The migrated \code{scope_object}. Lee matrices, clustering assignments, graphs,
grid layers, coordinates, and density results are not changed.
}
\description{
Upgrades a \code{scope_object} created by geneSCOPE v1.0.0/v1.0.1 so platform
and dataset labels are stored under \code{scope_obj@stats$object_metadata}
instead of being repeated in \code{scope_obj@meta.data}. The legacy scalar
locations \code{scope_obj@stats$platform} and
\code{scope_obj@stats$dataset} are retained for backward compatibility.
}
\examples{
\dontrun{
scope_obj <- readRDS("legacy_scope.rds")
scope_obj <- migrateScopeMetadata(scope_obj)
}
}
