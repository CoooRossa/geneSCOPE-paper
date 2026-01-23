\name{ScoreModuleQuality}
\alias{ScoreModuleQuality}
\title{Score module quality using Hotspot-style statistics}
\description{
Computes member- and module-level cohesion/separation metrics from a graph
and membership labels. Optionally runs a label-permutation null to derive
z/p/FDR statistics.
}
\usage{
ScoreModuleQuality(
  x,
  membership = NULL,
  grid_name = "grid30",
  stats_layer = "LeeStats_Xz",
  graph_slot = "g_consensus",
  edge_weight = "weight",
  n_perm = 0L,
  seed = 1,
  use_cpp = TRUE,
  compute_conductance = TRUE,
  writeback = FALSE,
  verbose = getOption("geneSCOPE.verbose", TRUE),
  ...
)
}
\arguments{
\item{x}{A \code{scope_object}, \code{igraph}, or list with \code{graph} and optional \code{membership}.}
\item{membership}{Optional membership vector named by gene.}
\item{grid_name}{Grid layer name for \code{scope_object} input.}
\item{stats_layer}{Stats layer name for \code{scope_object} input.}
\item{graph_slot}{Graph slot name for \code{scope_object} input.}
\item{edge_weight}{Edge attribute used as weight.}
\item{n_perm}{Number of label permutations for null.}
\item{seed}{Random seed for permutations.}
\item{use_cpp}{Use C++ permutation if available.}
\item{compute_conductance}{Whether to compute conductance.}
\item{writeback}{Write results to \code{scope_obj@stats[[grid_name]][[stats_layer]]$module_quality}.}
\item{verbose}{Emit progress messages when TRUE.}
\item{...}{Reserved for future extensions.}
}
\value{
A list with \code{member_scores}, \code{module_scores}, and \code{meta}.
The member table includes \code{score_member} (within-minus-between mean edge
weight), \code{z_member}, \code{p_member}, and \code{padj_member} plus
edge-count and weight summaries. The module table includes \code{score_module}
and \code{effect_size} (within-minus-between), \code{p_module}, \code{padj_module},
and null mean/sd when permutations are enabled.
}
\details{
Writeback occurs only when \code{writeback = TRUE} and the input is a
\code{scope_object}. Existing fields are not modified.
}
\section{Interpretation}{
Member scores quantify how strongly a gene connects within its module versus
outside it. Module scores summarize the same contrast at the module level,
with permutation-based z/p values providing an objective, quantifiable
quality assessment.
}
\examples{
\dontrun{
res <- ScoreModuleQuality(scope_obj, grid_name = "grid30", n_perm = 50)
}
}
