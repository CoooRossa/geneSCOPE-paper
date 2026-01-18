#' Compute spatial weights for a grid layer (public wrapper).
#' @description
#' Builds spatial neighbourhood structures and stores adjacency/weights into the
#' specified grid layer of a `scope_object`.
#' @param scope_obj A `scope_object` with grid layers.
#' @param grid_name Grid layer name (defaults to active layer).
#' @param style Style for `spdep::nb2listw()` (default `"B"`).
#'   Kernel weights (main-2 parity) are enabled by setting `style = "kernel_gaussian"`
#'   or `style = "kernel_flat"`. In that case, use `kernel_radius`/`kernel_sigma`
#'   to control the kernel shape. For backward compatibility, `max_order` and
#'   `min_neighbors` are only used as fallbacks when the kernel-specific
#'   parameters are not supplied.
#' @param store_mat Store the adjacency matrix in the grid layer when TRUE.
#' @param zero_policy Whether zero-neighbour cells are allowed.
#' @param store_listw Store the `listw` object when TRUE.
#' @param verbose Emit progress messages when TRUE.
#' @param topology Topology choice (`auto`, `queen`, `rook`, `hex`, `fuzzy_queen`, `fuzzy_hex`).
#' @param nb_order Optional higher-order neighbour expansion.
#' @param min_neighbors Minimum neighbour count when expanding (non-kernel styles).
#' @param max_order Maximum neighbour order (non-kernel styles).
#' @param kernel_radius Kernel radius for `style="kernel_*"` (overrides `max_order`).
#' @param kernel_sigma Kernel sigma for `style="kernel_*"` (overrides `min_neighbors`).
#' @param ncores Number of threads to use.
#' @return The modified `scope_object` (invisibly).
#' @examples
#' \dontrun{
#' scope_obj <- computeWeights(scope_obj, grid_name = "grid30", topology = "queen")
#' }
#' @seealso `computeL()`, `fuzzy_queen_jaccard()`
#' @export
computeWeights <- function(
    scope_obj,
    grid_name = NULL,
    style = "B",
    store_mat = TRUE,
    zero_policy = TRUE,
    store_listw = TRUE,
    verbose = TRUE,
    topology = c("auto", "queen", "rook", "hex", "fuzzy_queen", "fuzzy_hex"),
    nb_order = 1L,
    min_neighbors = NULL,
    max_order = NULL,
    kernel_radius = NULL,
    kernel_sigma = NULL,
    ncores = detectCores(logical = TRUE)) {
    .compute_weights(
        scope_obj = scope_obj,
        grid_name = grid_name,
        style = style,
        store_mat = store_mat,
        zero_policy = zero_policy,
        store_listw = store_listw,
        verbose = verbose,
        topology = topology,
        nb_order = nb_order,
        min_neighbors = min_neighbors,
        max_order = max_order,
        kernel_radius = kernel_radius,
        kernel_sigma = kernel_sigma,
        ncores = ncores
    )
}
