#' Attach histology assets to a scope object.
#' @description
#' Adds Visium-style histology images and scalefactors to an existing grid layer
#' within a `scope_object`.
#' @param scope_obj A `scope_object` to annotate.
#' @param grid_name Grid layer name to update.
#' @param png_path Path to the histology PNG.
#' @param json_path Optional Visium scalefactors JSON; used when `scalefactors` is not provided.
#' @param level Image resolution level (`lowres` or `hires`).
#' @param crop_bbox_px Optional pixel-space bounding box to crop the PNG before attaching.
#' @param roi_bbox Optional ROI bounding box metadata stored alongside the image.
#' @param coord_type Coordinate system for scaling (`visium` or `manual`).
#' @param scalefactors Optional list of scalefactors to use instead of the JSON file.
#' @param y_origin Coordinate origin policy (`auto`, `top-left`, `bottom-left`).
#' @return Updated `scope_object` with histology metadata.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir")
#' scope_obj <- addScopeHistology(scope_obj, grid_name = "grid55", png_path = "/path/to/spatial/tissue_hires_image.png", json_path = "/path/to/spatial/scalefactors_json.json")
#' }
#' @seealso `plotGridBoundary()`, `createSCOPE_visium()`
#' @export
addScopeHistology <- function(
    scope_obj,
    grid_name,
    png_path,
    json_path = NULL,
    level = c("lowres", "hires"),
    crop_bbox_px = NULL,
    roi_bbox = NULL,
    coord_type = c("visium", "manual"),
    scalefactors = NULL,
    y_origin = c("auto", "top-left", "bottom-left")) {
    if (missing(level)) {
        level <- NULL
    }
    .add_scope_histology(
        scope_obj = scope_obj,
        grid_name = grid_name,
        png_path = png_path,
        json_path = json_path,
        level = level,
        crop_bbox_px = crop_bbox_px,
        roi_bbox = roi_bbox,
        coord_type = coord_type,
        scalefactors = scalefactors,
        y_origin = y_origin
    )
}

#' Add single-cell spatial data.
#' @description
#' Attaches platform-specific single-cell layers (CosMx or Xenium) to an existing
#' `scope_object`.
#' @param scope_obj A `scope_object` to extend.
#' @param xenium_dir Optional Xenium run directory.
#' @param cosmx_root Optional CosMx project root.
#' @param platform Platform hint (`auto`, `Xenium`, or `CosMx`).
#' @param verbose Emit progress messages when TRUE.
#' @param ... Additional arguments passed to the platform-specific loader.
#'   For backward compatibility, `data_dir` can be supplied here as an alias for
#'   `xenium_dir` / `cosmx_root` (auto-detected when possible).
#' @return Updated `scope_object` with single-cell layers.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE(data_dir = "/path/to/output_dir", grid_length = 30)
#' scope_obj <- addSingleCells(scope_obj, data_dir = "/path/to/xenium_run")
#' }
#' @seealso `createSCOPE_xenium()`, `createSCOPE_cosmx()`
#' @export
addSingleCells <- function(
    scope_obj,
    xenium_dir = NULL,
    cosmx_root = NULL,
    platform = c("auto", "Xenium", "CosMx"),
    verbose = TRUE,
    ...) {
    .add_single_cells(
        scope_obj = scope_obj,
        xenium_dir = xenium_dir,
        cosmx_root = cosmx_root,
        platform = platform,
        verbose = verbose,
        ...
    )
}

#' Create a scope object from platform outputs.
#' @description
#' Auto-detects Xenium/CosMx/Visium inputs under `data_dir` (or uses an explicit
#' `prefer` override) and dispatches to the matching constructor.
#' @param data_dir Root directory containing Xenium/CosMx/Visium outputs.
#' @param prefer Override platform detection (`auto`, `xenium`, `cosmx`, `visium`).
#' @param verbose Emit progress messages when TRUE.
#' @param sctransform Run SCTransform when using the Visium Seurat builder.
#' @param ... Additional arguments (currently unused).
#' @return A constructed `scope_object`.
#' @examples
#' \dontrun{
#' scope_obj <- createSCOPE("/path/to/output_dir", grid_length = 30)
#' scope_obj <- createSCOPE("/path/to/visium_run")
#' }
#' @seealso `createSCOPE_visium()`, `createSCOPE_xenium()`, `createSCOPE_cosmx()`
#' @export
createSCOPE <- function(
    data_dir = NULL,
    prefer = c("auto", "xenium", "cosmx", "visium"),
    verbose = TRUE,
    sctransform = TRUE,
    ...) {
    .create_scope(
        data_dir = data_dir,
        prefer = prefer,
        verbose = verbose,
        sctransform = sctransform,
        ...
    )
}
