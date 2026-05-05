Eval <- function(object, pdf, length_grid = 1001) {
  if (!is.function(pdf)) {
    stop("pdf must be a density function.", call. = FALSE)
  }
  grid <- object$grid
  if (is.null(grid)) {
    stop("object must be a fitted DensOLog or kernel-method object with a grid.", call. = FALSE)
  }
  eval_grid <- .eval_grid_from_breaks(grid, length_grid = length_grid)
  dx <- eval_grid[2L] - eval_grid[1L]
  fhat <- if (inherits(object, "DensOLog")) {
    ddensolog(object, eval_grid)
  } else {
    .density_kernel(object, eval_grid)
  }
  sqrt(sum((pdf(eval_grid) - fhat)^2) * dx)
}
