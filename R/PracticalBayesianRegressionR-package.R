#' PracticalBayesianRegressionR: 'brms'-style regression with a Julia back end
#'
#' Write models as in 'brms' -- [brm()], [bf()], [set_prior()], families
#' such as [poisson()] or [negbinomial()] -- and fit them with
#' PracticalBayesianRegression.jl.
#'
#' Two back ends:
#' * **Julia** (`backend = "julia"`): a Julia session through 'JuliaCall',
#'   started on first use ([pbr_setup()]). Every model the Julia package
#'   supports.
#' * **Compiled models** (`backend = "binary"`): [compile_model()] builds a
#'   standalone executable for a model once (minutes) and keeps it in the
#'   model store ([model_store()]), shared with the Julia package. From then
#'   on [brm()] runs it on any data with the right columns without starting
#'   Julia; [compiled_models()] lists the store.
#'
#' The default, `backend = "auto"`, uses a compiled model when the store has
#' one for the model being fitted, and Julia otherwise.
#'
#' Options: `pbr.backend` (default back end), `pbr.juliac` (JuliaC command
#' line for compiling), `pbr.julia_dev` (path of a local
#' PracticalBayesianRegression.jl checkout). Environment variables:
#' `PBR_MODEL_STORE`, `PBR_JULIA_PROJECT`, `PBR_JULIAC`.
#' `options(pbr.julia = FALSE)` forbids starting Julia (only compiled models run).
#' @keywords internal
"_PACKAGE"

# stats families re-exported so `family = poisson()` works without stats:: in
# scripts that mask them
#' @importFrom stats gaussian binomial poisson Gamma inverse.gaussian
#' @export
stats::gaussian
#' @export
stats::binomial
#' @export
stats::poisson
#' @export
stats::Gamma
#' @export
stats::inverse.gaussian
