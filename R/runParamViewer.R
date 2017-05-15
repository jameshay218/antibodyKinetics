#' @export
paramViewer <- function() {
  appDir <- system.file("shiny-examples", "paramViewer", package = "antibodyKinetics")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `mypackage`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal",launch.browser = TRUE)
  
}