.onAttach <-function (lib, pkg) {
  packageStartupMessage(c("Package saemix, version ", as.character(packageVersion("saemix")), ", March 2024\n"),
    "  please direct bugs, questions and feedback to emmanuelle.comets@inserm.fr\n")
}
