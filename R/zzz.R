.onAttach <-function (lib, pkg) {
  packageStartupMessage(c("Package saemix, version ", as.character(packageVersion("saemix")), "\n"),
    "  please direct bugs, questions and feedback to emmanuelle.comets@inserm.fr\n")
}
