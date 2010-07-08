".onAttach" <- function(lib, pkg) {
  mylib <- dirname(system.file(package = pkg))
  title <- packageDescription(pkg, lib = mylib)$Title
  ver <- packageDescription(pkg, lib = mylib)$Version
  author <- packageDescription(pkg, lib = mylib)$Author
  cat(paste("\n", pkg, ": ", title, "\nVersion: ", ver, "\nAuthor: ", author, "\n\n", sep=""))
}

