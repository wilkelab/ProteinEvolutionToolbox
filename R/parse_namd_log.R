parse.log <- function(file.name, column) {
  log.file <- readLines(file.name)

  return.values <- c()

  for(line in log.file) {
    tmp.line <- unlist(strsplit(line, " "))
    tmp.line <- tmp.line[tmp.line != ""]

    if(length(tmp.line) >= 1 && tmp.line[1] == "ENERGY:")
      return.values <- append(return.values, tmp.line[column])
  }

  return(return.values)
}