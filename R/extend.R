extend <- function (baseargs, values) {
	nval <- sapply(values, length)
	nva <- prod(nval)
	args <-  rep(list(baseargs), nva)
	valuedf <- do.call(expand.grid, values)
	for (i in 1:nrow(valuedf)) args[[i]][['details']][names(values)] <- valuedf[i,]
	args
}

