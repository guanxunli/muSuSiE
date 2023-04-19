library("pcalg")
library("Hmisc")

set.seed(1)

#' Auxiliary function bringing targets in a standard format.
#'
#' At the same time, the function checks if the targets are valid; if not,
#' it throws an exception.
#'
#' @param 	p				number of vertices
#' @param 	targets			list of (unique) targets
#' @param 	target.index	vector of target indices, or NULL
#' @return  depends on arguments:
#'   if target.index == NULL: list of sorted targets
#'   if target.index != NULL: list with two entries, "targets" and "target.index"
.tidyTargets <- function(p, targets, target.index = NULL) {
  stopifnot((p <- as.integer(p)) > 0)

  # Check and convert targets
  if (!is.list(targets) || !all(sapply(targets, is.numeric))) {
    stop("Argument 'targets' must be a list of integer vectors.")
  }
  rawTargets <- lapply(targets, function(v) unique(sort(as.integer(v))))
  targets <- unique(rawTargets)
  if (length(targets) < length(rawTargets)) {
    stop("List of targets must be unique.")
  }
  allTargets <- unlist(targets)
  if (length(allTargets) > 0) {
    if (any(is.na(allTargets))) {
      stop("Argument 'targets' must not contain NAs.")
    }
    min.max <- range(allTargets)
    if (min.max[1] <= 0 || min.max[2] > p) {
      stop("Targets are out of range.")
    }
  }

  # Check validity of target index, if provided
  if (!is.null(target.index)) {
    if (!is.numeric(target.index)) {
      stop("Argument 'target.index' must be an integer vector.")
    }
    target.index <- as.integer(target.index)
    min.max <- range(target.index)
    if (min.max[1] <= 0 || min.max[2] > length(targets)) {
      stop("Target index is out of range.")
    }
    # target.index <- match(rawTargets, targets)[target.index]
  }

  # Return value
  if (is.null(target.index)) {
    targets
  } else {
    list(targets = targets, target.index = target.index)
  }
}

setRefClass("MultiGaussL0pen",
	contains = "GaussL0penObsScore",
	fields = list(
		.gauss.vec = "list"),

	methods = list(
		#' Constructor
		initialize = function(data = list(matrix(1,1,1)),
		nodes = colnames(data[[1]]),
		lambda = 0,
		intercept = FALSE,
		format = c("raw", "scatter"),
		use.cpp = FALSE,
		...) {
			#transform the data format
			.gauss.vec <<- lapply(data, function(x) 
				new("GaussL0penObsScore", 
					data =  x,
					lambda = lambda / length(data), 
					intercept = intercept,
					use.cpp = TRUE))
			data = do.call(rbind, data)

			#call super class
			callSuper(data = data,
				nodes = nodes,
				lambda = lambda,
				intercept = intercept,
				format = format,
				use.cpp = use.cpp,
				...)
		},

		#' Calculates the local score of a vertex and its parents
		local.score = function(vertex, parents, ...) {
			return(sum(sapply(.gauss.vec, function(x) x$local.score(vertex, parents))))
		},

		#' Calculates the local mle
		local.fit = function(vertex, parents, ...) {
			return(c(sapply(.gauss.vec, function(x) x$local.fit(vertex, parents))))
		}
	)
)
