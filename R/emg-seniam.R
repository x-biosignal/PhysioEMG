# Normalize a muscle name for matching (lower-case, letters only).
.seniam_norm <- function(s) gsub("[^a-z]", "", tolower(as.character(s)))

# Resolve a muscle name against the SENIAM catalog. Returns the catalog row
# index (or NA) plus a vector of near matches for messaging.
.seniam_resolve <- function(muscle, catalog) {
  idx <- match(.seniam_norm(muscle), .seniam_norm(catalog$muscle))
  near <- if (is.na(idx)) {
    agrep(muscle, catalog$muscle, ignore.case = TRUE, value = TRUE,
          max.distance = 0.35)
  } else {
    character(0)
  }
  list(index = idx, near = near)
}

.seniam_resolve_channels <- function(x, channel, n_channels) {
  if (is.null(channel)) return(seq_len(n_channels))
  if (is.numeric(channel)) {
    idx <- as.integer(channel)
  } else {
    labels <- tryCatch(as.character(SummarizedExperiment::colData(x)$label),
                       error = function(e) NULL)
    idx <- match(as.character(channel), labels)
  }
  if (any(is.na(idx)) || any(idx < 1L) || any(idx > n_channels)) {
    stop("channel out of range or not found", call. = FALSE)
  }
  idx
}

#' SENIAM Muscle Catalog
#'
#' Returns the bundled catalog of SENIAM-recommended surface-EMG electrode
#' locations (Hermens et al. 2000): for each muscle the recommended sensor
#' location, the electrode orientation (fibre direction) and the recommended
#' inter-electrode distance (mm).
#'
#' @return A data.frame with columns \code{muscle}, \code{location},
#'   \code{orientation} and \code{ied_mm}.
#' @seealso [setEMGElectrode()] to attach electrode metadata to a
#'   PhysioExperiment, [getEMGElectrode()] to retrieve it
#' @references Hermens, H.J., Freriks, B., Disselhorst-Klug, C. & Rau, G. (2000).
#'   "Development of recommendations for SEMG sensors and sensor placement
#'   procedures." Journal of Electromyography and Kinesiology, 10(5), 361-374.
#'   doi:10.1016/S1050-6411(00)00027-4
#' @export
#' @examples
#' head(seniamMuscles())
seniamMuscles <- function() {
  path <- system.file("extdata", "seniam_muscles.csv", package = "PhysioEMG")
  utils::read.csv(path, stringsAsFactors = FALSE)
}

#' Set SENIAM Electrode Metadata
#'
#' Attaches standardized SENIAM electrode/muscle metadata to the channels of a
#' PhysioExperiment, written as \code{seniam_*} columns of \code{colData}.
#' Muscle names are validated against [seniamMuscles()]; a recognized muscle
#' also back-fills its recommended electrode location, orientation and
#' inter-electrode distance when those are not supplied.
#'
#' @param x A PhysioExperiment object with EMG channels.
#' @param muscle Muscle name(s) (character). Length 1 (applied to all targeted
#'   channels) or one per targeted channel.
#' @param placement Optional electrode location text; back-filled from the
#'   catalog when \code{NULL}.
#' @param ied_mm Optional inter-electrode distance in mm; back-filled from the
#'   catalog (SENIAM default 20 mm) when \code{NULL}.
#' @param reference Optional reference-electrode description.
#' @param side Optional body side ("left"/"right"/"L"/"R"/...).
#' @param channel Channels to annotate: integer indices, channel labels, or
#'   \code{NULL} (default) for all channels.
#' @return The PhysioExperiment with \code{seniam_muscle},
#'   \code{seniam_placement}, \code{seniam_orientation}, \code{seniam_ied_mm},
#'   \code{seniam_reference} and \code{seniam_side} added to \code{colData}.
#' @seealso [getEMGElectrode()], [seniamMuscles()], [emgMFCV()] which can use the
#'   stored inter-electrode distance
#' @references Hermens, H.J. et al. (2000). "Development of recommendations for
#'   SEMG sensors and sensor placement procedures." Journal of Electromyography
#'   and Kinesiology, 10(5), 361-374. doi:10.1016/S1050-6411(00)00027-4
#' @export
#' @examples
#' pe <- PhysioExperiment(
#'   assays = list(raw = matrix(rnorm(200), ncol = 2)),
#'   colData = S4Vectors::DataFrame(label = c("ch1", "ch2")),
#'   samplingRate = 1000)
#' pe <- setEMGElectrode(pe, muscle = c("Biceps brachii", "Triceps brachii"),
#'                       side = "right")
#' getEMGElectrode(pe)
setEMGElectrode <- function(x, muscle, placement = NULL, ied_mm = NULL,
                            reference = NULL, side = NULL, channel = NULL) {
  stopifnot(inherits(x, "PhysioExperiment"))
  n_channels <- ncol(SummarizedExperiment::assay(x, defaultAssay(x)))
  ch <- .seniam_resolve_channels(x, channel, n_channels)
  m <- length(ch)

  recycle <- function(v, what) {
    if (is.null(v)) return(rep(NA, m))
    if (length(v) == 1L) return(rep(v, m))
    if (length(v) != m) {
      stop(sprintf("'%s' must have length 1 or %d (targeted channels)",
                   what, m), call. = FALSE)
    }
    v
  }
  muscle <- recycle(muscle, "muscle")
  placement <- recycle(placement, "placement")
  ied_mm <- recycle(ied_mm, "ied_mm")
  reference <- recycle(reference, "reference")
  side <- recycle(side, "side")

  catalog <- seniamMuscles()
  canon <- character(m); orient <- rep(NA_character_, m)
  loc <- as.character(placement); ied <- as.numeric(ied_mm)
  for (i in seq_len(m)) {
    if (is.na(muscle[i])) { canon[i] <- NA_character_; next }
    res <- .seniam_resolve(muscle[i], catalog)
    if (is.na(res$index)) {
      canon[i] <- as.character(muscle[i])
      msg <- sprintf("muscle '%s' is not in the SENIAM catalog", muscle[i])
      if (length(res$near)) {
        msg <- paste0(msg, "; did you mean: ",
                      paste(utils::head(res$near, 5), collapse = ", "), "?")
      } else {
        msg <- paste0(msg, "; see seniamMuscles() for valid names")
      }
      warning(msg, call. = FALSE)
    } else {
      row <- catalog[res$index, ]
      canon[i] <- row$muscle
      orient[i] <- row$orientation
      if (is.na(loc[i])) loc[i] <- row$location
      if (is.na(ied[i])) ied[i] <- row$ied_mm
    }
  }

  cd <- SummarizedExperiment::colData(x)
  ensure <- function(name, default) {
    if (is.null(cd[[name]])) cd[[name]] <<- rep(default, nrow(cd))
  }
  ensure("seniam_muscle", NA_character_)
  ensure("seniam_placement", NA_character_)
  ensure("seniam_orientation", NA_character_)
  ensure("seniam_ied_mm", NA_real_)
  ensure("seniam_reference", NA_character_)
  ensure("seniam_side", NA_character_)

  cd$seniam_muscle[ch] <- canon
  cd$seniam_placement[ch] <- loc
  cd$seniam_orientation[ch] <- orient
  cd$seniam_ied_mm[ch] <- ied
  cd$seniam_reference[ch] <- as.character(reference)
  cd$seniam_side[ch] <- as.character(side)
  SummarizedExperiment::colData(x) <- cd
  x
}

#' Get SENIAM Electrode Metadata
#'
#' Retrieves the SENIAM electrode/muscle metadata previously attached with
#' [setEMGElectrode()].
#'
#' @param x A PhysioExperiment object.
#' @return A data.frame with one row per channel and columns \code{channel},
#'   \code{muscle}, \code{placement}, \code{orientation}, \code{ied_mm},
#'   \code{reference} and \code{side} (\code{NA} where unset).
#' @seealso [setEMGElectrode()], [seniamMuscles()]
#' @export
#' @examples
#' pe <- PhysioExperiment(assays = list(raw = matrix(rnorm(100), ncol = 1)),
#'                        samplingRate = 1000)
#' pe <- setEMGElectrode(pe, muscle = "Tibialis anterior")
#' getEMGElectrode(pe)
getEMGElectrode <- function(x) {
  stopifnot(inherits(x, "PhysioExperiment"))
  cd <- SummarizedExperiment::colData(x)
  n <- nrow(cd)
  col <- function(name, default) {
    if (is.null(cd[[name]])) rep(default, n) else cd[[name]]
  }
  data.frame(
    channel = seq_len(n),
    muscle = col("seniam_muscle", NA_character_),
    placement = col("seniam_placement", NA_character_),
    orientation = col("seniam_orientation", NA_character_),
    ied_mm = col("seniam_ied_mm", NA_real_),
    reference = col("seniam_reference", NA_character_),
    side = col("seniam_side", NA_character_),
    stringsAsFactors = FALSE)
}
