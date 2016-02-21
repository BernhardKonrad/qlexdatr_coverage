library(lexr)
library(dplyr)


# Preparation
HOURLY_INTERVAL = 1L
QLEX <- lexr::input_lex_data("qlexdatr")
# Most time is used by lexr::make_detect_data
QLEX_DETECT <- lexr::make_detect_data(QLEX, hourly_interval = HOURLY_INTERVAL)


ptm <- proc.time()
# The actual function
add_consistency_info_by_fish <- function(detect_data){
  # Time information
  detect_data$detection <- detect_data$detection %>%
    dplyr::left_join(detect_data$interval, by = c("IntervalDetection" = "Interval"))

  # Lookup matrix for distance between sections
  SECTIONS_FROM_TO_MATRIX = matrix(0, nrow = 32, ncol = 32)
  for (i in 1:32){
    for (j in 1:32){
      SECTIONS_FROM_TO_MATRIX[i, j] <- detect_data$distance %>%
        dplyr::filter(as.integer(SectionFrom) == i,
                      as.integer(SectionTo) == j) %>%
        dplyr::select(Distance) %>%
        as.integer()
    }
  }

  # Helper function for distance between sections
  get_distance <- function(A, B) {
    res = vector(mode = "numeric", length = length(A))
    for (i in seq(1:length(A))){
      res[i] <- SECTIONS_FROM_TO_MATRIX[as.integer(A[i]), as.integer(B[i])]
    }
    res
  }

  # Calculate TimeDiff, LocationDist, and inconsistency_val
  detect_data$detection <- detect_data$detection %>%
    dplyr::group_by(Capture) %>%
    dplyr::mutate(TimeDiff = IntervalDetection - lag(IntervalDetection)) %>%
    dplyr::mutate(LocationDist = get_distance(Section, lag(Section))) %>%
    dplyr::ungroup() %>%
    dplyr::rowwise() %>%
    dplyr::mutate(inconsistent_val = max(0, as.integer(LocationDist - TimeDiff)))

  # Add simple summary statistics for each fish
  detect_data$detection <- detect_data$detection %>%
    dplyr::ungroup() %>%
    dplyr::group_by(Capture) %>%
    dplyr::mutate(n_obs = n(),
                     frac_inconsistent = sum(inconsistent_val > 0,
                                             na.rm = TRUE)/n_obs,
                     largest_inconsistent_val = max(inconsistent_val,
                                                    na.rm = TRUE)) %>%
    dplyr::ungroup()

  detect_data
}

QLEX_DETECT_with_consistency_info_by_fish <- add_consistency_info_by_fish(QLEX_DETECT)
print(proc.time() - ptm)
