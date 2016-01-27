library(qlexdatr)
library(lexr)
library(tidyr)
library(spa)
library(dplyr)
library(ggplot2)

QLEX <- input_lex_data("qlexdatr")
QLEX_DETECT_DATA <- lexr::make_detect_data(QLEX)
SECTIONS_FROM_TO <- QLEX_DETECT_DATA$distance

distance_sections <- function(sectionA, sectionB) {
  # Returns pairwise distance between sections
  sectionA <- as.character(sectionA)
  sectionB <- as.character(sectionB)
  stopifnot(length(sectionA) == length(sectionB))

  res <- vector(mode="integer", length=length(sectionA))

  for (i in 1:length(res)) {
    if (is.na(sectionA[i]) | is.na(sectionB[i])) {
      res[i] = NA
    } else {
      res[i] <-   SECTIONS_FROM_TO %>%
        filter(SectionFrom == sectionA[i], SectionTo == sectionB[i]) %>%
        select(Distance) %>%
        as.integer()
    }
  }
  res
}

stopifnot(0 == distance_sections("S01", "S01"))
stopifnot(c(0, 0) == distance_sections(c("S01", "S01"), c("S01", "S01")))


QLEX_DETECTION <- QLEX_DETECT_DATA$detection %>%
  dplyr::left_join(QLEX_DETECT_DATA$interval, by = c("IntervalDetection" = "Interval"))

QLEX_DETECTION %>%
  filter(as.character(Capture) == "F001") %>%
  ggplot(aes(x = DateTime, y = as.numeric(Section))) +
  geom_point() +
  geom_line()

glimpse(QLEX_DETECTION)

#  mutate(TimeDiff = DayteTime - lag(DayteTime)) %>%


QLEX_DETECTION %<>%
#  filter(as.character(Capture) < "F003") %>%
#  droplevels() %>%
  group_by(Capture) %>%
  mutate(TimeDiff = IntervalDetection - lag(IntervalDetection)) %>%
  mutate(LocationDist = distance_sections(Section, lag(Section))) %>%
  filter(!is.na(TimeDiff)) %>%
  ungroup() %>%
  mutate(consistent = LocationDist <= TimeDiff) %>%
  rowwise %>%
  mutate(consistent_val = max(-1, as.integer(LocationDist - TimeDiff)))

glimpse(QLEX_DETECTION)


QLEX_DETECTION %>%
  ggplot(aes(x = as.integer(TimeDiff), y = LocationDist)) +
  geom_abline(intercept = 0, slope = 1, size = 2, color = "red") +
  stat_sum(aes(size = ..n..), geom = "point") +
  scale_size_area(max_size = 10) +
  xlim(0, 15) +
  ylim(0, 15) +
  annotate("text", label = "Inconsistent movement", x = 3, y = 8, size=5)



QLEX_DETECTION %>%
  ggplot(aes(x = DateTime, y = as.numeric(Section))) +
  geom_line() +
  geom_point(aes(color = consistent), size = 5) +
  facet_wrap(~ Capture)

QLEX_DETECTION %>%
  summarise(frac_consistent = sum(consistent) / n()) %>%
  summarize(grand_mean = mean(frac_consistent))

QLEX_DETECTION %>%
  ggplot(aes(x = consistent_val, fill = consistent)) +
  geom_histogram(binwidth = 1)


QLEX_DETECTION %<>%
  group_by(Capture) %>%
  mutate(n_obs = length(DateTime),
         n_consistent = sum(consistent),
         perc_consistent = sum(consistent)/length(consistent),
         mean_consistent_val = mean(consistent_val),
         max_consistent_val = max(consistent_val),
         median_consistent_val = quantile(consistent_val, 0.5)
  ) %>%
  ungroup()

QLEX_DETECTION %>%
  ggplot(aes(x = perc_consistent, y = reorder(Capture, perc_consistent))) +
  geom_point(aes(size = n_obs, color = max_consistent_val))

QLEX_DETECTION %>%
  ggplot(aes(x = perc_consistent, y = reorder(Capture, max_consistent_val))) +
  geom_point(aes(size = n_obs, color = max_consistent_val))

QLEX_DETECTION %>%
  filter(perc_consistent == 1) %>%
  ggplot(aes(x = DateTime, y = as.numeric(Section))) +
  geom_line() +
  geom_point(size = 5) +
  facet_wrap(~ Capture)

#We see that these `r QLEX_DETECTION %>% filter(perc_consistent == 1) %>% nrow()` fish are perfectly consistent, great!


### Fish with largest `consistent_val`, that is, the most inconsistent

QLEX_DETECTION %>%
  filter(max_consistent_val >= max(QLEX_DETECTION$consistent_val) - 2) %>%
  ggplot(aes(x = DateTime, y = as.numeric(Section))) +
  geom_line(aes(size = consistent_val)) +
  geom_point(aes(size = 5, color = DateTime)) +
  facet_wrap(~ Capture, ncol = 2)


QLEX_DETECTION %>%
  filter(max_consistent_val >= max(QLEX_DETECTION$consistent_val) - 2) %>%
  rowwise() %>%
  mutate(my_str = if (consistent_val > 0) as.character(consistent_val) else "") %>%
  ggplot(aes(x = DateTime, y = as.numeric(Section))) +
    geom_line(aes(size = consistent_val)) +
    geom_point(aes(size = 5)) +
    geom_text(aes(label = my_str), hjust = 1.5, size = 8) +
    facet_wrap(~ Capture)
