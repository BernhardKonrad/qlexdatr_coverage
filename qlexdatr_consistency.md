# qlexdatr_movement_consistency
Bernhard Konrad  
January 26, 2016  

**UPDATE 2016-02-13**

We want to run a data quality analysis on the `qlexdatr` dataset to check for movement consistencies. We're mainly following the analysis done earlier for `klesdatr`. Eventually these scripts could be combined to work on either data set. However, unfortunately the data sets are not yet as coherent as they would need for be.

Later, in our mathematical model, we assume that, **at each time point, a fish can either stay in the current lake section or move to an adjacent section**. The goal of this analysis is to compute how accurate this simplication is by comparing to the observed data.


## Preparation

First, let's grab the data


```r
library(qlexdatr)
library(lexr)
library(tidyr)
library(spa)
library(plyr)
library(dplyr)
library(magrittr)
library(ggplot2)
library(ggrepel)  # devtools::install_github("slowkow/ggrepel")


# WARNING: THIS WILL TAKE UP TO AN HOUR THE FIRST TIME IT'S RUN. THEN INSTANT, AS RESULTS ARE SAVED
load_or_generate_data <- function(){

  tryCatch(load(file = "QLEX_DETECT_DATA.RData", .GlobalEnv),
           error = function(e) {
             QLEX <- input_lex_data("qlexdatr")
             HOURLY_INTERVALS = c(1L, 2L, 4L, 6L, 12L)
             QLEX_DETECT_DATA <- list("1" = lexr::make_detect_data(QLEX, hourly_interval = 1L))
             for (hourly_interval in HOURLY_INTERVALS) {
               if (hourly_interval == 1) { next }
               else {
                 QLEX_DETECT_DATA[[as.character(hourly_interval)]] <-
                   lexr::make_detect_data(QLEX, hourly_interval = hourly_interval)
               }
             }
             save(QLEX_DETECT_DATA, file = "QLEX_DETECT_DATA.RData")
             
             #####
             #####
             #####
             
             SECTIONS_FROM_TO <- QLEX_DETECT_DATA[[1]]$distance
             
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
                   res[i] <- SECTIONS_FROM_TO %>%
                     filter(SectionFrom == sectionA[i], SectionTo == sectionB[i]) %>%
                     select(Distance) %>%
                     as.integer()
                 }
               }
               res
             }
             
             stopifnot(c(0, 0) == distance_sections(c("S01", "S01"), c("S01", "S01")))
             
             #####
             #####
             #####
             
             add_time_information <- function(detect_data){
               detect_data$detection %<>%
                 dplyr::left_join(detect_data$interval, by = c("IntervalDetection" = "Interval"))
               detect_data
             }
             
             for (i in 1:length(QLEX_DETECT_DATA)) {
               QLEX_DETECT_DATA[[i]] <- add_time_information(QLEX_DETECT_DATA[[i]])
             }
             
             #####
             #####
             #####
             
             add_consistency_score <- function(detect_data){
               detect_data$detection %<>%
                 group_by(Capture) %>%
                 mutate(TimeDiff = IntervalDetection - lag(IntervalDetection)) %>%
                 mutate(LocationDist = distance_sections(Section, lag(Section))) %>%
                 filter(!is.na(TimeDiff)) %>%
                 ungroup() %>%
                 mutate(is_consistent = LocationDist <= TimeDiff) %>%
                 rowwise %>%
                 mutate(inconsistent_val = max(-1, as.integer(LocationDist - TimeDiff)))
               
               detect_data
             }
             
             for (i in 1:length(QLEX_DETECT_DATA)) {
               QLEX_DETECT_DATA[[i]] <- add_consistency_score(QLEX_DETECT_DATA[[i]])
             }
             
             save(QLEX_DETECT_DATA, file = "QLEX_DETECT_DATA.RData")
             load(file = "QLEX_DETECT_DATA.RData", .GlobalEnv)
           })
}

load_or_generate_data()
```
  
> By extension of our assumption a movement of a fish is **consistent** if the number of time intervals between detections does not exceed the distance between the sections that the fish was detected in.


The variable `inconsistent_val` shows how inconsistent the fish movement is. If `inconsistent_val > 0` then the fish movement is inconsistent with our assumption that fish travel at most one section per day. In this case `inconsistent_val` says how many times, at least, the fish swam more than one section per day. E.g. if the fish was detected 3 sections away 2 days later, `inconsistent_val = 3-2 = 1`, if the fish was detected 10 sections away on the next day, then `inconsistent_val = 10-1 = 9`.

In contrast, `inconsistent_val = 0` indicated that the fish swam *on average* one section per day, that is, the fish was detected exactly $n$ sections aways $n$ days later, $n \geq 1$.

Finally, if the fish swam less than one section per day we set `inconsistent_val = -1`. That is, the fish was detected $k$ sections away $n$ days later, with $k < n$. We include the case where the fish was detected in the same section. In fact, this is the reason why we truncate at $-1$: Often fish don't seem to move at all, leading to large negative `inconsistent_val` values otherwise, which are not as informative. If a fish doesn't move it shouldn't matter if we detect it in the same section every day or not.

Note that `inconsistent_val` and `is_consistent` are both an underestimate of the error we make with our model assumption. It is possible that fish swim more than one section per day and are still listed as consistent, due to the non-perfect coverage. An inconsistent movement may not always be detected if the receiver in the violating section did not register the fast fish.


Here is a visualization of our result so far:


```r
QLEX_DETECT_DATA[[1]]$detection %>%
  ggplot(aes(x = as.integer(TimeDiff), y = LocationDist)) +
  geom_abline(intercept = 0.5, slope = 1, size = 2, color = "red") +
  stat_sum(aes(size = ..n..), geom = "point") +
  scale_size_area(max_size = 10) +
  xlim(0, 15) +
  ylim(0, 15) +
  annotate("text", label = "Inconsistent movement", x = 3, y = 8, size=5)
```

```
## Warning: Removed 17369 rows containing non-finite values (stat_sum).
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-2-1.png) 

We see that the most observations are fish which are detected in the same section on the next day.


Here are a few more metrics that may be useful

```r
plot_grand_mean_frac_consistent <- function(){
  lapply(X = QLEX_DETECT_DATA, FUN = function(x) {
    x$detection %>% 
      dplyr::summarise(frac_consistent = sum(is_consistent) / n()) %>%
      summarize(grand_mean = mean(frac_consistent))
  }) %>% 
    ldply(.fun = data.frame) %>% 
    transmute(time_diff_in_hours = as.numeric(.id),
              grand_mean_frac_consistent = grand_mean) %>% 
    ggplot(aes(x = time_diff_in_hours, y = grand_mean_frac_consistent)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    xlab("Length of time intervals (hours)") +
    ylab("Grand mean of fractions of consistent movements per fish")
}
plot_grand_mean_frac_consistent()
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-3-1.png) 

```r
#QLEX_DETECT_DATA[[1]]$detection %>%
#  ggplot(aes(x = consistent_val, fill = consistent)) +
#  geom_histogram(binwidth = 1)
```



```r
plot_frac_totally_consistent <- function(){
  lapply(X = QLEX_DETECT_DATA, FUN = function(x) {
  x$detection %>% 
    group_by(Capture) %>% 
    summarise(totally_consistent = all(is_consistent)) %>% 
    summarise(frac_totally_consistent = sum(totally_consistent) / n())
}) %>% 
  ldply(.fun = data.frame) %>% 
  transmute(time_diff_in_hours = as.numeric(.id),
            frac_totally_consistent = frac_totally_consistent ) %>% 
  ggplot(aes(x = time_diff_in_hours, y = frac_totally_consistent)) +
    geom_point(size = 5) +
    geom_line(size = 2) +
    xlab("Length of time intervals (hours)") +
    ylab("Fractions of fish that are totally consistent")
}

plot_frac_totally_consistent()
```

```
## Warning: Grouping rowwise data frame strips rowwise nature
```

```
## Warning: Grouping rowwise data frame strips rowwise nature
```

```
## Warning: Grouping rowwise data frame strips rowwise nature
```

```
## Warning: Grouping rowwise data frame strips rowwise nature
```

```
## Warning: Grouping rowwise data frame strips rowwise nature
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-4-1.png) 

We see that reducing the length of the time intervals has huge benefits to the consistency of the data.



## Understanding inconsistencies on a fish-by-fish basis

Next we want to better understand the inconsistent movements, where the largest `inconsistent_val` values are most concerning. To this end we add a couple consistency-related statistics to each fish.


```r
QLEX_DETECT_DATA[[1]]$detection %<>%
  group_by(Capture) %>%
  mutate(n_obs = length(DateTime),
         n_consistent = sum(is_consistent),
         frac_consistent = sum(is_consistent)/length(is_consistent),
         mean_inconsistent_val = mean(inconsistent_val),
         max_inconsistent_val = max(inconsistent_val),
         median_inconsistent_val = quantile(inconsistent_val, 0.5)
  ) %>%
  ungroup()
```

```
## Warning: Grouping rowwise data frame strips rowwise nature
```

```r
QLEX_DETECT_DATA[[1]]$detection %>%
  ggplot(aes(x = frac_consistent, y = reorder(Capture, frac_consistent))) +
  geom_point(aes(size = n_obs, color = max_inconsistent_val))
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-5-1.png) 

Or, plotting the same data sorted by their maximal consistent value:


```r
QLEX_DETECT_DATA[[1]]$detection %>%
  ggplot(aes(x = frac_consistent, y = reorder(Capture, max_inconsistent_val))) +
  geom_point(aes(size = n_obs, color = max_inconsistent_val))
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-6-1.png) 

For starters, we see that there is a number of fish who are always consistent.

### Fish with perfect consistency


```r
QLEX_DETECT_DATA[[1]]$detection %>%
  filter(frac_consistent == 1) %>%
  ggplot(aes(x = DateTime, y = as.numeric(Section))) +
  geom_line() +
  geom_point(size = 5) +
  facet_wrap(~ Capture)
```

```
## geom_path: Each group consists of only one observation. Do you need to
## adjust the group aesthetic?
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-7-1.png) 


### Fish with largest `inconsistent_val`, that is, the most inconsistent


```r
QLEX_DETECT_DATA[[1]]$detection %>%
  filter(max_inconsistent_val > max(QLEX_DETECT_DATA[[1]][["detection"]][["inconsistent_val"]]) - 2) %>%
  ggplot(aes(x = DateTime, y = as.numeric(Section))) +
   geom_line() +
   geom_point(aes(size = 2, colour = Hour)) +
   facet_wrap(~ Capture, ncol = 3) +
   theme(legend.position="none")
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-8-1.png) 


```r
bad_fish <- QLEX_DETECT_DATA[[1]]$detection %>%
  filter(max_inconsistent_val > max(QLEX_DETECT_DATA[[1]][["detection"]][["inconsistent_val"]] - 1)) %>% 
  select(Capture) %>% 
  unique()


plot_inconsistency_timeline_per_fish <- function(df, cap) {
  temp <- df %>%
    filter(Capture == cap) %>%
    rowwise() %>%
    mutate(my_str = if (inconsistent_val > 0) as.character(inconsistent_val) else "")

  temp %>%
    ggplot(aes(x = DateTime, y = as.numeric(Section))) +
    geom_line(aes(size = inconsistent_val + 1)) +
    geom_point(aes(size = 2)) +
    geom_text_repel(data = subset(temp, my_str != ""), aes(label = my_str), size = 8, colour = "red") +
    #geom_text(data = subset(temp, my_str != ""), aes(label = my_str), size = 8) +
    ggtitle(cap) +
    theme(legend.position="none")
}

apply(bad_fish, 1, plot_inconsistency_timeline_per_fish, df = QLEX_DETECT_DATA[[1]]$detection)
```

```
## [[1]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-1.png) 

```
## 
## [[2]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-2.png) 

```
## 
## [[3]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-3.png) 

```
## 
## [[4]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-4.png) 

```
## 
## [[5]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-5.png) 

```
## 
## [[6]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-6.png) 

```
## 
## [[7]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-7.png) 

```
## 
## [[8]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-8.png) 

```
## 
## [[9]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-9.png) 

```
## 
## [[10]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-10.png) 

```
## 
## [[11]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-11.png) 

```
## 
## [[12]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-12.png) 

```
## 
## [[13]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-13.png) 

```
## 
## [[14]]
```

![](qlexdatr_consistency_files/figure-html/unnamed-chunk-9-14.png) 



Finally, we're writing our results to a csv file.


```r
QLEX_DETECT_DATA[[1]]$detection %>% 
  group_by(Capture) %>% 
  summarize(n_obs = n(),
            n_consistent = sum(is_consistent),
            n_inconsistent = sum(!is_consistent),
            frac_consistent = min(frac_consistent),
            largest_jump = max(inconsistent_val)
  ) %>% 
  join(
    QLEX_DETECT_DATA[[1]]$detection %>% 
      group_by(Capture) %>% 
      filter(TimeDiff < 24) %>% 
      summarize(max_jump_in_24h = max(LocationDist)) %>% 
      join(qlexdatr::capture %>% 
             select(Capture, Species))
  ) %>% 
  write.csv(file = "overview_per_fish_qlexdatr.csv",
            row.names = FALSE,
            quote = FALSE)
```

```
## Joining by: Capture
## Joining by: Capture
```
