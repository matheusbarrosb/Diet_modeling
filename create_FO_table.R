require(data.table)
require(flextable)
require(dplyr)
require(rlang)


DT <- as.data.frame(cbind(LAGRHO$Site, LAGRHO$Treatment, LAGRHO$SAV, LAGRHO$Amphipod, LAGRHO$Crustacean))
names(DT) <- c("site", "treat", "SAV", "amphipod", "crustacean")

# main function
xtabs3 <- function(data,
                   x,
                   y,
                   z) {
  
  # internal helper function, converts to factor
  not_a_factor <- function(x){
    !is.factor(x)
  }
  
  # capture variable names
  xlab <- rlang::as_name(rlang::enquo(x))
  ylab <- rlang::as_name(rlang::enquo(y))
  zlab <- rlang::as_name(z)
  
  # create temp local dataframe 
  data <-
    dplyr::select(
      .data = data,
      x = {{ x }},
      y = {{ y }},
      z = {{ z }}
    )
  
  # calculate counts and percents 
  
  # x, y and z need to be a factor or ordered factor
  # also drop the unused levels of the factors and NAs
  data <- data %>%
    dplyr::mutate_if(.tbl = ., not_a_factor, as.factor) %>%
    dplyr::mutate_if(.tbl = ., is.factor, droplevels) %>%
    dplyr::filter_all(.tbl = ., all_vars(!is.na(.))) %>%
    dplyr::as_tibble(x = .)
  
  # convert the data into percentages; group by x, y, z
  # DO NOT Drop zeroes  
  df <-
    data %>%
    dplyr::group_by(.data = ., x, y, z, .drop = TRUE) %>%
    dplyr::summarize(.data = ., counts = n()) %>%
    dplyr::mutate(.data = ., perc = (counts / sum(counts)) * 100) %>%
    dplyr::ungroup(x = .) %>%
    rename(!!xlab := x, !!ylab := y, "level" := z)
  
  return(df)
  
}

# Make a list of all the binary variables we want to use
# best if it's a named list variables can be bare or quoted
fff <- alist(SAV = SAV, amph = amphipod, crust = crustacean)


xxx <- purrr::map_dfr(.x = fff, ~ xtabs3(DT, site, treat, .x), .id = "prey")

xxx <- xxx[xxx$level != 0,]

level <- rep(1, length(xxx$prey))

xxx$level <- as.factor(level)

myft <- flextable(xxx, col_keys = c("prey", "treat", "site", "perc"))
myft <- theme_vanilla(myft)
myft <- merge_v(myft, j = c("treat", "site", "prey") )
myft <- autofit(myft)
myft <- colformat_num(x = myft, j = c("perc"), digits = 1, suffix = "%", na_str = "0")
# reprex won't let me make an html table
plot(myft)


