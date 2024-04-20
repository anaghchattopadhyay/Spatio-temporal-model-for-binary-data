# reading the data and reordering as required
ny_data = read.csv(file.choose())
ny_data$week = as.Date(ny_data$week, format = "%Y-%m-%d")
dff = ny_data
dff = cbind(dff, time = (as.numeric(as.POSIXct(dff$week)) - min(as.numeric(as.POSIXct(dff$week)))) / (86400 * 7))

dff = as.data.frame(dff)
dff = dff %>% group_by(lat, long) %>%
  mutate(
    lagged_death = lag(log(1 + death), default = 0),
    linear_trend = c(1:n()),
    quadratic_trend = c(1:n()) ^ 2,
    sine_trend = sin(c(1:n())),
    cosine_trend = cos(c(1:n()))
  )

# taking necessary keys(key here refers to the different sub-parts in the given area)
Loc_unique = unique(dff$id)

# finding unique sub-locations
Loc_unique = as.vector(Loc_unique)

# finding number of locations
num_loc = length(Loc_unique)

# fix the regression formula for GLM and SAR
reg_f <- "covid_status ~ prop_vaccinated + I(log(death+1)) + lagged_death + I(log(population)) +
linear_trend + quadratic_trend + sine_trend + cosine_trend"
