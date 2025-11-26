library(readr)



# Create a sample CSV file
csv_data <- "datetime\n2001-10-10 20:10"
readr::write_file(csv_data, "datetimes.csv")

# By default, readr will interpret the time in UTC
read_csv("datetimes.csv")
#> # A tibble: 1 x 1
#>   datetime           
#>   <dttm>             
#> 1 2001-10-10 20:10:00

# To specify a different time zone, use the locale() function
read_csv("datetimes.csv", locale = locale(tz = "Pacific/Auckland"))
#> # A tibble: 1 x 1
#>   datetime           
#>   <dttm>             
#> 1 2001-10-10 20:10:00

# The output for Pacific/Auckland shows the correct time zone with Daylight Saving Time
# which is "NZDT" (New Zealand Daylight Time)
