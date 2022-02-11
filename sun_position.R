


library(lubridate)
library(tidyverse)
library(ggplot2)

angle_convert <- function(angle, to = "radians") {
        if (to == "radians") {
                return(angle * pi / 180)
                
        } else if (to == "degrees") {
                return(angle * 180 / pi)
                
        }
        
}

solar_declination <- function(day_num) {
        angle_rad <- 360 / 365 * (day_num + 284) * pi / 180
        
        23.45 * sin(angle_rad)
        
}

solar_time <- function(day_hour) {
        (day_hour - 12) / 24 * 360
        
}


solar_zenith <-
        function(declination_deg,
                 latitude_deg,
                 solar_time_deg) {
                declination_rad <- declination_deg |> angle_convert(to = "radians")
                latitude_rad <-
                        latitude_deg |> angle_convert(to = "radians")
                solar_time_rad <-
                        solar_time_deg |> angle_convert(to = "radians")
                
                cos_zenith <-
                        sin(declination_rad) * sin(latitude_rad) + cos(declination_rad) * cos(latitude_rad) * cos(solar_time_rad)
                
                return(acos(cos_zenith))
        }

solar_altitude <-
        function(declination_deg,
                 latitude_deg,
                 solar_time_deg) {
                declination_rad <- declination_deg |> angle_convert(to = "radians")
                latitude_rad <-
                        latitude_deg |> angle_convert(to = "radians")
                solar_time_rad <-
                        solar_time_deg |> angle_convert(to = "radians")
                
                sin_altitude <-
                        sin(declination_rad) * sin(latitude_rad) + cos(declination_rad) * cos(latitude_rad) * cos(solar_time_rad)
                
                return(asin(sin_altitude))
        }

solar_azimuth <-
        function(solar_altitude_rad,
                 latitude_deg,
                 solar_declination_deg) {
                latitude_rad <- latitude_deg |> angle_convert(to = "radians")
                solar_declination_rad <-
                        solar_declination_deg |> angle_convert(to = "radians")
                
                
                cos_azimuth <-
                        (sin(solar_altitude_rad) * sin(latitude_rad) - sin(solar_declination_rad)) * sign(latitude_deg) / (cos(solar_altitude_rad) * cos(latitude_rad))
                
                return(acos(cos_azimuth))
        }

solar_position_calculator <- function(latitude_deg, gcr = 0.3) {
        analysis_df <- tibble(timestamp = seq(
                ymd_hm("1990-01-01 00:00"),
                ymd_hm("1991-01-01 00:00"),
                by = "5 mins"
        )) %>%
                mutate(
                        dia = yday(timestamp),
                        hor = hour(timestamp) + minute(timestamp) / 60,
                        solar_time_deg = solar_time(hor),
                        solar_time_rad = angle_convert(solar_time_deg,
                                                       to = "radians"),
                        solar_declination_deg = solar_declination(dia),
                        solar_zenith_rad = solar_zenith(
                                declination_deg = solar_declination_deg,
                                latitude_deg = latitude_deg,
                                solar_time_deg = solar_time_deg
                        ),
                        solar_zenith_deg = angle_convert(solar_zenith_rad,
                                                         to = "degrees"),
                        solar_altitude_rad = solar_altitude(
                                declination_deg = solar_declination_deg,
                                latitude_deg = latitude_deg,
                                solar_time_deg = solar_time_deg
                        ),
                        solar_altitude_deg = angle_convert(solar_altitude_rad,
                                                           to = "degrees"),
                        solar_azimuth_rad = solar_azimuth(
                                solar_altitude_rad = solar_altitude_rad,
                                latitude_deg = latitude_deg,
                                solar_declination_deg = solar_declination_deg
                        ) * sign(solar_time_deg),
                        solar_azimuth_deg = angle_convert(solar_azimuth_rad,
                                                          to = "degrees"),
                        optimal_tracking_angle_rad = if_else(solar_altitude_deg >= 0, 
                                                         optimal_tracking_angle(solar_altitude = solar_altitude_rad,
                                                                        solar_azimuth = solar_azimuth_rad),
                                                         0),
                        optimal_tracking_angle_deg = angle_convert(optimal_tracking_angle_rad, 
                                                                   to = "degrees"),
                        backtracking_correction_angle_rad = backtracking_correction_angle(optimal_tracking_angle_rad,
                                                                                      gcr = gcr) * -sign(solar_time_deg),
                        backtracking_correction_angle_deg = angle_convert(backtracking_correction_angle_rad,
                                                                          to = "degrees") %>%
                                replace_na(replace = 0),
                        backtracking_angle_deg = backtracking_correction_angle_deg + optimal_tracking_angle_deg,
                                                        
                )
        
        return(analysis_df)
}

position_plotter_usual <- function(latitude_deg, custom_yday) {
        special_data <- solar_position_calculator(latitude_deg = latitude_deg) %>%
                mutate(special_date = if_else(
                        solar_declination_deg == max(solar_declination_deg, 
                                                     na.rm = T),
                        "Summer Solstice",
                        if_else(
                                round(solar_declination_deg, 2) == 11.75,
                                "Half-Point Summer",
                                if_else(
                                        round(solar_declination_deg, 2) == 0,
                                        "Fall Equinox",
                                        if_else(
                                                round(solar_declination_deg, 2) == -11.75,
                                                "Half-Point Fall",
                                                if_else(
                                                        solar_declination_deg == min(solar_declination_deg, 
                                                                                     na.rm = T),
                                                        "Winter Solstice",
                                                        if_else(dia == custom_yday,
                                                                "Selected Day",
                                                                "Not Special")
                                                )
                                        )
                                )
                        )
                )) %>%
                filter(special_date != "Not Special") 
        
        
        p <- special_data %>%
                ggplot(
                        aes(x = solar_azimuth_deg, 
                            y = solar_altitude_deg, 
                            color = special_date)) + 
                geom_point() + 
                ylim(c(0, 90)) +
                xlim(c(-180, 180))
        
        p + ylab("Solar Elevation [Degrees]") +
                xlab("Solar Azimuth [Degrees]") +
                ggtitle("Solar Trajectories Across Selected Dates",
                        subtitle = paste("Calculated for latitude =",
                                         latitude_deg, "degrees"))
        
        
}

optimal_tracking_angle <- function(solar_altitude, solar_azimuth) {
        
        
        x <- cos(solar_altitude) * sin(solar_azimuth) 
        z <- sin(solar_altitude)
        
        tan_omega_ideal <- x / z
        
        return(atan(tan_omega_ideal))
        
}

backtracking_correction_angle <- function(optimal_tracking_angle_rad, gcr) {
        
        cos_theta_c <- cos(optimal_tracking_angle_rad) / gcr
        
        return(acos(cos_theta_c))
        
}