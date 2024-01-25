# ZollnerGravity
Zollner Pendulum
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 22:17:21 2024

@author: jerom
"""

# -*- coding: utf-8 -*-
"""
Created on Wed Jan 24 21:57:13 2024

@author: jerom
"""

from astropy.coordinates import SkyCoord, EarthLocation, AltAz, get_moon, get_sun
from astropy.time import Time, TimeDelta
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt
import datetime

# Constants
G = 6.67430e-11  # Gravitational constant, N(m^2)/(kg^2)
mass_pendulum = 0.5  # Mass of the pendulum in kg
mass_moon = 7.342e22  # Mass of the Moon in kg
mass_sun = 1.9885e30  # Mass of the Sun in kg
radius_earth = 6371.0  # average radius of the Earth in km

# Function to calculate gravitational force
def gravitational_force(m1, m2, r):
    return G * (m1 * m2) / r**2

# Function to convert RA and Dec to Altitude and Azimuth
def ra_dec_to_alt_az(ra, dec, observer_lat, observer_lon, time):
    # Define the coordinates in the sky
    sky_coord = SkyCoord(ra=ra * u.deg, dec=dec * u.deg, frame='icrs')

    # Define the location of the observer
    observer_location = EarthLocation(lat=observer_lat * u.deg, lon=observer_lon * u.deg)

    # Define the time of observation (UTC)
    observation_time = Time(time)

    # Define the AltAz frame at the observation time and location
    altaz_frame = AltAz(obstime=observation_time, location=observer_location)

    # Transform the coordinates to the AltAz frame
    altaz_coord = sky_coord.transform_to(altaz_frame)

    # Get the Altitude and Azimuth
    altitude = altaz_coord.alt.deg
    azimuth = altaz_coord.az.deg

    return altitude, azimuth

# Prescott  Location
observer_location = EarthLocation(lat=34.54 * u.deg, lon=-112.47 * u.deg)

# Define time range
start_time = Time('2024-01-24T00:00:00')  # Start of 2024
end_time = Time('2024-01-26T23:59:59')  # End of May 2024
time_interval = TimeDelta((1/24) * u.day)  # Interval of 1 day

times = []
forces_ew = []

current_time = start_time
while current_time <= end_time:
    # Get the position of the Moon at the current time
    moon = get_moon(current_time, observer_location)

    # Transform moon position to AltAz frame
    altaz_frame = AltAz(obstime=current_time, location=observer_location)
    moon_altaz = moon.transform_to(altaz_frame)

    # Get altitude and azimuth from the AltAz object
    moon_altitude = moon_altaz.alt.deg
    moon_azimuth = moon_altaz.az.deg

    # Geocentric distance to the Moon
    distance_to_moon_geo = moon.distance.to(u.km).value
    
    # Distance from the center of the Earth to Prescott
    distance_earth_prescott = radius_earth + observer_location.height.to(u.km).value
    
    # Spherical law of cosines
    angle = np.radians(90 - moon_altitude)  # 90 degrees - altitude of the Moon
    
    # Law of cosines to find the distance from Prescott to the Moon
    distance_prescott_moon = np.sqrt(distance_to_moon_geo**2 + distance_earth_prescott**2 -
                                     2 * distance_to_moon_geo * distance_earth_prescott * np.cos(angle))

    # Calculate gravitational force using the distance from Prescott to the Moon
    force = gravitational_force(mass_pendulum, mass_moon, distance_prescott_moon * 1000)  # Convert km to m
    
    # Convert azimuth from degrees from North to radians for calculation
    azimuth_radians = np.radians(moon_azimuth)

    # Calculate force components in the horizontal plane (NS and EW)
    force_ns = force * np.cos(azimuth_radians) * np.cos(np.radians(moon_altitude))
    force_ew = force * np.sin(azimuth_radians) * np.cos(np.radians(moon_altitude))

    times.append(current_time.iso)
    forces_ew.append(force_ew)
    print(f"Time: {current_time.iso}, Distance: {distance_prescott_moon:.2f} km, Force_NS: {force_ns:.6e} N, Force_EW: {force_ew:.6e} N")

    # Move to the next time point
    current_time += time_interval
    
    
# ===Plotting===

# Parse and reformat dates to show only day and month
# formatted_dates = [datetime.datetime.strptime(date, '%Y-%m-%d %H:%M:%S.%f').strftime('%m-%d') for date in times]

plt.figure(figsize=(10, 6))
# plt.plot(formatted_dates, forces_ew, label='Force_EW', color='blue')
plt.plot(times, forces_ew, label='Force_EW', color='blue')
plt.xlabel('Time (2024)')
plt.ylabel('East-West Force (N)')
plt.title('East-West Gravitational Force on the Pendulum from the Moon (2024)')

# Set x-ticks to every 5th day
# date_ticks = [formatted_dates[i] for i in range(0, len(formatted_dates), 5)]  # Get every 5th date
date_ticks = [times[i] for i in range(0, len(times), 5)]  # Get every 5th date
plt.xticks(date_ticks, rotation=45)  # Set x-ticks and rotate for better readability

plt.tight_layout()  # Adjust layout
plt.grid(True)  # Add grid
plt.legend()
plt.show()
