#!/usr/bin/env python
# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
from ipywidgets import interact, FloatSlider, IntSlider, SelectionSlider
from IPython.display import HTML

#constants
sigma = 5.67*10**(-8)  # W/m^2/K^4
albedo = 0.3
S_true = 1360          # W/m^2
T_obs = 288            # K
lat_Seattle = 47.6061  # degrees


def emission_temp(S=S_true, a=albedo):
    """
    Calculates emission temperature of a planet given its:
    S, total solar irradiance (tsi)
    a, albedo: a measure of reflectivity
    """
    T = ((S/4)*(1-a)/sigma)**0.25
    
    return T

def energy_budget(S=S_true, a=albedo):
    """
    Plots interactive energy budget bar graph, given its:
    S, total solar irradiance (tsi)
    a, albedo: a measure of reflectivity
    """
    T = emission_temp(S, a)
    print(f"Emission temperature: {T:.2f} K ({T - 273.15:.2f} °C)")
    print()
    print(f"Observed global mean temperature: {T_obs} K ({T_obs - 273.15:.2f} °C)")

    #plotting
    labels = ['Reflected Solar', 'Absorbed Solar', 'Outgoing Longwave']
    reflected = a * S / 4
    absorbed = (1 - a) * S / 4
    outgoing = sigma * T**4

    vals = [reflected, absorbed, outgoing]
    cols = ['gold', 'cornflowerblue', 'indianred']
    actual = [100,240,239]     # from textbook Fig. 2.4

    plt.figure(figsize=(6, 4))
    plt.bar(labels, vals, color=cols, width=0.8)
    plt.bar(labels, actual, color=cols, alpha=0.5, width=0.9)
    plt.title("Energy Flow in W/m$^2$, with actual values shaded")
    plt.ylim(0, 500)
    plt.grid()
    plt.show();
    
def show_animation(url):
    """
    Displays animation found at input url
    """
    display(HTML(f"""
    <video width=800 controls autoplay loop muted>
        <source src="{url}" type="video/mp4">
        Your browser does not support the video tag.
    </video>
    """))
    
def avg_daily_insolation(S=S_true, lat_deg=lat_Seattle, date="2003-10-13"):
    """
    Calculates average daily insolation, given:
    S, total solar irradiance (tsi)
    lat_deg, latitude in degrees
    date, date as a string "yyyy-mm-dd"
    """
    #changing date object to integer to match theta_d equation
    if type(date) is str:
        day = datetime.strptime(date, "%Y-%m-%d").timetuple().tm_yday - 1
    elif type(int(date)) is int:
        day = int(date) - 1
        date = (datetime(2025,1, 1) + timedelta(days=day)).strftime("%m-%d")
    else:
        raise Exception(f"Date is the wrong type! {type(date)} instead of a string")
        
    #pg 355 equation for theta_d
    theta_d = 2*np.pi*day/365
    
    n0 = 0.006918
    n1 = -0.399912*np.cos(theta_d) + 0.070257*np.sin(theta_d)
    n2 = -0.006758*np.cos(2*theta_d) + 0.000907*np.sin(2*theta_d)
    n3 = -0.002697*np.cos(3*theta_d) + 0.00148*np.sin(3*theta_d)
    dec = n0 + n1 + n2 + n3
    
    m0 = 1.00011
    m1 = 0.034221*np.cos(theta_d) + 0.00128*np.sin(theta_d)
    m2 = 0.000719*np.cos(2*theta_d) + 0.000077*np.sin(theta_d)
    d_ratio = m0 + m1 + m2
    
    #convert latitude to radians for formulae
    lat = lat_deg * np.pi/180
    
    #cos(hour angle) & special cases of no sunrise/set
    tans = -np.tan(lat)*np.tan(dec)
    if tans > 1:
        h0 = 0
    elif tans < -1:
        h0 = np.pi
    else:
        h0 = np.arccos(tans)
    
    #eq 2.17 in textbook
    Q = (S/np.pi) * d_ratio * (h0*np.sin(lat)*np.sin(dec) + np.cos(lat)*np.cos(dec)*np.sin(h0))

    return Q

#array of integers of days of yr for slider
days = np.arange(1, 366, dtype=int)

def plot_Q(S=S_true, lat_deg=lat_Seattle, date="2003-10-13"):
    """
    Plotting average daily insolation over the year
    """
    Q = avg_daily_insolation(S, lat_deg, date)
    
    #changing date object to integer to match theta_d equation
    if type(date) is str:
        day = datetime.strptime(date, "%Y-%m-%d").timetuple().tm_yday - 1
    elif type(int(date)) is int:
        day = int(date) - 1
        date = (datetime(2025,1, 1) + timedelta(days=day)).strftime("%d %B")
    else:
        raise Exception(f"Date is the wrong type! {type(date)} instead of a string")
    
    #printing values
    print(f"Average daily insolation: {Q:.2f} W/m^2 at latitude {lat_deg:.0f} degrees on {date}")
    
    #plotting
    plt.figure(figsize=(8,4))
    plt.axvline(day, color='lightseagreen', linestyle='--')
    
    #all avg daily insolation values over a year for given latitude for plot
    Q_array = np.zeros(len(days))
    for i, day in enumerate(days):
        Q_array[i] = avg_daily_insolation(S, lat_deg, date=day.astype(int))

    #plotting
    plt.plot(days, Q_array, color='darkslategrey')
    plt.xlabel("Day of year")
    plt.ylabel("Average daily insolation (W/m$^2$)")
    plt.title(f"Seasonal variation of insolation at latitude {lat_deg:.0f} degrees")
    plt.xlim(0,365)
    plt.ylim(0,600)
    plt.grid()
    plt.show();

def insol_by_lat_day():
    """
    Plots a heatmap of daily average insolation by latitude and day of year
    """
    #daily insolation heatmap
    days = np.arange(1, 366, dtype=int)
    lats = np.linspace(-90, 90, 181)
    Q_lat_day = np.zeros((len(lats), len(days)))
    
    for i, lat in enumerate(lats):
        for j, day in enumerate(days):
            Q_lat_day[i,j] = np.array(avg_daily_insolation(lat_deg=lat, date=day.astype(int)))

    #plotting
    plt.figure(figsize=(8,4))
    plt.imshow(Q_lat_day, extent=[1, 365, -90, 90], aspect='auto', origin='lower', cmap='coolwarm')
    plt.colorbar(label='Insolation (W/m$^2$)')
    plt.xlabel("Day of year")
    plt.ylabel("Latitude (degrees)")
    plt.title("Daily insolation by latitude and day of year")
    plt.show();