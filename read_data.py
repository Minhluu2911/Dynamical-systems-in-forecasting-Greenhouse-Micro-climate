import pandas as pd
import numpy as np


def read_data(meteo_file, climate_file):
    meteo = pd.read_csv(meteo_file)
    T_Out = meteo["Tout"][2:-1]
    wind_spped = meteo["Windsp"][2:-1]

    climate = pd.read_csv(climate_file)
    T_Air = climate["Tair"][2:-1]
    return np.asarray(T_Out), np.asarray(wind_spped), np.asarray(T_Air)