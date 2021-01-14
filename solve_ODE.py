import numpy as np


def runge_kutta4th(dx, t, CO2_Air, CO2_Top, h):
    CO2 = np.array([CO2_Air, CO2_Top])

    k1 = np.array(dx(t, CO2_Air, CO2_Top))
    k2 = np.array(dx(t + h / 2, CO2_Air + 0.5 * h * k1[0], CO2_Top + 0.5 * h * k1[1]))
    k3 = np.array(dx(t + h / 2, CO2_Air + 0.5 * h * k2[0], CO2_Top + 0.5 * h * k2[1]))
    k4 = np.array(dx(t + h, CO2_Air + h * k3[0], CO2_Top + h * k3[1]))

    return CO2 + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6


def euler(dx, t, CO2_Air, CO2_Top, h):
    CO2 = np.array([CO2_Air, CO2_Top])

    d_CO2 = np.array(dx(t, CO2_Air, CO2_Top))

    return CO2 + h * d_CO2
