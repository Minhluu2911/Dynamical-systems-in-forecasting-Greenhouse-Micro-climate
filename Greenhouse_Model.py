import math
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

meteo_file = "meteo.csv"
climate_file = "Greenhouse_climate.csv"
meteo = pd.read_csv(meteo_file)
data_Tout = np.asarray(meteo["Tout"][2:-1])
wind_speed = np.asarray(meteo["Windsp"][2:-1])
climate = pd.read_csv(climate_file)
data_Tair = np.asarray(climate["Tair"][2:-1])
data_CO2air = np.asarray(climate["CO2air"][2:-1])

GRAVITY = 9.81
ETA_ROOF_THR = 0.9
M_CH2O = 30e-3
EVAPORATION_LATENT_HEAT = 2.45E6
BOUNDARY_LAYER_RESISTANCE = 275
MIN_CANOPY_TRANSPIRATION_RESISTANCE = 82.0
C_PAIR = 1E3
GAMMA = 65.8
ETA_HEATCO2 = 0.057
ETA_HEATVAP = 4.43E-8
M_WATER = 18.01528
M_GAS = 8.314E3
DENSITY_AIR0 = 1.20
M_AIR = 28.96

# Program to calculate Euler algorithm
def euler(dx, t, CO2_Air, CO2_Top, h):
    CO2 = np.array([CO2_Air, CO2_Top])
    d_CO2 = np.array(dx(t, CO2_Air, CO2_Top))
    return CO2 + h * d_CO2


# Program to calculate Runge-Kutta order 4 algorithm
def runge_kutta4th(dx, t, CO2_Air, CO2_Top, h):
    CO2 = np.array([CO2_Air, CO2_Top])

    k1 = np.array(dx(t, CO2_Air, CO2_Top))
    k2 = np.array(dx(t + h / 2, CO2_Air + 0.5 * h * k1[0], CO2_Top + 0.5 * h * k1[1]))
    k3 = np.array(dx(t + h / 2, CO2_Air + 0.5 * h * k2[0], CO2_Top + 0.5 * h * k2[1]))
    k4 = np.array(dx(t + h, CO2_Air + h * k3[0], CO2_Top + h * k3[1]))

    return CO2 + (k1 + 2 * k2 + 2 * k3 + k4) * h / 6




class GreenHouse:
    def __init__(self, initfile):
        filein = pd.read_csv(initfile)

        self.CO2_Air = filein["CO2_Air"][0]
        self.CO2_Top = filein["CO2_Air"][0]
        self.CO2_Out = filein["CO2_Air"][0]

        self.VP_Air = filein["VP_Air"][0]
        self.VP_Top = filein["VP_Air"][0]
        self.VP_Out = filein["VP_Air"][0]

        self.cover_area = filein["cover_area"][0]
        self.floor_area = filein["floor_area"][0]
        self.greenhouse_height = filein["greenhouse_height"][0]
        self.cap_CO2Air = filein["cap_CO2Air"][0]
        self.roof_ventilation_area = filein["roof_ventilation_area"][0]
        self.side_ventilation_area = filein["side_ventilation_area"][0]
        self.elevation_height = filein["elevation_height"][0]
        self.c_HECin = filein["c_HECin"][0]
        self.cap_CO2Top = self.greenhouse_height - self.cap_CO2Air

        self.C_d = filein["C_d"][0]
        self.C_w = filein["C_w"][0]
        self.eta_InsScr = filein["eta_InsScr"][0]
        self.c_leakage = filein["c_leakage"][0]
        self.h_SideRoof = filein["h_SideRoof"][0]
        self.h_Vent = filein["h_Vent"][0]

        self.P_Blow = filein["P_Blow"][0]
        self.phi_ExtCO2 = filein["phi_ExtCO2"][0]
        self.phi_Pad = filein["phi_Pad"][0]
        self.phi_Fog = filein["phi_Fog"][0]
        self.COP_MechCool = filein["COP_MechCool"][0]
        self.P_MechCool = filein["P_MechCool"][0]

        self.K_ThScr = filein["K_ThScr"][0]

        self.eta_CO2_AirStom = filein["eta_CO2_AirStom"][0]

    cover_area: float  # Surface of the cover including side-walls
    floor_area: float  # Surface of the greenhouse floor
    greenhouse_height: float
    cap_CO2Air: float
    cap_CO2Top: float
    roof_ventilation_area: float
    side_ventilation_area: float
    elevation_height: float  # The altitude of the greenhouse
    c_HECin: float  # Convective heat exchange parameter between cover and outdoor air that depends on the greenhouse shape

    C_d: float
    C_w: float
    eta_InsScr: float
    c_leakage: float
    h_SideRoof: float
    h_Vent: float

    P_Blow: float
    phi_ExtCO2: float
    phi_Pad: float
    phi_Fog: float
    COP_MechCool: float
    P_MechCool: float

    K_ThScr: float

    eta_CO2_AirStom: float

    U_Fog = 1
    U_Blow = 1
    U_ExtAir = 0.11
    U_Pad = 1
    U_ThScr = 0.863
    U_Roof = 1
    U_Side = 1
    U_MechCool = 0

    CO2_Air: float
    CO2_Top: float
    T_Air: float
    T_Top: float
    T_ThScr: float
    T_Cov_in: float
    T_MechCool: float
    VP_Air: float
    VP_Top: float

    CO2_Out: float
    T_Out: float
    VP_Out: float
    v_Wind: float

    T_Can: float
    VP_Can: float
    LAI = 2.8

    def update(self, T_air, T_Out, WindSpeed):
        self.T_Air = T_air
        self.T_Top = T_air + 1.5
        self.T_ThScr = T_air + 1
        self.T_Cov_in = T_air + 1
        self.T_MechCool = T_air + 1
        self.T_Can = T_air + 2
        self.VP_Can = GreenHouse.saturation_vapor_pressure(self.T_Can)
        self.T_Out = T_Out
        self.VP_Out = GreenHouse.saturation_vapor_pressure(T_Out)
        self.v_Wind = WindSpeed


    # ================ CO2 fluxes ===================
    def MC_BlowAir(self):  # this Assignment assuming do not have MC_BLowAir
        # Equation  8.54 and 8.53
        # same as equation 3 in pdf
        H_BlowAir = self.H_BlowAir()
        return ETA_HEATCO2 * H_BlowAir

    def MC_ExtAir(self):
        # Equation 8.77
        return self.U_ExtAir * self.phi_ExtCO2 \
               / self.floor_area

    def MC_PadAir(self):
        # do not have PadAir with our model
        f_Pad = self.f_Pad()
        return f_Pad * (self.CO2_Out - self.CO2_Air)

    def MC_AirTop(self):
        f_ThScr = self.f_ThScr()
        return self.air_flux(f_ThScr, self.CO2_Air, self.CO2_Top)

    def MC_AirOut(self):
        f_VentSide = self.f_VentSide()
        f_VentForced = self.f_VentForced()
        f_AirOut = f_VentSide + f_VentForced

        return self.air_flux(f_AirOut, self.CO2_Air, self.CO2_Out)

    def MC_TopOut(self):
        f_VentRoof = self.f_VentRoof()
        return self.air_flux(f_VentRoof, self.CO2_Top, self.CO2_Out)

    def MC_AirCan(self):
        # Equation 9.10
        h_AirCan = 1  # Consider Photosynthesis always happen
        P = self.photosynthesis_rate()
        R = self.photorespiration()
        return M_CH2O * h_AirCan * (P - R)


    # ================ Vapour fluxes ==================
    def MV_BlowAir(self):
        # Equation 8.55
        U_Blow = self.U_Blow
        P_Plow = self.P_Blow
        A_Flr = self.floor_area
        return ETA_HEATVAP * U_Blow * P_Plow / A_Flr

    def MV_AirThScr(self):
        # Equation 8.43 and table 8.4/page239

        # =============== HEC_AirThScr
        U_ThScr = self.U_ThScr
        T_ThScr = self.T_ThScr
        T_Air = self.T_Air
        HEC_AirThScr = 1.7 * U_ThScr * abs(T_Air - T_ThScr) ** 0.33
        # ===========

        VP_Air = self.VP_Air
        VP_ThScr = self.saturation_vapor_pressure(T_ThScr)
        if VP_Air < VP_ThScr:
            return 0
        else:
            return 6.4E-9 * HEC_AirThScr * (VP_Air - VP_ThScr)

    def MV_AirTop(self):
        # Equation 8.45
        f_ThScr = self.f_ThScr()
        VP_Air = self.VP_Air
        T_Air = self.T_Air
        VP_Top = self.VP_Top
        T_Top = self.T_Top
        return self.vapour_flux(f_ThScr, VP_Air, T_Air, VP_Top, T_Top)

    def MV_AirOut(self):
        # Equation 8.45
        f_VentSide = self.f_VentSide()
        f_VentForced = self.f_VentForced()
        f_AirOut = f_VentSide + f_VentForced
        VP_Air = self.VP_Air
        T_Air = self.T_Air
        VP_Out = self.VP_Out
        T_Out = self.T_Out
        return self.vapour_flux(f_AirOut, VP_Air, T_Air, VP_Out, T_Out)

    def MV_AirOutPad(self):
        # Equation 8.59 and 8.62
        U_Pad = self.U_Pad
        phi_Pad = self.phi_Pad
        A_Flr = self.floor_area
        VP_Air = self.VP_Air
        T_Air = self.T_Air
        return U_Pad * phi_Pad / A_Flr * M_WATER / M_GAS * VP_Air / (T_Air + 273.15)

    def MV_AirMech(self):
        #  System does not use
        return 0

    def MV_TopCovin(self):
        # Equation 8.43 and table 8.4/page239

        #  HEC_TopCovin
        c_HECin = self.c_HECin
        T_Top = self.T_Top
        T_Covin = self.T_Cov_in
        A_Cov = self.cover_area
        A_Flr = self.floor_area
        HEC_TopCovin = c_HECin * (T_Top - T_Covin) ** 0.3 * A_Cov / A_Flr
        # ===============
        VP_Top = self.VP_Top
        VP_Covin = self.saturation_vapor_pressure(T_Covin)
        if VP_Top < VP_Covin:
            return 0
        else:
            return 6.4E-9 * HEC_TopCovin * (VP_Top - VP_Covin)

    def MV_CanAir(self):
        # Equation 8.47 and 8.48

        #  VEC Can Air
        p_Air = self.air_density()
        LAI = self.LAI
        VEC_CanAir = 2 * p_Air * C_PAIR * LAI / (
                EVAPORATION_LATENT_HEAT * GAMMA * (
                BOUNDARY_LAYER_RESISTANCE + MIN_CANOPY_TRANSPIRATION_RESISTANCE))
        # =================

        VP_Air = self.VP_Air
        VP_Can = self.saturation_vapor_pressure(self.T_Can)
        return VEC_CanAir * (VP_Can - VP_Air)

    def MV_PadAir(self):
        # Because the system do not use so do not include it
        return 0

    def MV_FogAir(self):
        # Equation 8.64
        U_Fog = self.U_Fog
        phi_Fog = self.phi_Fog
        A_Flr = self.floor_area
        return U_Fog * phi_Fog / A_Flr

    def MV_TopOut(self):
        # Equation 8.45
        f_VentRoof = self.f_VentRoof()
        VP_Top = self.VP_Top
        T_Top = self.T_Top
        VP_Out = self.VP_Out
        T_Out = self.T_Out
        return self.vapour_flux(f_VentRoof, VP_Top, T_Top, VP_Out, T_Out)

    def cap_VPAir(self):
        # Equation 8.25
        h_Air = self.cap_CO2Air
        T_Air = self.T_Air
        return M_WATER * h_Air / (M_GAS * (T_Air + 273.15))

    def cap_VPTop(self):
        # Apply equation 8.25 but with Top
        h_Top = self.cap_CO2Top
        T_Top = self.T_Top
        return M_WATER * h_Top / (M_GAS * (T_Top + 273.15))


    # ================ use to calculate the MC and MV ==================
    def H_BlowAir(self):
        # Equation  8.53
        U_Blow = self.U_Blow
        P_Plow = self.P_Blow
        A_Flr = self.floor_area
        return U_Blow * P_Plow / A_Flr

    def air_flux(self, f, CO2_from, CO2_to):
        # Equation 8.46
        return f * (CO2_from - CO2_to)

    def vapour_flux(self, f, VP_from, T_from, VP_to, T_to):
        # Equation 8.45
        return M_WATER / M_GAS * f * (VP_from / (T_from + 273.15) - VP_to / (T_to + 273.15))

    def f_ThScr(self):
        # Equation 8.41
        U_ThScr = self.U_ThScr
        K_ThScr = self.K_ThScr
        T_Air = self.T_Air
        T_Out = self.T_Out
        p_Air = self.air_density()
        pressure = 101325 * (1 - 2.5577e-5 * self.elevation_height) ** 5.25588
        p_Out = M_AIR * pressure / ((self.T_Top + 273.15) * M_GAS)
        p_Mean = (p_Air + p_Out) / 2

        return U_ThScr * K_ThScr * abs(T_Air - T_Out) ** 0.66 + (1 - U_ThScr) / p_Mean * math.sqrt(
            0.5 * p_Mean * (1 - U_ThScr) * GRAVITY * abs(p_Air - p_Out))

    def f_VentSide(self):
        # Equation 8.73
        eta_Side = 0
        eta_Roof = 1
        eta_InsScr = self.eta_InsScr_()
        d2_f_VentSide = self.d2_f_VentSide()
        d2_f_VentRoofSide = self.d2_f_VentRoofSide()
        f_leakage = self.f_leakage()
        U_ThScr = self.U_ThScr

        if eta_Side >= ETA_ROOF_THR:
            return eta_InsScr * d2_f_VentSide + 0.5 * f_leakage
        else:
            return eta_InsScr * (
                    U_ThScr * d2_f_VentSide + (1 - U_ThScr) * d2_f_VentRoofSide * eta_Side) + 0.5 * f_leakage

    def f_VentForced(self):
        # Do not use
        return 0

    def f_VentRoof(self):
        # Equation 8.72
        eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
        eta_InsScr = self.eta_InsScr_()
        d2_f_VentRoof = self.d2_f_VentRoof()
        d2_f_VentRoofSide = self.d2_f_VentRoofSide()
        f_leakage = self.f_leakage()
        U_ThScr = self.U_ThScr
        if eta_Roof >= ETA_ROOF_THR:
            return eta_InsScr * d2_f_VentRoof + 0.5 * f_leakage
        else:
            return eta_InsScr * (
                    U_ThScr * d2_f_VentRoof + (1 - U_ThScr) * d2_f_VentRoofSide * eta_Roof) + 0.5 * f_leakage

    def f_Pad(self):
        U_Pad = self.U_Pad
        phi_Pad = self.phi_Pad
        A_Flr = self.floor_area
        return U_Pad * phi_Pad / A_Flr

    def air_density(self):
        # Equation 8.24
        return DENSITY_AIR0 * math.exp(
            GRAVITY * M_AIR * self.elevation_height / (
                    293.15 * M_GAS))

    def eta_InsScr_(self):
        # Equation 8.70
        eta_InsScr = self.eta_InsScr
        return eta_InsScr * (1 - eta_InsScr)

    def f_leakage(self):
        # Equation 8.71
        c_leakage = self.c_leakage
        v_Wind = self.v_Wind
        return c_leakage * max(0.25, v_Wind)

    def d2_f_VentRoof(self):
        # Equation 8.65
        U_Roof = self.U_Roof
        A_Roof = self.roof_ventilation_area
        A_Flr = self.floor_area
        C_d = self.C_d
        C_w = self.C_w
        T_Air = self.T_Air
        T_Out = self.T_Out
        T_Mean = (T_Air + T_Out) / 2
        v_Wind = self.v_Wind
        return 0.5 * U_Roof * A_Roof * C_d / A_Flr * math.sqrt(
            0.5 * GRAVITY * (T_Air - T_Out) / (T_Mean + 273.15) + C_w * v_Wind ** 2)

    def d2_f_VentRoofSide(self):
        # Equation 8.66
        U_Roof = self.U_Roof
        U_Side = self.U_Side
        A_Roof = self.roof_ventilation_area
        A_Side = self.side_ventilation_area
        A_Flr = self.floor_area
        h_SideRoof = self.h_SideRoof
        C_d = self.C_d
        C_w = self.C_w
        T_Air = self.T_Air
        T_Out = self.T_Out
        T_Mean = (T_Air + T_Out) / 2
        v_Wind = self.v_Wind
        AU_Roof = A_Roof * U_Roof
        AU_Side = A_Side * U_Side
        return C_d / A_Flr * math.sqrt((AU_Roof * AU_Side / math.sqrt(AU_Roof ** 2 + AU_Side ** 2)) ** 2 * (
                2 * GRAVITY * h_SideRoof * (T_Air - T_Out) / (T_Mean + 273.15)) + 0.25 * (
                                               AU_Roof + AU_Side) ** 2 * C_w * v_Wind ** 2)

    def d2_f_VentSide(self):
        # Equation 8.67
        U_Side = self.U_Side
        A_Side = self.side_ventilation_area
        A_Flr = self.floor_area
        C_d = self.C_d
        C_w = self.C_w
        v_Wind = self.v_Wind
        AU_Side = A_Side * U_Side
        return 0.5 * C_d * AU_Side * v_Wind / A_Flr * math.sqrt(C_w)

    @staticmethod
    def saturation_vapor_pressure(temp):
        return 610.78 * math.exp(temp / (temp + 238.3) * 17.2694)

    def photosynthesis_rate(self):
        J = self.electron_transport_rate()
        CO2_Stom = self.CO2_Stom()
        gamma = self.gamma()
        return 0.25 * J * (CO2_Stom - gamma) / (CO2_Stom + 2 * gamma)

    def photorespiration(self):
        P = self.photosynthesis_rate()
        gamma = self.gamma()
        CO2_Stom = self.CO2_Stom()
        return P * gamma / CO2_Stom

    def electron_transport_rate(self):
        # Equation 9.14
        theta = 0.7
        alpha = 0.385
        # ================ JPot (9.15) ==============
        R = 8.314
        S = 710
        T_CanK = self.T_Can + 273.15
        T_25K = 298.15
        Ej = 37E3
        H = 22E4
        LAI = self.LAI
        J_max_25Leaf = 210
        J_max_25Can = J_max_25Leaf * LAI  # Equation 9.16
        J_POT = J_max_25Can * math.exp(Ej * (T_CanK - T_25K) / (R * T_CanK * T_25K)) * (
                1 + math.exp((S * T_25K - H) / (R * T_25K))) / (1 + math.exp((S * T_CanK - H) / (R * T_CanK)))
        # ================ PAR_Can (9.17) ==============
        # 9.19
        K1 = 0.7
        K2 = 0.7
        p_Can = 0.07
        p_Flr = 0.5
        light_transmission = 0.78
        eta_Glob_PAR = 2.3
        I_Glob = 225
        PAR_Gh = light_transmission * eta_Glob_PAR * I_Glob
        PAR_GhCan = PAR_Gh * (1 - p_Can) * (1 - math.exp(-K1 * LAI))
        PAR_FlrCan = p_Flr * PAR_Gh * (1 - p_Can) * math.exp(-K1 * LAI) * (1 - math.exp(-K2 * LAI))
        PAR_Can = PAR_FlrCan + PAR_GhCan

        return 0.5 * (J_POT + alpha * PAR_Can - math.sqrt(
            (J_POT + alpha * PAR_Can) ** 2 - 4 * theta * J_POT * alpha * PAR_Can)) / theta

    def gamma(self):
        # Equation 9.23
        LAI = self.LAI
        J_max_25Leaf = 210
        J_max_25Can = J_max_25Leaf * LAI  # Equation 9.16
        c_gamma = 1.7
        T_Can = self.T_Can

        return J_max_25Leaf / J_max_25Can * c_gamma * T_Can + 20 * c_gamma * (1 - J_max_25Leaf / J_max_25Can)

    def CO2_Stom(self):
        eta_CO2_AirStom = self.eta_CO2_AirStom
        return eta_CO2_AirStom * self.CO2_Air


    # ================ ODE ====================
    def dx_CO2(self, t, CO2_Air, CO2_Top):
        # Update environment variable at time t (s)
        self.update(data_Tair[int(t / 300)], data_Tout[int(t / 300)], wind_speed[int(t / 300)])

        self.CO2_Air = CO2_Air
        self.CO2_Top = CO2_Top

        return (
                       self.MC_BlowAir() + self.MC_ExtAir() + self.MC_PadAir() - self.MC_AirCan() - self.MC_AirTop() - self.MC_AirOut()) / self.cap_CO2Air, (
                       self.MC_AirTop() - self.MC_TopOut()) / self.cap_CO2Top

    def dx_VP(self, t, VP_Air, VP_Top):
        # Update environment variable at time t (s)
        self.update(data_Tair[int(t / 300)], data_Tout[int(t / 300)], wind_speed[int(t / 300)])
        self.VP_Air = VP_Air
        self.VP_Top = VP_Top

        return (self.MV_CanAir() + self.MV_PadAir() + self.MV_FogAir() + self.MV_BlowAir() -
                self.MV_AirThScr() - self.MV_AirTop() - self.MV_AirOut() -
                self.MV_AirOutPad() - self.MV_AirMech()) / self.cap_VPAir(), \
               (self.MV_AirTop() - self.MV_TopCovin() - self.MV_TopOut()) / self.cap_VPTop()


if __name__ == '__main__':
    model = GreenHouse("init.csv")
    data_VPAir = [GreenHouse.saturation_vapor_pressure(t) for t in data_Tair]

    predict_CO2air = []
    predict_CO2top = []
    predict_VPair = []
    predict_VPtop = []

    h = 5  # step(s)
    num_day = 7
    for t in range(0, num_day * 86400, h):
        CO2 = euler(model.dx_CO2, t, model.CO2_Air, model.CO2_Top, h)
        model.CO2_Air, model.CO2_Top = CO2[0], CO2[1]
        VP = euler(model.dx_VP, t, model.VP_Air, model.VP_Top, h)
        model.VP_Air, model.VP_Top = VP[0], VP[1]
        if t % 300 == 0:
            predict_CO2air.append(model.CO2_Air)
            predict_CO2top.append(model.CO2_Top)
            predict_VPair.append(model.VP_Air)
            predict_VPtop.append(model.VP_Top)

    # Calculate the error
    N = len(predict_CO2air)
    RMSE_CO2 = 0
    for co2_pre, co2 in zip(predict_CO2air, data_CO2air[:N]):
        RMSE_CO2 += (co2 - co2_pre) ** 2
    RMSE_CO2 = (RMSE_CO2 / N) ** (1 / 2)
    print("Root Mean Squared Error of CO2: ", RMSE_CO2)

    RMSE_VP = 0
    for VP_pre, VP in zip(predict_VPair, data_CO2air[:N]):
        RMSE_VP += (VP - VP_pre) ** 2
    RMSE_VP = (RMSE_VP / N) ** (1 / 2)
    print("Root Mean Squared Error of Vapor Pressure: ", RMSE_VP)

    # plot CO2
    plt.subplot(121)
    plt.title("CO2")
    plt.xlabel("Step : 5 min")
    plt.ylabel("Concentration of C02 (mg m^-3)")
    plt.plot(data_CO2air[:N], color="#4aff02")
    plt.plot(predict_CO2air, color="#ff024a", linewidth=2)
    plt.legend(["Reality", "Predict"])

    # plot VP
    plt.subplot(122)
    plt.title("Vapor Pressure")
    plt.ylabel("Vapor Pressure (Pa)")
    plt.xlabel("Step : 5 min")
    plt.plot(data_VPAir[:N], color="#4aff02")
    plt.plot(predict_VPair, color="#ff024a", linewidth=2)
    plt.legend(["Reality", "Predict"])

    plt.show()
