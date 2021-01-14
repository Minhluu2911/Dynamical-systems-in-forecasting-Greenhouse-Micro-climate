import math
import pandas as pd
import read_data

# ===============================================================
# Environment data
# ===============================================================
data_Tout, wind_speed, data_Tair = read_data.read_data("meteo.csv", "Greenhouse_climate.csv")


# ===============================================================
# Greenhouse Model
# ===============================================================
class GreenHouse:
    def __init__(self, initfile):
        filein = pd.read_csv(initfile)

        # at t0 Air = Top = Out
        self.ClimateStates.CO2_Air = filein["CO2_Air"][0]
        self.ClimateStates.CO2_Top = filein["CO2_Air"][0]
        self.Weather.CO2_Out = filein["CO2_Air"][0]

        self.ClimateStates.VP_Air = filein["VP_Air"][0]
        self.ClimateStates.VP_Top = filein["VP_Air"][0]
        self.Weather.VP_Out = filein["VP_Air"][0]

        self.Coefficients.Construction.cover_area = filein["cover_area"][0]
        self.Coefficients.Construction.floor_area = filein["floor_area"][0]
        self.Coefficients.Construction.greenhouse_height = filein["greenhouse_height"][0]
        self.Coefficients.Construction.cap_CO2Air = filein["cap_CO2Air"][0]
        self.Coefficients.Construction.roof_ventilation_area = filein["roof_ventilation_area"][0]
        self.Coefficients.Construction.side_ventilation_area = filein["side_ventilation_area"][0]
        self.Coefficients.Construction.elevation_height = filein["elevation_height"][0]
        self.Coefficients.Construction.c_HECin = filein["c_HECin"][0]
        self.Coefficients.Construction.cap_CO2Top = self.Coefficients.Construction.greenhouse_height - self.Coefficients.Construction.cap_CO2Air

        self.Coefficients.Ventilation.C_d = filein["C_d"][0]
        self.Coefficients.Ventilation.C_w = filein["C_w"][0]
        self.Coefficients.Ventilation.eta_InsScr = filein["eta_InsScr"][0]
        self.Coefficients.Ventilation.c_leakage = filein["c_leakage"][0]
        self.Coefficients.Ventilation.h_SideRoof = filein["h_SideRoof"][0]
        self.Coefficients.Ventilation.h_Vent = filein["h_Vent"][0]

        self.Coefficients.ActiveClimateControl.P_Blow = filein["P_Blow"][0]
        self.Coefficients.ActiveClimateControl.phi_ExtCO2 = filein["phi_ExtCO2"][0]
        self.Coefficients.ActiveClimateControl.phi_Pad = filein["phi_Pad"][0]
        self.Coefficients.ActiveClimateControl.phi_Fog = filein["phi_Fog"][0]
        self.Coefficients.ActiveClimateControl.COP_MechCool = filein["COP_MechCool"][0]
        self.Coefficients.ActiveClimateControl.P_MechCool = filein["P_MechCool"][0]

        self.Coefficients.Thermalscreen.K_ThScr = filein["K_ThScr"][0]

        self.Coefficients.Photosynthesis.eta_CO2_AirStom = filein["eta_CO2_AirStom"][0]

    class Constant:
        # The acceleration of gravity. Unit: m s^-2
        GRAVITY = 9.81
        # The ratio between the roof vent area and total ventilation area above no chimney effect was assumed.
        # Unit: [-].
        # Ref: Assumed by Vanthoor
        ETA_ROOF_THR = 0.9

        M_CH2O = 30e-3

        # Latent heat of evaporation.
        # Unit: J kg^-1 {water}
        EVAPORATION_LATENT_HEAT = 2.45E6

        # Boundary layer resistance of the canopy for vapor transport.
        # Unit: s m^-1.
        # Ref: Mean value based on Stanghellini (Stanghellini, 1987)
        BOUNDARY_LAYER_RESISTANCE = 275

        # The minimum canopy resistance for transpiration.
        # Unit: s m^-1.
        # Ref: Stanghellini (1987)
        MIN_CANOPY_TRANSPIRATION_RESISTANCE = 82.0

        # Specific heat capacity of the air.
        # Unit: J K^-1 kg^-1
        C_PAIR = 1E3

        # Psychrometric constant.
        # Unit: Pa K^-1
        GAMMA = 65.8

        # Amount of CO2 which is released when 1 Joule of sensible energy is produced by the heat blower.
        # Unit: mg {CO2} J^-1
        # MC_BlowAir
        ETA_HEATCO2 = 0.057

        # Amount of vapor which is released when 1 Joule of sensible energy is produced by the heat blower.
        # Unit: kg {vapour} J-1
        ETA_HEATVAP = 4.43E-8

        # Molar mass of water
        # Unit: kg kmol^-1
        M_WATER = 18.01528

        # Molar gas constant.
        # Unit: J Kmol^-1 K^-1
        M_GAS = 8.314E3

        # Density of the air at sea level.
        # Unit: kg m^-3
        DENSITY_AIR0 = 1.20

        # Molar mass of air.
        # Unit: kg kmol^-1
        M_AIR = 28.96

    class Coefficients:
        class Construction:
            cover_area: float  # Surface of the cover including side-walls
            floor_area: float  # Surface of the greenhouse floor
            greenhouse_height: float
            cap_CO2Air: float
            cap_CO2Top: float
            roof_ventilation_area: float
            side_ventilation_area: float
            elevation_height: float  # The altitude of the greenhouse
            c_HECin: float  # Convective heat exchange parameter between cover and outdoor air that depends on the greenhouse shape

        class Ventilation:
            C_d: float
            C_w: float
            eta_InsScr: float
            c_leakage: float
            h_SideRoof: float
            h_Vent: float

        class ActiveClimateControl:
            P_Blow: float
            phi_ExtCO2: float
            phi_Pad: float
            phi_Fog: float
            COP_MechCool: float
            P_MechCool: float

        class Thermalscreen:
            K_ThScr: float

        class Photosynthesis:
            eta_CO2_AirStom: float

    class Assuming:
        U_Fog = 1
        U_Blow = 1
        U_ExtAir = 0.11
        U_Pad = 1
        U_ThScr = 0.863
        U_Roof = 1
        U_Side = 1
        U_MechCool = 0

    class ClimateStates:
        CO2_Air: float
        CO2_Top: float
        T_Air: float
        T_Top: float
        T_ThScr: float
        T_Cov_in: float
        T_MechCool: float
        VP_Air: float
        VP_Top: float

    class Weather:
        CO2_Out: float
        T_Out: float
        VP_Out: float
        v_Wind: float

    class Tomato:
        T_Can: float
        VP_Can: float
        LAI = 2.8

    # =====================================
    # ================ Property
    # =====================================
    @property
    def CO2_Air(self):
        return self.ClimateStates.CO2_Air

    @CO2_Air.setter
    def CO2_Air(self, CO2_Air):
        self.ClimateStates.CO2_Air = CO2_Air

    @property
    def CO2_Top(self):
        return self.ClimateStates.CO2_Top

    @CO2_Top.setter
    def CO2_Top(self, CO2_Top):
        self.ClimateStates.CO2_Top = CO2_Top

    @property
    def VP_Air(self):
        return self.ClimateStates.VP_Air

    @VP_Air.setter
    def VP_Air(self, VP_Air):
        self.ClimateStates.VP_Air = VP_Air

    @property
    def VP_Top(self):
        return self.ClimateStates.VP_Top

    @VP_Top.setter
    def VP_Top(self, VP_Top):
        self.ClimateStates.VP_Top = VP_Top

    def update(self, T_air, T_Out, WindSpeed):
        # Assuming difference temperature is constant
        self.ClimateStates.T_Air = T_air
        self.ClimateStates.T_Top = T_air + 1.5
        self.ClimateStates.T_ThScr = T_air + 1
        self.ClimateStates.T_Cov_in = T_air + 1
        self.ClimateStates.T_MechCool = T_air + 1
        self.Tomato.T_Can = T_air + 2

        self.Tomato.VP_Can = GreenHouse.saturation_vapor_pressure(self.Tomato.T_Can)

        self.Weather.T_Out = T_Out
        self.Weather.VP_Out = GreenHouse.saturation_vapor_pressure(T_Out)

        self.Weather.v_Wind = WindSpeed

    # =====================================
    # ================ CO2 fluxes
    # =====================================
    def MC_BlowAir(self):  # this Assignment assuming do not have MC_BLowAir
        # Equation  8.54 and 8.53
        # same as equation 3 in pdf
        H_BlowAir = self.H_BlowAir()
        return self.Constant.ETA_HEATCO2 * H_BlowAir

    def MC_ExtAir(self):
        # Equation 8.77
        U_ExtCO2 = self.Assuming.U_ExtAir
        phi_ExtCO2 = self.Coefficients.ActiveClimateControl.phi_ExtCO2
        A_Flr = self.Coefficients.Construction.floor_area

        return U_ExtCO2 * phi_ExtCO2 / A_Flr

    def MC_PadAir(self):
        # Mo hinh dang xet khong co he thong quat gio
        f_Pad = self.f_Pad()
        return f_Pad * (self.Weather.CO2_Out - self.CO2_Air)

    def MC_AirTop(self):
        f_ThScr = self.f_ThScr(1)
        return self.air_flux(f_ThScr, self.CO2_Air, self.CO2_Top)

    def MC_AirOut(self):
        f_VentSide = self.f_VentSide()
        f_VentForced = self.f_VentForced()
        f_AirOut = f_VentSide + f_VentForced

        return self.air_flux(f_AirOut, self.CO2_Air, self.Weather.CO2_Out)

    def MC_TopOut(self):
        f_VentRoof = self.f_VentRoof()
        return self.air_flux(f_VentRoof, self.CO2_Top, self.Weather.CO2_Out)

    def MC_AirCan(self):
        # Equation 9.10
        h_AirCan = 1  # Consider Photosynthesis always happen
        P = self.photosynthesis_rate()
        R = self.photorespiration()
        return self.Constant.M_CH2O * h_AirCan * (P - R)

    # =====================================
    # ================ Vapour fluxes
    # =====================================

    def MV_BlowAir(self):
        # Equation 8.55
        U_Blow = self.Assuming.U_Blow
        P_Plow = self.Coefficients.ActiveClimateControl.P_Blow
        A_Flr = self.Coefficients.Construction.floor_area
        return self.Constant.ETA_HEATVAP * U_Blow * P_Plow / A_Flr

    def MV_AirThScr(self):
        # Equation 8.43 and table 8.4/page239

        # =============== HEC_AirThScr
        U_ThScr = self.Assuming.U_ThScr
        T_ThScr = self.ClimateStates.T_ThScr
        T_Air = self.ClimateStates.T_Air
        HEC_AirThScr = 1.7 * U_ThScr * abs(T_Air - T_ThScr) ** 0.33
        # ===========

        VP_Air = self.ClimateStates.VP_Air
        VP_ThScr = self.saturation_vapor_pressure(T_ThScr)
        if VP_Air < VP_ThScr:
            return 0
        else:
            return 6.4E-9 * HEC_AirThScr * (VP_Air - VP_ThScr)

    def MV_AirTop(self):
        # Equation 8.45
        f_ThScr = self.f_ThScr(0)
        VP_Air = self.ClimateStates.VP_Air
        T_Air = self.ClimateStates.T_Air
        VP_Top = self.ClimateStates.VP_Top
        T_Top = self.ClimateStates.T_Top
        return self.vapour_flux(f_ThScr, VP_Air, T_Air, VP_Top, T_Top)

    def MV_AirOut(self):
        # Equation 8.45
        f_VentSide = self.f_VentSide()
        f_VentForced = self.f_VentForced()
        f_AirOut = f_VentSide + f_VentForced
        VP_Air = self.ClimateStates.VP_Air
        T_Air = self.ClimateStates.T_Air
        VP_Out = self.Weather.VP_Out
        T_Out = self.Weather.T_Out
        return self.vapour_flux(f_AirOut, VP_Air, T_Air, VP_Out, T_Out)

    def MV_AirOutPad(self):
        # Equation 8.59 and 8.62
        U_Pad = self.Assuming.U_Pad
        phi_Pad = self.Coefficients.ActiveClimateControl.phi_Pad
        A_Flr = self.Coefficients.Construction.floor_area
        VP_Air = self.ClimateStates.VP_Air
        T_Air = self.ClimateStates.T_Air
        return U_Pad * phi_Pad / A_Flr * self.Constant.M_WATER / self.Constant.M_GAS * VP_Air / (T_Air + 273.15)

    def MV_AirMech(self):
        # ============== System does not use
        return 0

    def MV_TopCovin(self):
        # Equation 8.43 and table 8.4/page239

        # =============== HEC_TopCovin
        c_HECin = self.Coefficients.Construction.c_HECin
        T_Top = self.ClimateStates.T_Top
        T_Covin = self.ClimateStates.T_Cov_in
        A_Cov = self.Coefficients.Construction.cover_area
        A_Flr = self.Coefficients.Construction.floor_area
        HEC_TopCovin = c_HECin * (T_Top - T_Covin) ** 0.3 * A_Cov / A_Flr
        # ===============
        VP_Top = self.ClimateStates.VP_Top
        VP_Covin = self.saturation_vapor_pressure(T_Covin)
        if VP_Top < VP_Covin:
            return 0
        else:
            return 6.4E-9 * HEC_TopCovin * (VP_Top - VP_Covin)

    def MV_CanAir(self):
        # Equation 8.47 and 8.48

        # ================= VEC Can Air
        p_Air = self.air_density()
        LAI = self.Tomato.LAI
        VEC_CanAir = 2 * p_Air * self.Constant.C_PAIR * LAI / (
                self.Constant.EVAPORATION_LATENT_HEAT * self.Constant.GAMMA * (
                self.Constant.BOUNDARY_LAYER_RESISTANCE + self.Constant.MIN_CANOPY_TRANSPIRATION_RESISTANCE))
        # =================

        VP_Air = self.ClimateStates.VP_Air
        VP_Can = self.saturation_vapor_pressure(self.Tomato.T_Can)
        return VEC_CanAir * (VP_Can - VP_Air)

    def MV_PadAir(self):
        # Because the system do not use so do not include it
        return 0

    def MV_FogAir(self):
        # Equation 8.64
        U_Fog = self.Assuming.U_Fog
        phi_Fog = self.Coefficients.ActiveClimateControl.phi_Fog
        A_Flr = self.Coefficients.Construction.floor_area
        return U_Fog * phi_Fog / A_Flr

    def MV_TopOut(self):
        # Equation 8.45
        f_VentRoof = self.f_VentRoof()
        VP_Top = self.ClimateStates.VP_Top
        T_Top = self.ClimateStates.T_Top
        VP_Out = self.Weather.VP_Out
        T_Out = self.Weather.T_Out
        return self.vapour_flux(f_VentRoof, VP_Top, T_Top, VP_Out, T_Out)

    def cap_VPAir(self):
        # Equation 8.25
        h_Air = self.Coefficients.Construction.cap_CO2Air
        T_Air = self.ClimateStates.T_Air
        return self.Constant.M_WATER * h_Air / (self.Constant.M_GAS * (T_Air + 273.15))

    def cap_VPTop(self):
        # Apply equation 8.25 but with Top
        h_Top = self.Coefficients.Construction.cap_CO2Top
        T_Top = self.ClimateStates.T_Top
        return self.Constant.M_WATER * h_Top / (self.Constant.M_GAS * (T_Top + 273.15))

    # ==================================================================
    # ================ use to calculate the MC and MV ==================
    # ==================================================================
    def H_BlowAir(self):
        # Equation  8.53
        U_Blow = self.Assuming.U_Blow
        P_Plow = self.Coefficients.ActiveClimateControl.P_Blow  # zero
        A_Flr = self.Coefficients.Construction.floor_area
        return U_Blow * P_Plow / A_Flr

    def air_flux(self, f, CO2_from, CO2_to):
        # Equation 8.46
        return f * (CO2_from - CO2_to)

    def vapour_flux(self, f, VP_from, T_from, VP_to, T_to):
        # Equation 8.45
        return self.Constant.M_WATER / self.Constant.M_GAS * f * (VP_from / (T_from + 273.15) - VP_to / (T_to + 273.15))

    def f_ThScr(self,input_):
        # Equation 8.41
        U_ThScr = self.Assuming.U_ThScr
        K_ThScr = self.Coefficients.Thermalscreen.K_ThScr
        T_Air = self.ClimateStates.T_Air
        T_Out = self.Weather.T_Out
        p_Air = self.air_density()
        pressure = 101325 * (1 - 2.5577e-5 * self.Coefficients.Construction.elevation_height) ** 5.25588
        p_Out = self.Constant.M_AIR * pressure / ((self.ClimateStates.T_Top + 273.15) * self.Constant.M_GAS)
        p_Mean = (p_Air + p_Out) / 2

        return (U_ThScr * K_ThScr * abs(T_Air - T_Out) ** (2 / 3)) + (1 - U_ThScr) * math.sqrt(
                (0.5 * (1 - U_ThScr) * self.Constant.GRAVITY * abs(p_Air - p_Out)) / p_Mean)


    def f_VentSide(self):
        # Equation 8.73
        eta_Side = 0  # Note: line 611 / setGlAux / GreenLight
        eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
        eta_InsScr = self.eta_InsScr()
        d2_f_VentSide = self.d2_f_VentSide()
        d2_f_VentRoofSide = self.d2_f_VentRoofSide()
        f_leakage = self.f_leakage()
        U_ThScr = self.Assuming.U_ThScr

        if eta_Side >= self.Constant.ETA_ROOF_THR:
            return eta_InsScr * d2_f_VentSide + 0.5 * f_leakage
        else:
            return eta_InsScr * (
                    U_ThScr * d2_f_VentSide + (1 - U_ThScr) * d2_f_VentRoofSide * eta_Side) + 0.5 * f_leakage

    def f_VentForced(self):
        return 0

    def f_VentRoof(self):
        # Equation 8.72
        eta_Roof = 1  # Note: line 606 / setGlAux / GreenLight
        eta_InsScr = self.eta_InsScr()
        d2_f_VentRoof = self.d2_f_VentRoof()
        d2_f_VentRoofSide = self.d2_f_VentRoofSide()
        f_leakage = self.f_leakage()
        U_ThScr = self.Assuming.U_ThScr
        if eta_Roof >= self.Constant.ETA_ROOF_THR:
            return eta_InsScr * d2_f_VentRoof + 0.5 * f_leakage
        else:
            return eta_InsScr * (
                    U_ThScr * d2_f_VentRoof + (1 - U_ThScr) * d2_f_VentRoofSide * eta_Roof) + 0.5 * f_leakage

    def f_Pad(self):
        U_Pad = self.Assuming.U_Pad
        phi_Pad = self.Coefficients.ActiveClimateControl.phi_Pad
        A_Flr = self.Coefficients.Construction.floor_area
        return U_Pad * phi_Pad / A_Flr

    def air_density(self):
        # Equation 8.24
        return self.Constant.DENSITY_AIR0 * math.exp(
            self.Constant.GRAVITY * self.Constant.M_AIR * self.Coefficients.Construction.elevation_height / (
                    293.15 * self.Constant.M_GAS))

    def eta_InsScr(self):
        # Equation 8.70
        eta_InsScr = self.Coefficients.Ventilation.eta_InsScr
        return eta_InsScr * (2 - eta_InsScr)

    def f_leakage(self):
        # Equation 8.71
        c_leakage = self.Coefficients.Ventilation.c_leakage
        v_Wind = self.Weather.v_Wind
        return c_leakage * max(0.25, v_Wind)

    def d2_f_VentRoof(self):
        # Equation 8.65
        U_Roof = self.Assuming.U_Roof
        A_Roof = self.Coefficients.Construction.roof_ventilation_area
        A_Flr = self.Coefficients.Construction.floor_area
        C_d = self.Coefficients.Ventilation.C_d
        C_w = self.Coefficients.Ventilation.C_w
        T_Air = self.ClimateStates.T_Air
        T_Out = self.Weather.T_Out
        T_Mean = (T_Air + T_Out) / 2
        v_Wind = self.Weather.v_Wind
        return 0.5 * U_Roof * A_Roof * C_d / A_Flr * math.sqrt(
            0.5 * self.Constant.GRAVITY * (T_Air - T_Out) / (T_Mean + 273.15) + C_w * v_Wind ** 2)

    def d2_f_VentRoofSide(self):
        # Equation 8.66
        U_Roof = self.Assuming.U_Roof
        U_Side = self.Assuming.U_Side
        A_Roof = self.Coefficients.Construction.roof_ventilation_area
        A_Side = self.Coefficients.Construction.side_ventilation_area
        A_Flr = self.Coefficients.Construction.floor_area
        h_SideRoof = self.Coefficients.Ventilation.h_SideRoof
        C_d = self.Coefficients.Ventilation.C_d
        C_w = self.Coefficients.Ventilation.C_w
        T_Air = self.ClimateStates.T_Air
        T_Out = self.Weather.T_Out
        T_Mean = (T_Air + T_Out) / 2
        v_Wind = self.Weather.v_Wind
        AU_Roof = A_Roof * U_Roof
        AU_Side = A_Side * U_Side
        return C_d / A_Flr * math.sqrt((AU_Roof * AU_Side / math.sqrt(AU_Roof ** 2 + AU_Side ** 2)) ** 2 * (
                2 * self.Constant.GRAVITY * h_SideRoof * (T_Air - T_Out) / (T_Mean + 273.15)) + 0.25 * (
                                               AU_Roof + AU_Side) ** 2 * C_w * v_Wind ** 2)

    def d2_f_VentSide(self):
        # Equation 8.67
        U_Side = self.Assuming.U_Side
        A_Side = self.Coefficients.Construction.side_ventilation_area
        A_Flr = self.Coefficients.Construction.floor_area
        C_d = self.Coefficients.Ventilation.C_d
        C_w = self.Coefficients.Ventilation.C_w
        v_Wind = self.Weather.v_Wind
        AU_Side = A_Side * U_Side
        return 0.5 * C_d * AU_Side * v_Wind / A_Flr * math.sqrt(C_w)

    @staticmethod
    def saturation_vapor_pressure(temp):
        # Calculation based on
        # http://www.conservationphysics.org/atmcalc/atmoclc2.pdf
        return 610.78 * math.exp(temp / (temp + 238.3) * 17.2694)  # Pascal

    # =====================================
    # ================ Photosynthesis
    # =====================================
    def photosynthesis_rate(self):
        # Equation 9.12
        J = self.electron_transport_rate()
        CO2_Stom = self.CO2_Stom()
        gamma = self.gamma()
        return 0.25 * J * (CO2_Stom - gamma) / (CO2_Stom + 2 * gamma)

    def photorespiration(self):
        # Equation 9.13
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
        T_CanK = self.Tomato.T_Can + 273.15
        T_25K = 298.15
        Ej = 37E3
        H = 22E4
        LAI = self.Tomato.LAI
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
        LAI = self.Tomato.LAI
        J_max_25Leaf = 210
        J_max_25Can = J_max_25Leaf * LAI  # Equation 9.16
        c_gamma = 1.7
        T_Can = self.Tomato.T_Can

        return J_max_25Leaf / J_max_25Can * c_gamma * T_Can + 20 * c_gamma * (1 - J_max_25Leaf / J_max_25Can)

    def CO2_Stom(self):
        eta_CO2_AirStom = self.Coefficients.Photosynthesis.eta_CO2_AirStom
        return eta_CO2_AirStom * self.CO2_Air

    # =====================================
    # ================ ODE
    # =====================================
    def dx_CO2(self, t, CO2_Air, CO2_Top):
        # Update environment variable at time t (s)
        self.update(data_Tair[int(t / 300)], data_Tout[int(t / 300)], wind_speed[int(t / 300)])

        self.CO2_Air = CO2_Air
        self.CO2_Top = CO2_Top

        cap_CO2Air = self.Coefficients.Construction.cap_CO2Air
        cap_CO2Top = self.Coefficients.Construction.cap_CO2Top
        MC_BlowAir = self.MC_BlowAir()
        MC_ExtAir = self.MC_ExtAir()
        MC_PadAir = self.MC_PadAir()
        MC_AirTop = self.MC_AirTop()
        MC_AirOut = self.MC_AirOut()
        MC_TopOut = self.MC_TopOut()
        MC_AirCan = self.MC_AirCan()
        return (MC_BlowAir + MC_ExtAir + MC_PadAir - MC_AirCan - MC_AirTop - MC_AirOut) / cap_CO2Air, (
                MC_AirTop - MC_TopOut) / cap_CO2Top

    def dx_VP(self, t, VP_Air, VP_Top):
        # Update environment variable at time t (s)
        self.update(data_Tair[int(t / 300)], data_Tout[int(t / 300)], wind_speed[int(t / 300)])

        self.VP_Air = VP_Air
        self.VP_Top = VP_Top

        cap_VPAir = self.cap_VPAir()
        cap_VPTop = self.cap_VPTop()
        MV_CanAir = self.MV_CanAir()
        MV_PadAir = self.MV_PadAir()
        MV_FogAir = self.MV_FogAir()
        MV_BlowAir = self.MV_BlowAir()
        MV_AirThScr = self.MV_AirThScr()
        MV_AirTop = self.MV_AirTop()
        MV_AirOut = self.MV_AirOut()
        MV_AirOutPad = self.MV_AirOutPad()
        MV_AirMech = self.MV_AirMech()
        MV_TopCovin = self.MV_TopCovin()
        MV_TopOut = self.MV_TopOut()
        return (MV_CanAir + MV_PadAir + MV_FogAir + MV_BlowAir - MV_AirThScr - MV_AirTop -
                MV_AirOut - MV_AirOutPad - MV_AirMech) / cap_VPAir, (
                       MV_AirTop - MV_TopCovin - MV_TopOut) / cap_VPTop
