from Greenhouse_Model import GreenHouse
from solve_ODE import runge_kutta4th, euler
import matplotlib.pyplot as plt
import pandas as pd
import csv

# Permission to write on file
Permission_Write_Predict_File = False

# ============================================================
# Data from the real life
T_air = pd.read_csv("Greenhouse_climate.csv")["Tair"][2:-1]
data_CO2Air = pd.read_csv("Greenhouse_climate.csv")["CO2air"][2:-1]
data_VPAir = [GreenHouse.saturation_vapor_pressure(t) for t in T_air]
# =============================================================
initFileName = "init.csv"
model = GreenHouse(initFileName)


# ============================================================
# Simulate

predict_CO2air = []
predict_CO2top= []
predict_VPair = []
predict_VPtop = []

h = 10 # step(s)
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
# ============================================================


# ============================================================
# Calculate the error
N = len(predict_CO2air)
RMSE_CO2 = 0
for co2_pre, co2 in zip(predict_CO2air, data_CO2Air[:N]):
    RMSE_CO2 += (co2 - co2_pre)**2
RMSE_CO2 = (RMSE_CO2 / N)**(1/2)
print("Root Mean Squared Error of CO2: ", RMSE_CO2)

RMSE_VP = 0
for VP_pre, VP in zip(predict_VPair, data_CO2Air[:N]):
    RMSE_VP += (VP -VP_pre)**2
RMSE_VP = (RMSE_VP / N)**(1/2)
print("Root Mean Squared Error of Vapor Pressure: ", RMSE_VP)


# ============================================================
# ============================================================

# visual CO2==================
plt.subplot(121)
plt.title("CO2")
plt.xlabel("Step : 5 min")
plt.ylabel("Concentration of C02 (mg m^-3)")
plt.plot(data_CO2Air[:N], color="#01FB00")
plt.plot(predict_CO2air, color="#FB0000", linewidth=2)
plt.legend(["Reality","Predict"])

# visual VP==================
plt.subplot(122)
plt.title("Vapor Pressure")
plt.ylabel("Vapor Pressure (Pa)")
plt.xlabel("Step : 5 min")
plt.plot(data_VPAir[:N], color="#01FB00")
plt.plot(predict_VPair, color="#FB0000", linewidth=2)
plt.legend(["Reality","Predict"])

plt.show()


# ============================================================
# ========================================
# tao file de ghi ket qua
if Permission_Write_Predict_File:
    writer = csv.writer(open("predict.csv", mode='w'), delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
    writer.writerow(['t(h)', 'CO2_air(mg m^-3)', 'CO2_top(mg m^-3)', 'VP_air(Pa)', 'Vp_top(Pa)'])
    for i in range(0, N, 12):
        writer.writerow([i//12,predict_CO2air[i],predict_CO2top[i],predict_VPair[i],predict_VPtop[i]])
# ========================================