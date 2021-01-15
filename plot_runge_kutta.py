import numpy as np
import matplotlib.pyplot as plt


def dydt(t, y):
    return -t / (y - 3)

def rungeKutta(t0, t, y0, h):
    # Count number of iterations using step size or
    # step height h
    n = int((t - t0) / h)
    # Iterate for number of iterations
    t_result = []
    y_result = []
    t_result.append(t0)
    y_result.append(y0)
    for i in range(1, n + 1):
        k1 = dydt(t0, y0)
        k2 = dydt(t0 + 0.5 * h, y0 + 0.5 * k1 * h)
        k3 = dydt(t0 + 0.5 * h, y0 + 0.5 * k2 * h)
        k4 = dydt(t0 + h, y0 + k3 * h)

        # Update next value of y
        y0 = y0 + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4) * h
        y_result.append(y0)
        # Update next value of t
        t0 = t0 + h
        t_result.append(t0)
    return t_result, y_result

t_exact = np.linspace(0,2,100)
y_exact = 3 - (4-t_exact**2)**0.5


t_1, y_1 = rungeKutta(0,2,1,0.4)
t_2, y_2 = rungeKutta(0,2,1,0.2)
t_3, y_3 = rungeKutta(0,2,1,0.1)


plt.plot(t_exact,y_exact,'-', linewidth=2)
plt.plot(t_1,y_1,'-',linewidth=2)
plt.plot(t_2,y_2,'-',linewidth=2)
plt.plot(t_3,y_3,'-',linewidth=2)
plt.xlabel('t')
plt.ylabel('y')
plt.title('Differential Equation: dy/dt=-t/(y-3), y(0) = 1 ')
plt.legend(["Exact Solution","h = 0.4","h = 0.2","h = 0.1"])
plt.show()
