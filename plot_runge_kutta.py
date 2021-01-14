import numpy as np
import matplotlib.pyplot as plt


def dydx(x, y):
    return -x / (y - 3)


def rungeKutta(x0, x, y0, h):
    # Count number of iterations using step size or
    # step height h
    n = (int)((x - x0) / h)
    # Iterate for number of iterations
    y = y0
    x_result = []
    y_result = []
    x_result.append(x0)
    y_result.append(y)
    for i in range(1, n + 1):
        k1 = dydx(x0, y)
        k2 = dydx(x0 + 0.5 * h, y + 0.5 * k1 * h)
        k3 = dydx(x0 + 0.5 * h, y + 0.5 * k2 * h)
        k4 = dydx(x0 + h, y + k3 * h)

        # Update next value of y
        y = y + (1.0 / 6.0) * (k1 + 2 * k2 + 2 * k3 + k4) * h
        y_result.append(y)
        # Update next value of x
        x0 = x0 + h
        x_result.append(x0)
    return x_result, y_result

x_exact = np.linspace(0,2,100)
y_exact = 3 - (4-x_exact**2)**0.5


x_1, y_1 = rungeKutta(0,2,1,0.4)
x_2, y_2 = rungeKutta(0,2,1,0.2)
x_3, y_3 = rungeKutta(0,2,1,0.1)


plt.plot(x_exact,y_exact,'-', linewidth=2)
plt.plot(x_1,y_1,'-',linewidth=2)
plt.plot(x_2,y_2,'-',linewidth=2)
plt.plot(x_3,y_3,'-',linewidth=2)
plt.xlabel('xlabel')
plt.ylabel('ylabel')
plt.title('Differential Equation: dy/dx=x/(y-3), y(0) = 1 ')
plt.legend(["Exact Solution","h = 0.4","h = 0.2","h = 0.1"])
plt.show()
