import numpy as np
import matplotlib.pyplot as plt

# Python Code to find approximation 
# of a ordinary differential equation 
# using euler method. 

def func( x, y ): 
	return -x/(y-3)
	
# Function for euler formula 
def euler( x0, x, y0, h ): 
    temp = -0
    
    x_result = []
    y_result = []
    x_result.append(x0)
    y_result.append(y0)
    # Iterating till the point at which we 
    # need approximation 
    while x0 - x < 0.0000000001: 
        temp = y0 
        y0 = y0 + h * func(x0, y0) 
        y_result.append(y0)
        x0 = x0 + h 
        x_result.append(x0)
    
    return x_result,y_result

x_exact = np.linspace(0,2,100)
y_exact = 3 - (4-x_exact**2)**0.5


x_1, y_1 = euler(0,2,1,0.4)
x_2, y_2 = euler(0,2,1,0.2)
x_3, y_3 = euler(0,2,1,0.1)

plt.plot(x_exact,y_exact,'-', linewidth=2)
plt.plot(x_1,y_1,'-',linewidth=2)
plt.plot(x_2,y_2,'-',linewidth=2)
plt.plot(x_3,y_3,'-',linewidth=2)
plt.xlabel('xlabel')
plt.ylabel('ylabel')
plt.title('Differential Equation: dy/dx=x/(y-3), y(0) = 1 ')
plt.legend(["Exact Solution","h = 0.4","h = 0.2","h = 0.1"])
plt.show()