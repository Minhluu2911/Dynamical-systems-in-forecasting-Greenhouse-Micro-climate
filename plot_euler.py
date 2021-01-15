import numpy as np
import matplotlib.pyplot as plt

# Python Code to find approximation 
# of a ordinary differential equation 
# using euler method. 

def func( t, y ): 
	return -t/(y-3)
	
# Function for euler formula 
def euler( t0, t, y0, h ):
    t_result = []
    y_result = []
    t_result.append(t0)
    y_result.append(y0)
    # Iterating till the point at which we 
    # need approximation 
    for _ in range(int(t/h)):
        y0 = y0 + h * func(t0, y0) 
        y_result.append(y0)
        t0 = t0 + h 
        t_result.append(t0)
    
    return t_result,y_result

t_exact = np.linspace(0,2,100)
y_exact = 3 - (4-t_exact**2)**0.5


t_1, y_1 = euler(0,2,1,0.4)
t_2, y_2 = euler(0,2,1,0.2)
t_3, y_3 = euler(0,2,1,0.1)

plt.plot(t_exact,y_exact,'-', linewidth=2)
plt.plot(t_1,y_1,'-',linewidth=2)
plt.plot(t_2,y_2,'-',linewidth=2)
plt.plot(t_3,y_3,'-',linewidth=2)
plt.xlabel('t')
plt.ylabel('y')
plt.title('Differential Equation: dy/dt=-t/(y-3), y(0) = 1 ')
plt.legend(["Exact Solution","h = 0.4","h = 0.2","h = 0.1"])
plt.show()
