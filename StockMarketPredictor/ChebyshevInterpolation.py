from sympy import *
import math
import numpy as np
import matplotlib.pyplot as plt
from NumericalMethods import *

#set x as a variable
x = symbols('x')

#interpolation_points = [(960.0, 38472.19), (1080.0, 38538.91), (1200.0, 38686.55), (1560.0, 38818.78), (1680.0, 38543.51), (1800.0, 38953.84), (1920.0, 38680.46), (2040.0, 38371.27), (2160.0, 38548.26), (2400.0, 38528.11), (3480.0, 38112.55), (3600.0, 37727.88), (3840.0, 37748.01), (3960.0, 37886.57), (4080.0, 38026.2), (4560.0, 39003.07), (4680.0, 38760.81), (4800.0, 38962.93), (4920.0, 39085.77), (5040.0, 39822.46), (5160.0, 39759.5), (5280.0, 39699.27), (5400.0, 39644.54), (5520.0, 39723.92), (6480.0, 36304.72), (6600.0, 36484.32), (6720.0, 36573.18), (6840.0, 36316.97), (7080.0, 36494.05), (7200.0, 36370.44), (7680.0, 36035.4), (7920.0, 36013.06), (8040.0, 35983.9), (8160.0, 36042.5), (9240.0, 35904.09), (9360.0, 35885.67), (9720.0, 34884.45), (9840.0, 34505.46), (9960.0, 34451.61)]
#interpolation_points = [(960.0, 38472.19), (1080.0, 38538.91), (1200.0, 38686.55), (1560.0, 38818.78)]
interpolation_points = [(960.0, 38472.19), (1080.0, 38538.91), (1200.0, 38686.55), (1560.0, 38818.78), (1680.0, 38543.51), (1800.0, 38953.84), (1920.0, 38680.46), (2040.0, 38371.27), (2160.0, 38548.26), (2400.0, 38528.11)]
#parse the interpolation points into x and y values
xi = []
fxi = []
for i in interpolation_points:
    xi.append(i[0])
    fxi.append(i[1])

print("x values: {0}".format(xi))
print("y values: {0}".format(fxi))
print("====="*10)

def chebyshev(n, x):
    if n == 0:
        return 1
    elif n == 1:
        return x
    else:
        return 2*x*chebyshev(n-1, x) - chebyshev(n-2, x)

def nodes(a, b, n):
    x_nodes = []
    for k in range(n+1):
        #this method only calculates nodes from [-1, 1]
        #node = math.cos(((2*k+1)*math.pi)/(2*(n+1)))

        #this calculates nodes from [a, b]
        x_k = math.cos(((2*k+1)*math.pi)/(2*(n+1)))
        node = a+0.5*(b-a)*(x_k+1)
        x_nodes.append(node)

    return x_nodes

"""
#print(nodes(4))
print(nodes(-1, 1, 1))
def f(x):
    return math.exp(x)
cheb_nodes = nodes(1, 5, 4)
eqn, x_vals, y_vals = divided_differences(cheb_nodes, [f(i) for i in cheb_nodes], True, False)
print(simplify(eqn))
"""

plt.scatter(xi, fxi, label="Raw Data")
x_lim = plt.xlim()
y_lim = plt.ylim()

eqn, x_vals, y_vals = divided_differences(xi, fxi, False, False)
print("Lagrange equation:\n{0}".format(simplify(eqn)))
plt.plot(x_vals, y_vals, label="Lagrange", color='blue')
print("====="*10)

n=len(fxi)-1
cheb_nodes = nodes(xi[0], xi[-1], n)
cheb_nodes.sort()
eqn, x_vals, y_vals = divided_differences(cheb_nodes, [calc(eqn, i) for i in cheb_nodes], False, False)
print("Chebyshev equation:\n{0}".format(simplify(eqn)))
plt.plot(x_vals, y_vals, label="Chebyshev", color='green')

plt.xlim(x_lim)
plt.ylim(y_lim)
plt.title("Price of Bitcoin 5/1/2022 8:00 am - 5/7/2022 4:00 am")
plt.ylabel("Price USD")
plt.xlabel("Mintues Since 5/1/2022 8:00 am")
plt.grid(which='major', axis='both')
plt.legend(loc="upper left")
plt.show()



"""
eqn = chebyshev(3, x)
eqn = simplify(eqn)
print(eqn)


#Chebyshev interpolate with numpy
cheb = np.polynomial.chebyshev.chebfit(xi, fxi, len(xi)-1)
print(cheb)
#print("Chebyshev Polynomial: " + str(cheb[3]) + "x^3 + " + str(cheb[2]) + "x^2 + " + str(cheb[1]) + "x + " + str(cheb[0]))

string = ""
for x in range(len(cheb)):
    coeff = str(cheb[x])
    
    if "e" in coeff:
        coeff = coeff.split("e")[0]+"^{"+coeff.split("e")[1]+"}"

    #print(coeff)
    #string += coeff+"x^{"+str(len(cheb)-x-1)+"}+"
    string += coeff+"x^{"+str(x)+"}*cos("+str(x)+"arccos(x))+"

print(string)
"""
