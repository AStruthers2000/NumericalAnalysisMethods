import sys

import matplotlib.pyplot as plt

from NumericalMethods import *


def run():
    methods, xi, fxi = gather_data()

    plt.scatter(xi, fxi, label="Scatter of Data", color='black')
    x_lim = plt.xlim()
    y_lim = plt.ylim()

    equation = divided_differences(xi, fxi, False, False)[0]

    for method in methods:
        if method == "m":
            p, r, e = taylor(equation, xi, fxi, 0, 3)
            plt.plot(xi, e, label="Maclaurin")
            pass

        elif method == "t":
            print("Taylor method not implemented")

        elif method == "l":
            print("Lagrange not implemented yet")

        elif method == "ldd":
            print("Divided Differences Method:")
            equation, x_vals, y_vals = divided_differences(xi, fxi, True, False)
            plt.plot(x_vals, y_vals, label="Lagrange - Divided Difference", color='red')
            print("=-=-=-"*6+"=")

        elif method == "h":
            print("Hermite Method:")
            equation, x_vals, y_vals = hermite(xi, fxi)
            plt.plot(x_vals, y_vals, label="Hermite", color='green')
            print("=-=-=-"*6+"=")

        elif method == "s":
            print("Cubic Spline:")
            spline(xi, fxi)
            print("=-=-=-" * 6 + "=")

        elif method == "nd":
            print("Numerical Derivative:")
            fprime_est = estimate(fxi, h=xi[1]-xi[0])
            plt.scatter(xi, fprime_est)
            print(fprime_est)
            

        else:
            print("Method not defined")

    #plt.xlim(x_lim)
    #plt.ylim(y_lim)
    plt.xlim((-0.5,9))
    plt.ylim((-3,2))
    plt.grid(which='major', axis='both')
    plt.legend(loc="upper left")
    img = plt.imread("WhaleFlukeData.png")
    plt.imshow(img, extent=[-0.5+0.25, 9+0.1, -3-0.15, 2+0.35])
    plt.show()


def gather_data():
    xi = get_x_values()
    fxi = get_y_values(xi)

    if not xi:
        xi = [i for i in range(len(fxi))]

    desired_methods = input("Which methods do you want? (m,t,nd,l,ldd,h,s) ").replace(" ", "").split(",")

    data = [[] for _ in range(len(xi))]
    for i in range(len(xi)):
        data[i] = [xi[i], fxi[i]]
    print(tabulate(data, headers=["x", "F(x)"]))

    print("desired methods: {0}".format(desired_methods))
    return desired_methods, xi, fxi


def get_x_values():
    try:
        xi = input("Enter x values separated by comma, path and name of datafile, or nothing: ")
        if xi == "":
            pass
        elif "," in xi:
            xi = [float(f) for f in xi.split(",")]
        else:
            with open(xi, 'r') as file:
                xi = []
                for line in file:
                    xi.append(float(line.strip("\n").replace(" ", "")))
    except FileNotFoundError:
        print("File not found")
        sys.exit()
    except ValueError:
        print("Bad dataset")
        sys.exit()
    return xi


def get_y_values(xi):
    try:
        fxi = input("Enter y values separated by comma, or path and name of datafile: ")
        # f = sp.sympify(fxi)
        if "," in fxi:
            fxi = [float(f) for f in fxi.split(",")]
        else:
            try:
                with open(fxi, 'r') as file:
                    fxi = []
                    for line in file:
                        if line.strip("\n"):
                            fxi.append(float(line.strip("\n").replace(" ", "")))
                
            except FileNotFoundError:
                f = sp.sympify(fxi)
                x = sp.symbols('x')
                fxi = [f.subs(x, i) for i in xi]
                print(f)
                print(fxi)
    except FileNotFoundError:
        print("File not found")
        sys.exit()
    except ValueError:
        print("Bad dataset")
        sys.exit()
    return fxi


if __name__ == "__main__":
    run()
