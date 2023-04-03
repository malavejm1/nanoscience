from mpl_toolkits import mplot3d
import numpy as np
import matplotlib.pyplot as plt



def f(x, y):
    return np.sqrt(x ** 2) + np.sqrt(y ** 2)

x = np.linspace(0, 1, 30)
y = np.linspace(0, 1, 30)

X, Y = np.meshgrid(x, y)
Z = f(X, Y)

fig = plt.figure()
ax = plt.axes(projection='3d')
ax.contour3D(X, Y, Z, 50, cmap='binary')
ax.set_xlabel('x')
ax.set_ylabel('y')
ax.set_zlabel('z')
plt.show()



