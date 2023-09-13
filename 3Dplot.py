import numpy as np
import matplotlib.pyplot as plt

plt.figure()

ax = plt.axes(projection ='3d')

x = np.linspace(.1, 3, 201)
y = np.linspace(.1, 3, 201)
x,y = np.meshgrid (x,y)
z = ((1-x**2)**2+x**2/(y**2))**(-1/2)

ax.plot_surface( x,y,z, )

plt.xlabel('$x$')
plt.ylabel('$Q$')


plt.show()
