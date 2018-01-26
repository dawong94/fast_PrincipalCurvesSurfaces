import numpy as np
import matplotlib.pyplot as plt
import Principal_Curves as pc
#import princurve2 as pc2
import time
"testing file for PCS "

def curve(n):
    "define random points around a sin curve"
    t = np.linspace(0.0, 4*np.pi, num = n)
    x = np.sin(t)+ np.random.normal(scale = 0.1, size = n)
    return t, x   

t, x = curve(300)
data = np.vstack([t, x])
points = np.copy(data) # so that manipulations don't change data
t1 = np.copy(data)
t2 = np.copy(data)
tol = (10.0)**(-2); d = 1
n, N = data.shape
h = N**(-1./(n+4)) # Scott's Rule

start_time = time.time()

curve= pc.principal_curve_surface(data, points, tol, d, h)


print("--- %s seconds ---" % (time.time() - start_time))

#start_time2 = time.time()
#curve2 = pc2.princurve(t1, t2, tol, d, h)
#print("--- %s seconds ---" % (time.time() - start_time2))

#error = np.linalg.norm(curve - curve2)
##print error
tmin = t.min(); tmax = t.max()
xmin = x.min(); xmax = x.max()

fig = plt.figure()
ax = fig.add_subplot(111)
#ax.plot(t, x, 'k.', markersize=5)
#ax.plot(curve2[0, :], curve2[1, :], 'b')
ax.set_xlim([-2, 14])
ax.set_ylim([-2, 2])
#ax.plot(t, x, 'k.', markersize=5)
ax.plot(curve[0, :], curve[1, :], 'r')

ax.set_xlim([-2, 14])
ax.set_ylim([-2, 2])
plt.show()
