import matplotlib.pyplot as plt
import numpy as np

y4 = data4[:,1]
data5 = np.loadtxt("bs4_tq_plk1/FE_bs4_tq.dat")
x5 = data5[:,0]
y5 = data5[:,1]
data6 = np.loadtxt("bs5_tq/FE_bs5_tq.dat")
x6 = data6[:,0]
y6 = data6[:,1]
data7 = np.loadtxt("bs6_tq/FE_bs6-tq.dat")
x7 = data7[:,0]
y7 = data7[:,1]
data8 = np.loadtxt("bs7_tq/FE_bs7_tq.dat")
x8 = data8[:,0]
y8 = data8[:,1]
data3 = np.loadtxt("../bs3_tq1/FE_bs3_tq.dat")
x3 = data3[:,0]
y3 = data3[:,1]
#plt.plot(x,y, linewidth=2, label = '0.1')
#plt.plot(x2,y2, linewidth=2, label = '0.2')
#plt.plot(x3,y3, linewidth=2, label = 'wt_0.1')
plt.plot(x4,y4, linewidth=1.5, label = 'BS1', color = 'red')
plt.plot(x3,y3, linewidth=1.5, label = 'BS3', color = 'green')
plt.plot(x5,y5, linewidth=1.5, label = 'BS4', color = 'orange')
plt.plot(x6,y6, linewidth=1.5, label = 'BS5', color = 'cyan')
plt.plot(x7,y7, linewidth=1.5, label = 'BS6', color = 'brown')
plt.plot(x8,y8, linewidth=1.5, label = 'BS7', color = 'blue')
#plt.plot(x10,y10, linewidth=2, label = '10ns')
plt.xlabel('COM-Distance (Å)')
plt.legend()
plt.xlim(0,30)
#plt.ylim(-8,0)
plt.grid(visible=True, which='major', axis='both', color='#808080', linestyle='--')
plt.ylabel('PMF (kcal/mol)')
plt.savefig('FE.png', dpi=200)
plt.show()
