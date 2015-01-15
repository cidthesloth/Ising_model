from pylab import *
f1=genfromtxt('/home/justin/IsingEavg1.txt')
f2=genfromtxt('/home/justin/IsingEavg2.txt')
f3=genfromtxt('/home/justin/IsingEavg3.txt')
f4=genfromtxt('/home/justin/IsingEavg4.txt')
f = f1+f2+f3+f4
Eavg = f

f1=genfromtxt('/home/justin/IsingMavg1.txt')
f2=genfromtxt('/home/justin/IsingMavg2.txt')
f3=genfromtxt('/home/justin/IsingMavg3.txt')
f4=genfromtxt('/home/justin/IsingMavg4.txt')
f = f1+f2+f3+f4
Mavg = f


temp = arange(0.1,5.1,.1)

figure(1)
title('temp VS Eavg')
plot(temp,Eavg,'.')


figure(2)
title('temp VS Mavg')
plot(temp,Mavg,'.')
show()