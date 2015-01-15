from pylab import *
f=genfromtxt('/home/justin/IsingEavg.txt')
Eavg = f

f=genfromtxt('/home/justin/IsingMavg.txt')
Mavg = f


temp = arange(0.,5.1,.1)
print size(temp)
figure(1)
title('temp VS Eavg')
plot(temp,Eavg,'.')

temp = arange(0.0,5.1,.1)

figure(2)
title('temp VS Mavg')
plot(temp,Mavg,'.')
show()