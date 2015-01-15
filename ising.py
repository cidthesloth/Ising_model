from pylab import *
import time

def bc(x,N):
    if x+1> N-1:
        return 0
    if x-1<0:
        return N-1
    else:
        return x
def createarr(N): #creates random array of 1 and -1
    arr = zeros((N,N))
    for y in range(0, N):
        for x in range(0, N):
            r = randint(0, 2)
            if (r == 0):
                arr[x,y]= -1
            else:
                arr[x,y]= 1
    return arr
    
def copy(arr,N): #copys the arr to a new array
    arrnew = zeros((N,N))
    for y in range(0, N):
        for x in range(0, N):
            arrnew[x,y]=arr[x,y]
    return arrnew
def calcE(arr,N):#calculates energy
    E = 0
    m = 0
    for y in range(0, N):
        for x in range(0, N):
            m +=arr[x,y]
            E += -arr[x,y]*arr[bc(x+1,N),y]-arr[x,y]*arr[bc(x-1,N),y]-arr[x,y]*arr[x,bc(y+1,N)]-arr[x,y]*arr[x,bc(y-1,N)]
    return E/2, m
    
def relcalcE(arr,x,y,N):#calculates relative energy
    E = -arr[x,y]*arr[bc(x+1,N),y]-arr[x,y]*arr[bc(x-1,N),y]-arr[x,y]*arr[x,bc(y+1,N)]-arr[x,y]*arr[x,bc(y-1,N)]
    return E
    
def flip(arr,N,T,kb):
        for y in range(0, N):
            for x in range(0, N):
                E0 = relcalcE(arr,x,y,N)
                arrnew = copy(arr,N)
                if(arrnew[x,y] ==1.):
                    arrnew[x,y] = -1.
                    Enew = relcalcE(arrnew,x,y,N)
                else:
                    arrnew[x,y] = 1.
                    Enew = relcalcE(arrnew,x,y,N)
                if(Enew<= E0):
                    arr = copy(arrnew,N)
                else:
                    p = exp(-(Enew-E0)/(kb*T))
                    x = random()
                    if(p>=x):
                        arr = copy(arrnew,N)

        E0, m = calcE(arr,N) 
        return E0,arr,m
                
              
#constants
kb = 1 #1.38*10**-23    
                       
#declaring variables

N = input('Input the size of square array: ')
I = input('Input the number of iterations: ')
timestart = time.clock()
Eavg = zeros(50)
Mavg = zeros(50)
arr = createarr(N)
i = 0
for y in range(1,51,1):
    T = y/10.
    for x in range(0, I):       
        E0, arr,M = flip(arr,N,T,kb)
        Eavg[i] +=E0
        Mavg[i] +=M

    Mavg[i] = (Mavg[i]/I)/(N**2)
    Eavg[i] = (Eavg[i]/I)/(N**2)
    i +=1

temp = arange(0.1,5.1,.1)

figure(1)
title('temp VS Eavg')
plot(temp,Eavg,'.')
figure(2)
title('temp VS Mavg')
plot(temp,Mavg,'.')


show()

timeend = time.clock()


runtime = timeend-timestart
print runtime


