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
    
def calcE(arr,N):#calculates energy
    E = 0
    m = 0
    for y in range(0, N):
        for x in range(0, N):
            '''
            	if(x<N-1):
			E=E -arr[x,y]*arr[x+1,y]
		elif(x==N-1): 
			E=E -arr[x,y]*arr[1,y]
		if(y<N-1):
			E=E -arr[x,y]*arr[x,y+1]
		elif(y==N-1):
			E=E -arr[x,y]*arr[x,1]
	'''
            m +=arr[x,y]
            E += -arr[x,y]*arr[bc(x+1,N),y]-arr[x,y]*arr[bc(x-1,N),y]-arr[x,y]*arr[x,bc(y+1,N)]-arr[x,y]*arr[x,bc(y-1,N)]
            
    return E/2, m
    
def relcalcE(arr,x,y,N):#calculates relative energy
    E = -arr[x,y]*arr[bc(x+1,N),y]-arr[x,y]*arr[bc(x-1,N),y]-arr[x,y]*arr[x,bc(y+1,N)]-arr[x,y]*arr[x,bc(y-1,N)]
    return E
    
def flip(arr,N,T):
        for y in range(0, N):
            for x in range(0, N):
                E0 = relcalcE(arr,x,y,N)
                Enew = -1*E0
                if(Enew<= E0):
                    arr[x,y] *= -1
                else:
                    p = exp(-(Enew-E0)/(T))
                    ran = random()
                    if(p>=ran):
                        arr[x,y] *= -1

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
        E0, arr,M = flip(arr,N,T)
        Eavg[i] +=E0
        Mavg[i] +=M

    Mavg[i] = (Mavg[i]/I)/(N**2)
    Eavg[i] = (Eavg[i]/I)/(N**2)
    i +=1

temp = arange(0.1,5.1,.1)

figure(1)
title('temp VS Eavg')
plot(temp,Eavg,'.',label='N = 50')
xlabel('Temperature')
ylabel('Eavg')
legend(loc = 'best')
figure(2)
title('temp VS Mavg')
plot(temp,Mavg,'.',label='N = 50')
xlabel('Temperature')
ylabel('Mavg')
legend(loc = 'best')


show()

timeend = time.clock()


runtime = timeend-timestart
print 'runtime: ', runtime



