import sys

orig_stdout = sys.stdout
fp = open('Output.txt', 'w')
sys.stdout = fp

with open('Input_file.txt') as f:       #taking input data
    step = []
    for line in f:
        step+=([float(x) for x in line.split()])
k=0.06
# step=[0.1, 0.3, 0.5, 0.7, 0.8]
def f(y):
    return (-k*((y)**0.5))

def fd(y):
    return (0.5*k*k)

def fdd(y):
    return 0

def t_an(yi,yf):
    return 2*((yi)**(0.5)-(yf)**(0.5))/k
print("Euler's Method\n")
print("dt\ttime\t\trelative error\t\t\ttrun_error\t\tfun_evaluations\n")    
for dt in step: 
    count=0
    t=0
    y_old=3.0
    y=y_old
    while y>=0:
        y=y_old+f(y_old)*dt
        count+=2
        
        t+=dt 
        y_old=y
    print(dt,"\t",round(t,6),"\t\t",round(((t-t_an(3,0))*100/t_an(3,0)),6),"\t\t\t",round((fd(y)*(dt**2)/2),6),"\t\t\t",count)
print("\n\n")
print("Heun's Method\n")
print("dt\ttime\t\trelative error\t\t\ttrun_error\t\tfun_evaluations\n")  
for dt in step:   
    count=0
    tt=0
    y_old0=3.0
    y_hn=y_old0
    while y_hn>=0:
        y0=y_old0+f(y_old0)*dt
        count+=2
        if y0<0:
            tt+=dt
            break
        
        y_hn=y_old0+((f(y_old0)+f(y0))*dt/2.0)
        count+=4
        y_old0=y_hn
        tt+=dt
       
    print(dt,"\t",round(tt,6),"\t\t",round(((tt-t_an(3,0))*100/t_an(3,0)),6),"\t\t\t",round(fdd(y)*(dt**3)/6,6),"\t\t\t",count)
    
sys.stdout = orig_stdout
fp.close()

