import numpy as np



def euler(func, x0, xf, y0, h):
    x = x0
    #y = y0
    y = y0
    
    i = 0
    while x < xf:
        y = (y + h*func(x,y))
        i += 1
        x += h
        print(x,y)
    return y

def dydx(x,y):
    return 5*y

print(euler(dydx, 0, 1,0,0.25))

    



    
