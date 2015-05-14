from scipy.integrate import simps
from scipy.integrate import trapz
from scipy.integrate import romb
import numpy as np

def f(x):
    """
    """
    return 1

def g(x):
    """
    Function taken from https://www.math.duke.edu/vigre/pruv/studentwork/atwood.nearsing.pdf
    """

    a = 10e-3

    return f(x)* np.log(x**2 + a**2)

def simpson(f, a, b, n):
    """
    Simpson composite rule
    """
    h=(b-a)/n
    k=0.0
    x=a + h
    for i in range(1,n/2 + 1):
        k += 4.*f(x)
        x += 2.*h

    x = a + 2.*h
    for i in range(1,n/2):
        k += 2.*f(x)
        x += 2.*h
    return (h/3.)*(f(a)+f(b)+k)



# Computing integral using less interval
# Correct answer is -3.99372

I = -3.99372

i1 = simpson(g, -1., 1., 2**5+ 1)

i2 = simpson(g, -1., 1., 2**10+1)


print 'Reference answer is ' + str(I)
print 'Previous code answer is ' + str(i1)

print 'Current code answer is ' + str(i2)

print 'Previous code absolute error is ' + str(np.abs(I-i1))

print 'Current code absolute error is ' + str(np.abs(I-i2))


