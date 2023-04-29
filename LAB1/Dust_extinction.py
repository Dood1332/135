import numpy as np
import matplotlib.pyplot as plt
import astropy.units as u
from pylab import *

from dust_extinction.averages import GCC09_MWAvg


# generate the curves and plot them
x = np.arange(0.3,10.0,0.1)/u.micron
ext_model = GCC09_MWAvg()

plot(1./x,ext_model(x))
plot(0.550,1, 'o',label='A($\lambda$) = A(V)')



title('Normalized Extinction')
xlabel('$\lambda$ [$\mu m$]')
ylabel('$A(\lambda)/A(V)$')
xscale('log')
xlim(0.09,4.0)
legend(loc='best')
tight_layout()
show()