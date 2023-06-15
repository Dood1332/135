import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import argrelextrema
x = np.linspace(-10, 10, 1000)
# Example y array representing a double Gaussian curve
y = np.exp(-0.5 * ((x - 2) / 1.5) ** 2) + np.exp(-0.5 * ((x + 2) / 1.5) ** 2)

# Plotting the curve
plt.plot(x, y, label='Double Gaussian')
plt.legend()
plt.xlabel('x')
plt.ylabel('y')
plt.show()

# Find the indices of local maxima
maxima_indices = argrelextrema(y, np.greater)

# Get the x-values corresponding to the maxima indices
x_maxima = x[maxima_indices]

# Sort the x-values in descending order of y-values
sorted_indices = np.argsort(y[maxima_indices])[::-1]
x_peaks = x_maxima[sorted_indices][:2]

print("Peak x-values:", x_peaks)
