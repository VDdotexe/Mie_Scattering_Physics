from mie_model import mie_scattering

import matplotlib.pyplot as plt

# Example usage
wavelength = 0.5  # micrometers (visible for now)
radius = 0.1  # micrometers
refractive_index = 1.5
num_angles = 180

angles, intensity = mie_scattering(wavelength, radius, refractive_index, num_angles)


plt.polar(angles, intensity)
plt.title('Mie Scattering Intensity')
plt.show()