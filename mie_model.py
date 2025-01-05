import numpy as np
import scipy.special as sp

def mie_scattering(wavelength, radius, refractive_index, num_angles):
    
    # Constants
    k = 2 * np.pi / wavelength
    x = k * radius
    m = refractive_index

    # Size parameter
    size_parameter = x

    # Scattering coefficients
    a_n = []
    b_n = []

    # Calculate scattering coefficients
    for n in range(1, 100):  # Adjust the range as needed
        psi_n = sp.spherical_jn(n, size_parameter)
        psi_n_prime = sp.spherical_jn(n, size_parameter, derivative=True)
        xi_n = psi_n + 1j * sp.spherical_yn(n, size_parameter)
        xi_n_prime = psi_n_prime + 1j * sp.spherical_yn(n, size_parameter, derivative=True)

        a_n.append((m * psi_n * psi_n_prime - psi_n * xi_n_prime) / (m * psi_n * xi_n_prime - xi_n * psi_n_prime))
        b_n.append((psi_n * psi_n_prime - m * psi_n * xi_n_prime) / (psi_n * xi_n_prime - m * xi_n * psi_n_prime))

    # Scattering intensity
    angles = np.linspace(0, np.pi, num_angles)
    intensity = np.zeros_like(angles)

    for i, theta in enumerate(angles):
        S1 = 0
        S2 = 0
        for n in range(1, len(a_n)):
            pi_n = sp.lpmv(0, n, np.cos(theta))
            tau_n = np.sin(theta) * sp.lpmv(1, n, np.cos(theta))

            S1 += (2 * n + 1) / (n * (n + 1)) * (a_n[n] * pi_n + b_n[n] * tau_n)
            S2 += (2 * n + 1) / (n * (n + 1)) * (a_n[n] * tau_n + b_n[n] * pi_n)

        intensity[i] = (np.abs(S1)**2 + np.abs(S2)**2) / 2

    return angles, intensity
