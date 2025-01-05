### Mie_Scattering
#### Understanding the Physics of Mie Scattering

0. Mie scattering describes the scattering of electromagnetic waves (like light) by spherical particles. It is a solution to Maxwell's equations for the scattering of a plane wave by a homogeneous sphere. This theory is particularly useful when the size of the scattering particles is comparable to the wavelength of the incident light.

Key points:
Mie Theory: provides a comprehensive solution for scattering by spheres of any size.
Applications: Used in meteorology, astronomy, and optical engineering to understand phenomena like the color of the sky, radar detection etc.
Mathematical Formulation: Involves solving Maxwell's equations using boundary conditions at the surface of the sphere, resulting in an infinite series of spherical harmonics




1. Define constants and parameters:
   - Wavelength of incident light (λ)
   - Radius of the sphere (r)
   - Refractive index of the sphere (m)
   - Number of angles for scattering calculation (N)

2. Calculate size parameter (x):
   x = 2 * π * r / λ

3. Initialize arrays for scattering coefficients (a_n, b_n)

4. For each order n (from 1 to a maximum value):
   - Calculate spherical Bessel functions (j_n, y_n)
   - Calculate Riccati-Bessel functions (ψ_n, ξ_n)
   - Compute scattering coefficients a_n and b_n using Mie theory formulas

5. Initialize arrays for scattering intensity (I)

6. For each angle θ (from 0 to 180 degrees):
   - Calculate scattering amplitude functions (S1, S2)
   - Compute scattering intensity I(θ) using S1 and S2

7. Output scattering intensity as a function of angle