use std::f64::consts::PI;
use ndarray::Array1;
use special::SphericalBessel;
use num_complex::Complex;

fn mie_scattering(wavelength: f64, radius: f64, refractive_index: f64, num_angles: usize) -> (Array1<f64>, Array1<f64>) {
    // Constants
    let k = 2.0 * PI / wavelength;
    let x = k * radius;
    let m = refractive_index;

    // Size parameter
    let size_parameter = x;

    // Scattering coefficients
    let mut a_n = Vec::new();
    let mut b_n = Vec::new();

    // Calculate scattering coefficients
    for n in 1..100 {
        let psi_n = SphericalBessel::jn(n, size_parameter);
        let psi_n_prime = SphericalBessel::jn_prime(n, size_parameter);
        let xi_n = Complex::new(psi_n, SphericalBessel::yn(n, size_parameter));
        let xi_n_prime = Complex::new(psi_n_prime, SphericalBessel::yn_prime(n, size_parameter));

        let a_n_value = (m * psi_n * psi_n_prime - psi_n * xi_n_prime) / (m * psi_n * xi_n_prime - xi_n * psi_n_prime);
        let b_n_value = (psi_n * psi_n_prime - m * psi_n * xi_n_prime) / (psi_n * xi_n_prime - m * xi_n * psi_n_prime);

        a_n.push(a_n_value);
        b_n.push(b_n_value);
    }

    // Scattering intensity
    let angles = Array1::linspace(0.0, PI, num_angles);
    let mut intensity = Array1::zeros(num_angles);

    for (i, &theta) in angles.iter().enumerate() {
        let mut S1 = Complex::new(0.0, 0.0);
        let mut S2 = Complex::new(0.0, 0.0);
        for n in 1..a_n.len() {
            let pi_n = special::legendre::P(n as i32, 0, theta.cos());
            let tau_n = theta.sin() * special::legendre::P(n as i32, 1, theta.cos());

            S1 += (2 * n + 1) as f64 / (n * (n + 1)) as f64 * (a_n[n] * pi_n + b_n[n] * tau_n);
            S2 += (2 * n + 1) as f64 / (n * (n + 1)) as f64 * (a_n[n] * tau_n + b_n[n] * pi_n);
        }
        intensity[i] = (S1.norm_sqr() + S2.norm_sqr()) / 2.0;
    }

    (angles, intensity)
}