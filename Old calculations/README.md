# harmonic_oscillator_s-wave_contact
Various functions for calculating and plotting s-wave contact for two atoms in a HO


BuschFunc consists of various functions:
- Convert E, a0, B
- Calculate pseudopotential + HO wavefunctions
- integrands for numerical integrals
- interpolated contact in upper+ lower branch vs. B

Contact.py calculates contact through adiabatic relation and fourier transform.
- Fourier transform is very slow
- Change energies E_C to calculate for different branches
- Outputs to file C_array.csv

C_array.csv and C_array_lower.csv have tabulated contact values:
- Column 1: E (units: hbar*omega)
- Column 2: C from adiatic relation (numerical derivative)
- Column 3: C from fourier transform
- Column 4: lim r-> (r*psi(r)). C = 16*pi^2*(lim r-> (r*psi(r)))

In general, all units are harmonic oscillator units, unless converted to B-field, which is in Gauss, etc.

