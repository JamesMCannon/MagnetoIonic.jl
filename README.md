# MagnetoIonic

[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://JamesMCannon.github.io/MagnetoIonic.jl/stable/)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://JamesMCannon.github.io/MagnetoIonic.jl/dev/)
[![Build Status](https://github.com/JamesMCannon/MagnetoIonic.jl/actions/workflows/CI.yml/badge.svg?branch=main)](https://github.com/JamesMCannon/MagnetoIonic.jl/actions/workflows/CI.yml?query=branch%3Amain)
[![Coverage](https://codecov.io/gh/JamesMCannon/MagnetoIonic.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/JamesMCannon/MagnetoIonic.jl)

MagnetoIonic.jl serves to aid in the computation of radio wave propagation through plasma. 

As it currently stands, this package can calculate the squared index of refraction using the Appleton-Hartree equation in either Hertz or radians from `appleton_hartree_rad()` and `appleton_hartree_hz()`. 

Multiple input types are supported for both but all functions need a provided electron density `n_e` in m⁻³, effective electron-neutral collision frequency in s⁻¹, magnetic field and wave propogation as either as a 3 component vector `B_vec` in T with associated wave `k` vector or as magentitude `B` in T and angle `θ` between the wave normal and magnetic field in radians and radio frequency in either radians s⁻¹ or Hertz.

When supplying `B_vec`, `k` is optional and, if not included, the assumed direction of wave propagation is vertical [0, 0, 1].

Additional functions are provided for calculating the associated absorption coeficients for both O and X mode waves given the same information using `radio_absorption_rad()` and `radio_absorption_hz()`. Returned coeficients, obtained from the imaginary part of the refractive index `n = √(n²)` are in units of dB (loss) m⁻¹

Convenience functions for common plasma calculations are exported including `plasma_frequency_(rad/hz)()`, `plasma_frequency_(rad/hz)_sq()`, `gyrofrequency_(rad/hz)()`, and `gyrofrequency_(rad/hz)_sq()`.

Under the hood, all computations are converted to radians first. For many repeated calls, consider working natively in radians. 

See Zawdie2017 for a discussion on the use of different collision frequencies in application to ionospheric plasma calculations.

# References
- [Appleton1932]: Appleton, E. V. (1932), Wireless studies of the ionosphere, Int. J. Electr. Eng., 71, 642–650.
- [Budden1985]: Budden, K. G. (1985), The propagation of radio waves, Cambridge Univ. Press, Cambridge.
- [Haselgrove1960]: Haselgrove, C. B., and J. Haselgrove (1960), Twisted ray paths in the ionosphere, Proc. Phys. Soc. London, 75, 357–363.
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.
