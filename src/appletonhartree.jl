"""
    _appleton_hartree_kernel(X, Y, Z, sinθ, cosθ)

Evaluate the Appleton-Hartree equation from precomputed dimensionless
parameters and return the squared refractive indices of the two magnetoionic
modes.

# Arguments
- `X`: squared plasma-to-wave frequency ratio `ωp²/ω²`
- `Y`: gyro-to-wave frequency ratio `ωH/ω`
- `Z`: collision-to-wave frequency ratio `ν/ω`
- `sinθ`, `cosθ`: sine and cosine of the angle between the wave normal and the magnetic field

# Returns
- Tuple `(O, X)` of the complex squared refractive indices for the `+` and `-`
  roots of the discriminant, respectively

# Notes
- Internal helper; not part of the public API.
- The collision term enters through `U = 1 + iZ`, so the result is complex and
  its imaginary part carries the collisional damping.
- `YL = Y cosθ` and `YT = Y sinθ` are the longitudinal and transverse components
  of the magnetoionic `Y` parameter.
- The `O` and `X` fields correspond to the `+` and `-` signs of the discriminant.
  Their identification with the physical ordinary and extraordinary modes is
  regime-dependent and can interchange across the reflection condition `X = 1`;
  callers needing a guaranteed physical labeling should check the regime
  explicitly.
"""
@inline function _appleton_hartree_kernel(X, Y, Z, sinθ, cosθ)
    YL  = Y * cosθ
    YT  = Y * sinθ
    YT² = YT * YT

    U  = complex(one(Z), Z)              
    UX = U - X
    Δ  = sqrt(YT² * YT² / 4 + YL^2 * UX^2)   

    term = U - YT² / (2 * UX)            
    disc = Δ / UX                        
    n²_O = 1 - X / (term + disc)
    n²_X = 1 - X / (term - disc)
    return (O = n²_O, X = n²_X)
end


"""
    _magnetoionic_parameters(Ne, ν, B, ω)

Return the dimensionless magnetoionic parameters `(X, Y, Z)` used by the
Appleton-Hartree equation.

# Arguments
- `Ne`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B`: magnetic field strength in T
- `ω`: wave angular frequency in rad s⁻¹

# Returns
- Tuple `(X, Y, Z)` where
    - `X = ωp² / ω²` is the squared plasma-to-wave frequency ratio
    - `Y = ωH / ω` is the gyro-to-wave frequency ratio
    - `Z = ν / ω` is the collision-to-wave frequency ratio

# Notes
- Internal helper; not part of the public API.
"""
function _magnetoionic_parameters(Ne, ν, B, ω)
    X = plasma_frequency_rad_sq(Ne) / ω^2   
    Y = gyrofrequency_rad(B) / ω
    Z = ν / ω
    return X, Y, Z
end

"""
    appleton_hartree_rad(n_e, ν, B, θ, ω)

Return the squared index of refraction for the O and X modes as calculated using the Appleton-Hartree equation.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B`: magnetic field strength in T
- `θ`: angle between the wave normal and the magnetic field in radians
- `ω`: wave angular frequency in rad s⁻¹

# Returns
- Tuple `(O, X)` where `O` and `X` are the squared indices of refraction for the O and X modes, respectively

# Notes
- This is the ''full'' Appleton-Hartree equation, including collisions and valid for any angle θ.
- The squared indices of refraction are returned directly, as they are the natural output of the core calculation.
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Appleton1932]: Appleton, E. V. (1932), Wireless studies of the ionosphere, Int. J. Electr. Eng., 71, 642–650.
- [Budden1985]: Budden, K. G. (1985), The propagation of radio waves, Cambridge Univ. Press, Cambridge.
- [Haselgrove1960]: Haselgrove, C. B., and J. Haselgrove (1960), Twisted ray paths in the ionosphere, Proc. Phys. Soc. London, 75, 357–363.
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`appleton_hartree_hz`](@ref): equivalent function for Hz input
"""
function appleton_hartree_rad(n_e, ν, B::Real, θ::Real, ω::Real)
    X, Y, Z = _magnetoionic_parameters(n_e, ν, B, ω)
    s, c = sincos(θ)                      # one call, not two
    return _appleton_hartree_kernel(X, Y, Z, s, c)
end

"""
    appleton_hartree_rad(n_e, ν, B_vec, k, ω)

Return the squared index of refraction for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B_vec`: magnetic field as a 3-vector in T
- `k`: wave normal direction as a 3-vector (need not be normalized)
- `ω`: wave angular frequency in rad s⁻¹

# Returns
- Tuple `(O, X)` where `O` and `X` are the squared indices of refraction for the O and X modes, respectively

# Notes
- This is the ''full'' Appleton-Hartree equation, including collisions and valid for any angle θ.
- The squared indices of refraction are returned directly, as they are the natural output of the core calculation.
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Appleton1932]: Appleton, E. V. (1932), Wireless studies of the ionosphere, Int. J. Electr. Eng., 71, 642–650.
- [Budden1985]: Budden, K. G. (1985), The propagation of radio waves, Cambridge Univ. Press, Cambridge.
- [Haselgrove1960]: Haselgrove, C. B., and J. Haselgrove (1960), Twisted ray paths in the ionosphere, Proc. Phys. Soc. London, 75, 357–363.
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`appleton_hartree_hz`](@ref): equivalent function for Hz input
"""
function appleton_hartree_rad(n_e, ν, B_vec::AbstractVector, k::AbstractVector, ω::Real)
    B, s, c = field_magnitude_and_dircos(B_vec, k)
    X, Y, Z = _magnetoionic_parameters(n_e, ν, B, ω)
    return _appleton_hartree_kernel(X, Y, Z, s, c)
end

"""
    appleton_hartree_rad(n_e, ν, B_vec, ω)

Return the squared index of refraction for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B_vec`: magnetic field as a 3-vector in T
- `ω`: wave angular frequency in rad s⁻¹

# Returns
- Tuple `(O, X)` where `O` and `X` are the squared indices of refraction for the O and X modes, respectively

# Notes
- This assumes a vertical wave normal (parallel to the local vertical), so the angle θ is determined solely by the orientation of the magnetic field vector.
- This is the ''full'' Appleton-Hartree equation, including collisions.
- The squared indices of refraction are returned directly, as they are the natural output of the core calculation.
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Appleton1932]: Appleton, E. V. (1932), Wireless studies of the ionosphere, Int. J. Electr. Eng., 71, 642–650.
- [Budden1985]: Budden, K. G. (1985), The propagation of radio waves, Cambridge Univ. Press, Cambridge.
- [Haselgrove1960]: Haselgrove, C. B., and J. Haselgrove (1960), Twisted ray paths in the ionosphere, Proc. Phys. Soc. London, 75, 357–363.
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`appleton_hartree_hz`](@ref): equivalent function for Hz input
"""
function appleton_hartree_rad(n_e, ν, B_vec::AbstractVector, ω::Real)
    B, s, c = field_magnitude_and_dircos(B_vec)
    X, Y, Z = _magnetoionic_parameters(n_e, ν, B, ω)
    return _appleton_hartree_kernel(X, Y, Z, s, c)
end

"""
    appleton_hartree_hz(n_e, ν, B, θ, f)

Return the squared index of refraction for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B`: magnetic field strength in T
- `θ`: angle between the wave normal and the magnetic field in radians
- `f`: wave frequency in Hz

# Returns
- Tuple `(O, X)` where `O` and `X` are the squared indices of refraction for the O and X modes, respectively

# Notes
- This is the ''full'' Appleton-Hartree equation, including collisions and valid for any angle θ.
- The squared indices of refraction are returned directly, as they are the natural output of the core calculation.
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Appleton1932]: Appleton, E. V. (1932), Wireless studies of the ionosphere, Int. J. Electr. Eng., 71, 642–650.
- [Budden1985]: Budden, K. G. (1985), The propagation of radio waves, Cambridge Univ. Press, Cambridge.
- [Haselgrove1960]: Haselgrove, C. B., and J. Haselgrove (1960), Twisted ray paths in the ionosphere, Proc. Phys. Soc. London, 75, 357–363.
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`appleton_hartree_rad`](@ref): equivalent function for rad s⁻¹ input
"""
function appleton_hartree_hz(n_e, ν, B::Real, θ::Real, f::Real)
    ω = 2π * f
    return appleton_hartree_rad(n_e, ν, B, θ, ω)
end

"""
    appleton_hartree_hz(n_e, ν, B_vec, k, f)

Return the squared index of refraction for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B_vec`: magnetic field as a 3-vector in T
- `k`: wave normal direction as a 3-vector (need not be normalized)
- `f`: wave frequency in Hz

# Returns
- Tuple `(O, X)` where `O` and `X` are the squared indices of refraction for the O and X modes, respectively

# Notes
- This is the ''full'' Appleton-Hartree equation, including collisions and valid for any angle θ.
- The squared indices of refraction are returned directly, as they are the natural output of the core calculation.
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Appleton1932]: Appleton, E. V. (1932), Wireless studies of the ionosphere, Int. J. Electr. Eng., 71, 642–650.
- [Budden1985]: Budden, K. G. (1985), The propagation of radio waves, Cambridge Univ. Press, Cambridge.
- [Haselgrove1960]: Haselgrove, C. B., and J. Haselgrove (1960), Twisted ray paths in the ionosphere, Proc. Phys. Soc. London, 75, 357–363.
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`appleton_hartree_rad`](@ref): equivalent function for rad s⁻¹ input
"""
function appleton_hartree_hz(n_e, ν, B_vec::AbstractVector, k::AbstractVector, f::Real)
    ω = 2π * f
    return appleton_hartree_rad(n_e, ν, B_vec, k, ω)
end

"""
    appleton_hartree_hz(n_e, ν, B_vec, f)

Return the squared index of refraction for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B_vec`: magnetic field as a 3-vector in T
- `f`: wave frequency in Hz

# Returns
- Tuple `(O, X)` where `O` and `X` are the squared indices of refraction for the O and X modes, respectively

# Notes
- This assumes a vertical wave normal (parallel to the local vertical), so the angle θ is determined solely by the orientation of the magnetic field vector.
- This is the ''full'' Appleton-Hartree equation, including collisions.
- The squared indices of refraction are returned directly, as they are the natural output of the core calculation.
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Appleton1932]: Appleton, E. V. (1932), Wireless studies of the ionosphere, Int. J. Electr. Eng., 71, 642–650.
- [Budden1985]: Budden, K. G. (1985), The propagation of radio waves, Cambridge Univ. Press, Cambridge.
- [Haselgrove1960]: Haselgrove, C. B., and J. Haselgrove (1960), Twisted ray paths in the ionosphere, Proc. Phys. Soc. London, 75, 357–363.
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.


# See Also
- [`appleton_hartree_rad`](@ref): equivalent function for rad s⁻¹ input
"""
function appleton_hartree_hz(n_e, ν, B_vec::AbstractVector, f::Real)
    ω = 2π * f
    return appleton_hartree_rad(n_e, ν, B_vec, ω)
end