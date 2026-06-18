"""
    radio_absorption_rad(n_e, ν, B_vec, k, ω)

Return the power absorption coefficient for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B_vec`: magnetic field as a 3-vector in T
- `k`: wave normal direction as a 3-vector (need not be normalized)
- `ω`: wave angular frequency in rad s⁻¹

# Returns
- NamedTuple `(O, X)` where `O` and `X` are the absorption coefficients in dB m⁻¹ for the O and X modes, respectively

# Notes
- The absorption coefficient describes the spatial decay of wave power, `P(z) = P₀ · 10^(-α z / 10)` with `α` in dB m⁻¹.
- It is obtained from the imaginary part of the refractive index `n = √(n²)`, where the squared index `n²` is returned by the Appleton-Hartree equation, as `α = 2 · imag(n) · ω / c`, a power absorption coefficient in Np m⁻¹, then converted to dB m⁻¹ via the factor `10 · log₁₀(e)`.
- Absorption is non-zero when the medium is lossy (collisions present) or the wave is evanescent.
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`appleton_hartree_rad`](@ref): the refractive index this is derived from, including references
- [`radio_absorption_hz`](@ref): equivalent function for Hz input
"""
function radio_absorption_rad(n_e, ν, B_vec::AbstractVector, k::AbstractVector, ω::Real)
    n²_O,n²_X = appleton_hartree_rad(n_e, ν, B_vec, k, ω)

    αO = 2 * imag(sqrt(n²_O)) * ω / C_LIGHT
    αX = 2 * imag(sqrt(n²_X)) * ω / C_LIGHT 

    Absorp_O = NEPER_TO_DB * αO # dB/m
    Absorp_X = NEPER_TO_DB * αX # dB/m
    return (O=Absorp_O, X=Absorp_X)
end

"""
    radio_absorption_rad(n_e, ν, B_vec, ω)

Return the power absorption coefficient for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B_vec`: magnetic field as a 3-vector in T
- `ω`: wave angular frequency in rad s⁻¹

# Returns
- NamedTuple `(O, X)` where `O` and `X` are the absorption coefficients in dB m⁻¹ for the O and X modes, respectively

# Notes
- This assumes a vertical wave normal (parallel to the local vertical), so the angle θ is determined solely by the orientation of the magnetic field vector.
- See the five-argument method for the definition of the absorption coefficient.
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`appleton_hartree_rad`](@ref): the refractive index this is derived from, including references
- [`radio_absorption_hz`](@ref): equivalent function for Hz input
"""
function radio_absorption_rad(n_e, ν, B_vec::AbstractVector, ω::Real)
    n²_O,n²_X = appleton_hartree_rad(n_e, ν, B_vec, ω)

    αO = 2 * imag(sqrt(n²_O)) * ω / C_LIGHT
    αX = 2 * imag(sqrt(n²_X)) * ω / C_LIGHT 

    Absorp_O = NEPER_TO_DB * αO # dB/m
    Absorp_X = NEPER_TO_DB * αX # dB/m
    return (O=Absorp_O, X=Absorp_X)
end

"""
    radio_absorption_rad(n_e, ν, B, θ, ω)

Return the power absorption coefficient for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B`: magnetic field strength in T
- `θ`: angle between the wave normal and the magnetic field in radians
- `ω`: wave angular frequency in rad s⁻¹

# Returns
- NamedTuple `(O, X)` where `O` and `X` are the absorption coefficients in dB m⁻¹ for the O and X modes, respectively

# Notes
- See the `B_vec, k` method for the definition of the absorption coefficient.
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`appleton_hartree_rad`](@ref): the refractive index this is derived from, including references
- [`radio_absorption_hz`](@ref): equivalent function for Hz input
"""
function radio_absorption_rad(n_e, ν, B::Real, θ::Real, ω::Real)
    n²_O,n²_X = appleton_hartree_rad(n_e, ν, B, θ, ω)

    αO = 2 * imag(sqrt(n²_O)) * ω / C_LIGHT
    αX = 2 * imag(sqrt(n²_X)) * ω / C_LIGHT 

    Absorp_O = NEPER_TO_DB * αO # dB/m
    Absorp_X = NEPER_TO_DB * αX # dB/m
    return (O=Absorp_O, X=Absorp_X)
end

"""
    radio_absorption_hz(n_e, ν, B_vec, k, f)

Return the power absorption coefficient for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B_vec`: magnetic field as a 3-vector in T
- `k`: wave normal direction as a 3-vector (need not be normalized)
- `f`: wave frequency in Hz

# Returns
- NamedTuple `(O, X)` where `O` and `X` are the absorption coefficients in dB m⁻¹ for the O and X modes, respectively

# Notes
- Converts the frequency to angular frequency (`ω = 2πf`) and forwards to [`radio_absorption_rad`](@ref).
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`radio_absorption_rad`](@ref): equivalent function for rad s⁻¹ input
- [`appleton_hartree_hz`](@ref): the refractive index this is derived from, including references
"""
function radio_absorption_hz(n_e, ν, B_vec::AbstractVector, k::AbstractVector, f::Real)
    ω = 2π * f
    return radio_absorption_rad(n_e, ν, B_vec, k, ω)
end

"""
    radio_absorption_hz(n_e, ν, B_vec, f)

Return the power absorption coefficient for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B_vec`: magnetic field as a 3-vector in T
- `f`: wave frequency in Hz

# Returns
- NamedTuple `(O, X)` where `O` and `X` are the absorption coefficients in dB m⁻¹ for the O and X modes, respectively

# Notes
- This assumes a vertical wave normal (parallel to the local vertical), so the angle θ is determined solely by the orientation of the magnetic field vector.
- Converts the frequency to angular frequency (`ω = 2πf`) and forwards to [`radio_absorption_rad`](@ref).
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`radio_absorption_rad`](@ref): equivalent function for rad s⁻¹ input
- [`appleton_hartree_hz`](@ref): the refractive index this is derived from, including references
"""
function radio_absorption_hz(n_e, ν, B_vec::AbstractVector, f::Real)
    ω = 2π * f
    return radio_absorption_rad(n_e, ν, B_vec, ω)
end

"""
    radio_absorption_hz(n_e, ν, B, θ, f)

Return the power absorption coefficient for the O and X modes.

# Arguments
- `n_e`: electron number density in m⁻³
- `ν`: effective collision frequency in s⁻¹
- `B`: magnetic field strength in T
- `θ`: angle between the wave normal and the magnetic field in radians
- `f`: wave frequency in Hz

# Returns
- NamedTuple `(O, X)` where `O` and `X` are the absorption coefficients in dB m⁻¹ for the O and X modes, respectively

# Notes
- Converts the frequency to angular frequency (`ω = 2πf`) and forwards to [`radio_absorption_rad`](@ref).
- See [Zawdie2017] for a detailed discussion on the choice of collision frequency.

# References
- [Zawdie2017]: Zawdie, K. A., D. P. Drob, D. E. Siskind, and C. Coker (2017), Calculating the absorption of HF radio waves in the ionosphere, Radio Sci., 52, 767–783, doi:10.1002/2017RS006256.

# See Also
- [`radio_absorption_rad`](@ref): equivalent function for rad s⁻¹ input
- [`appleton_hartree_hz`](@ref): the refractive index this is derived from, including references
"""
function radio_absorption_hz(n_e, ν, B::Real, θ::Real, f::Real)
    ω = 2π * f
    return radio_absorption_rad(n_e, ν, B, θ, ω)
end