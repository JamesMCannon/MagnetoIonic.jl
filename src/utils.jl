"""
    plasma_frequency_rad_sq(n_e)

Return the squared plasma angular frequency ωp² in rad² s⁻².

# Arguments
- `n_e`: electron number density in m⁻³

# Returns
- ωp² in rad² s⁻²

# Notes
The plasma frequency is defined as the natural oscillation
frequency of the electron gas in the absence of a magnetic field:
`ωp² = n_e * e² / (m_e * ε₀)`
"""
function plasma_frequency_rad_sq(n_e::Real)
    n_e < 0 && throw(DomainError(n_e, "electron density must be non-negative"))
    return (n_e * Q_ELECTRON^2) / (M_ELECTRON * EPSILON_0)
end

"""
    plasma_frequency_rad(n_e)

Return the plasma angular frequency ωp in rad s⁻¹.

# Arguments
- `n_e`: electron number density in m⁻³

# Returns
- ωp in rad s⁻¹

# See Also
- [`plasma_frequency_rad_sq`](@ref): the squared form used internally
- [`plasma_frequency_hz`](@ref): equivalent in Hz
"""
function plasma_frequency_rad(n_e)
    return sqrt(plasma_frequency_rad_sq(n_e))
end

"""
    plasma_frequency_hz_sq(n_e)

Return the squared plasma frequency fp² in Hz².

Computed as ωp² / (4π²) to defer the Hz conversion until after
the core rad² s⁻² calculation, minimizing floating point operations.

# Arguments
- `n_e`: electron number density in m⁻³

# Returns
- fp² in Hz²

# See Also
- [`plasma_frequency_rad_sq`](@ref): the primary internal computation
- [`plasma_frequency_hz`](@ref): the unsquared form in Hz
"""
function plasma_frequency_hz_sq(n_e)
    return plasma_frequency_rad_sq(n_e) / (4π^2)
end

"""
    plasma_frequency_hz(n_e)

Return the plasma frequency fp in Hz.

Computed as sqrt(ωp²) / (2π) to defer the Hz conversion until after
the square root.

# Arguments
- `n_e`: electron number density in m⁻³

# Returns
- fp in Hz

# Notes
The plasma frequency represents the natural oscillation frequency of
the electron gas and is a key parameter in the Appleton-Hartree equation.
For the D-region ionosphere, typical values range from a few kHz at night
to several MHz during daytime solar maximum conditions.

# See Also
- [`plasma_frequency_rad`](@ref): equivalent in rad s⁻¹
- [`plasma_frequency_hz_sq`](@ref): the squared form in Hz²
"""
function plasma_frequency_hz(n_e)
    return sqrt(plasma_frequency_rad_sq(n_e)) / (2π)
end

"""
    gyrofrequency_rad(B)

Return the electron gyrofrequency ωH in rad s⁻¹.

# Arguments
- `B`: geomagnetic field magnitude in Tesla as a scalar
"""
function gyrofrequency_rad(B::Real)
    B < 0 && throw(DomainError(B, "magnetic field magnitude must be non-negative"))
    return abs(Q_ELECTRON) * B / M_ELECTRON
end

"""
    gyrofrequency_rad(B_vec)

Return the electron gyrofrequency ωH in rad s⁻¹.

Computes the gyrofrequency from the magnitude of a geomagnetic
field vector.

# Arguments
- `B_vec`: geomagnetic field vector in Tesla
"""
function gyrofrequency_rad(B_vec::AbstractVector)
    return gyrofrequency_rad(norm(B_vec))
end

"""
    gyrofrequency_hz(B)

Return the electron gyrofrequency fH in Hz.

# Arguments
- `B`: geomagnetic field magnitude in Tesla, or field vector in Tesla

# See Also
- [`gyrofrequency_rad`](@ref): equivalent in rad s⁻¹
"""
gyrofrequency_hz(B) = gyrofrequency_rad(B) / (2π)

"""
    gyrofrequency_rad_sq(B)

Return the squared electron gyrofrequency ωH² in rad² s⁻².

# Arguments
- `B`: geomagnetic field magnitude in Tesla, or field vector in Tesla

# See Also
- [`gyrofrequency_rad`](@ref): base calculation in rad s⁻¹
- [`gyrofrequency_hz_sq`](@ref): equivalent in Hz²
"""
gyrofrequency_rad_sq(B) = gyrofrequency_rad(B)^2

"""
    gyrofrequency_hz_sq(B)

Return the squared electron gyrofrequency fH² in Hz².

# Arguments
- `B`: geomagnetic field magnitude in Tesla, or field vector in Tesla

# See Also
- [`gyrofrequency_rad`](@ref): base calculation in rad s⁻¹
- [`gyrofrequency_rad_sq`](@ref): equivalent in rad² s⁻²
"""
gyrofrequency_hz_sq(B) = gyrofrequency_rad_sq(B) / (4π^2)

"""
    wave_normal_angle(k, B_vec)

Return the angle θ in radians between wave normal vector `k`
and geomagnetic field vector `B`.

# Arguments
- `k`: wave normal direction as a 3-vector (need not be normalized)
- `B_vec`: geomagnetic field as a 3-vector in Tesla

# Notes
This is clamped at ±1 before applying `acos` to avoid NaN results from 
floating point imprecision when `k` and `B_vec` are nearly parallel or antiparallel.
"""
function wave_normal_angle(k::AbstractVector, B_vec::AbstractVector)
    cos_θ = dot(k, B_vec) / (norm(k) * norm(B_vec))
    return acos(clamp(cos_θ, -1.0, 1.0))
end

"""
    field_magnitude_and_angle(B_vec, k)

Return a tuple (|B|, θ) given a geomagnetic field vector and wave normal.

# Arguments
- `B_vec`: geomagnetic field as a 3-vector in Tesla
- `k`: wave normal direction as a 3-vector (need not be normalized)

# Returns
- Tuple `(B, θ)` where `B` is magnitude in Tesla and `θ` is in radians
"""
function field_magnitude_and_angle(B_vec::AbstractVector, k::AbstractVector)
    B = norm(B_vec)
    θ = wave_normal_angle(k, B_vec)
    return B, θ
end

const VERTICAL = [0.0, 0.0, 1.0]

"""
    field_magnitude_and_angle(B_vec)

Return a tuple (|B|, θ) given a geomagnetic field vector and assumed vertical wave normal.

# Arguments
- `B_vec`: geomagnetic field as a 3-vector in Tesla

# Returns
- Tuple `(B, θ)` where `B` is magnitude in Tesla and `θ` is in radians
"""
function field_magnitude_and_angle(B_vec::AbstractVector)
    return field_magnitude_and_angle(B_vec, VERTICAL)
end

function field_magnitude_and_dircos(B_vec::AbstractVector, k::AbstractVector)
    B = norm(B_vec)
    c = clamp(dot(k, B_vec) / (norm(k) * B), -1.0, 1.0)
    s = sqrt(1 - c * c)        # θ ∈ [0,π] ⇒ sinθ ≥ 0
    return B, s, c
end

function field_magnitude_and_dircos(B_vec::AbstractVector)
    return field_magnitude_and_dircos(B_vec, VERTICAL)
end

"""
    vertical_wave_normal_angle(B_vec)

Return θ for vertical incidence propagation given a geomagnetic
field vector, i.e. the angle between vertical and B.

Equivalent to the magnetic dip angle measured from vertical
rather than from horizontal.

# Arguments  
- `B_vec`: geomagnetic field as a 3-vector in any consistent
  coordinate system where the third component is vertical
"""
function vertical_wave_normal_angle(B_vec::AbstractVector)
    return wave_normal_angle(VERTICAL, B_vec)
end