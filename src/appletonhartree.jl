@inline function _appleton_hartree_kernel(X, Y, Z, sinθ, cosθ)
    YL  = Y * cosθ
    YT  = Y * sinθ
    YT² = YT * YT

    U  = complex(one(Z), Z)              # 1 + im*Z, always ComplexF64
    UX = U - X
    Δ  = sqrt(YT² * YT² / 4 + YL^2 * UX^2)   # already complex since UX is

    term = U - YT² / (2 * UX)            # computed once
    disc = Δ / UX                        # computed once
    n²_O = 1 - X / (term + disc)
    n²_X = 1 - X / (term - disc)
    return (O = n²_O, X = n²_X)
end


function _magnetoionic_parameters(Ne, ν, B, ω)
    X = plasma_frequency_rad_sq(Ne) / ω^2   # no sqrt ever taken
    Y = gyrofrequency_rad(B) / ω
    Z = ν / ω
    return X, Y, Z
end

"""
    appleton_hartree_rad(n_e, ν, B, θ, ω)

Return the squared index of refraction for the O and X modes.

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

# See Also
- [`appleton_hartree_rad`](@ref): equivalent function for rad s⁻¹ input
"""
function appleton_hartree_hz(n_e, ν, B_vec::AbstractVector, k::AbstractVector, f)
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

# See Also
- [`appleton_hartree_rad`](@ref): equivalent function for rad s⁻¹ input
"""
function appleton_hartree_hz(n_e, ν, B_vec::AbstractVector, f)
    ω = 2π * f
    return appleton_hartree_rad(n_e, ν, B_vec, ω)
end