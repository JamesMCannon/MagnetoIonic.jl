"""
    MagnetoIonic

Tools for magnetoionic wave propagation in a cold, magnetized, collisional
plasma such as the ionosphere.

The package computes the complex refractive index of the ordinary (O) and
extraordinary (X) modes from the Appleton-Hartree equation (returned as n²), the characteristic
plasma and electron gyrofrequencies, and the resulting radio wave absorption.
Inputs are accepted either as scalars (field magnitude and wave normal angle)
or as vectors (geomagnetic field and wave normal direction), and in either
angular frequency (rad s⁻¹) or ordinary frequency (Hz).

# Exported functions
- [`plasma_frequency_rad`](@ref), [`plasma_frequency_hz`](@ref), [`plasma_frequency_rad_sq`](@ref), [`plasma_frequency_hz_sq`](@ref)
- [`gyrofrequency_rad`](@ref), [`gyrofrequency_hz`](@ref), [`gyrofrequency_rad_sq`](@ref), [`gyrofrequency_hz_sq`](@ref)
- [`appleton_hartree_rad`](@ref), [`appleton_hartree_hz`](@ref)
- [`radio_absorption_rad`](@ref), [`radio_absorption_hz`](@ref)
"""
module MagnetoIonic

using LinearAlgebra
import PhysicalConstants.CODATA2022: c_0, ElementaryCharge, m_e, ε_0, ustrip

const C_LIGHT  = Float64(ustrip(c_0))              # speed of light (m/s)
const Q_ELECTRON = Float64(ustrip(ElementaryCharge)) # elementary charge (C)
const M_ELECTRON = Float64(ustrip(m_e))            # electron mass (kg)
const EPSILON_0  = Float64(ustrip(ε_0))            # vacuum permittivity (F/m)
const NEPER_TO_DB = 10 * log10(ℯ)                   # conversion factor from Np to dB

include("utils.jl")
include("appletonhartree.jl")
include("radioabsorption.jl")

export plasma_frequency_rad, plasma_frequency_hz, plasma_frequency_rad_sq, plasma_frequency_hz_sq
export gyrofrequency_rad, gyrofrequency_hz, gyrofrequency_rad_sq, gyrofrequency_hz_sq
export appleton_hartree_rad, appleton_hartree_hz
export radio_absorption_rad, radio_absorption_hz

end
