module MagnetoIonic

using LinearAlgebra
import PhysicalConstants.CODATA2022: c_0, ElementaryCharge, m_e, ε_0, ustrip

const C_LIGHT  = Float64(ustrip(c_0))              # speed of light (m/s)
const Q_ELECTRON = Float64(ustrip(ElementaryCharge)) # elementary charge (C)
const M_ELECTRON = Float64(ustrip(m_e))            # electron mass (kg)
const EPSILON_0  = Float64(ustrip(ε_0))            # vacuum permittivity (F/m)

include("utils.jl")
include("appletonhartree.jl")
include("radioabsorption.jl")

# Write your package code here.
export plasma_frequency_rad, plasma_frequency_hz, plasma_frequency_rad_sq, plasma_frequency_hz_sq
export appleton_hartree_rad, appleton_hartree_hz
export radio_absorption_rad, radio_absorption_hz

end
