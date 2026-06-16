
function radio_absorption_rad(n_e, ν, B_vec::AbstractVector, k::AbstractVector, ω::Real)
    n²_O,n²_X = appleton_hartree_rad(n_e, ν, B_vec, k, ω)

    αO = 2 * imag(sqrt(n²_O)) * ω / C_LIGHT
    αX = 2 * imag(sqrt(n²_X)) * ω / C_LIGHT 

    Absorp_O = NEPER_TO_DB * αO # dB/m
    Absorp_X = NEPER_TO_DB * αX # dB/m
    return (O=Absorp_O, X=Absorp_X)
end

function radio_absorption_rad(n_e, ν, B_vec::AbstractVector, ω::Real)
    n²_O,n²_X = appleton_hartree_rad(n_e, ν, B_vec, ω)

    αO = 2 * imag(sqrt(n²_O)) * ω / C_LIGHT
    αX = 2 * imag(sqrt(n²_X)) * ω / C_LIGHT 

    Absorp_O = NEPER_TO_DB * αO # dB/m
    Absorp_X = NEPER_TO_DB * αX # dB/m
    return (O=Absorp_O, X=Absorp_X)
end

function radio_absorption_rad(n_e, ν, B::Real, θ::Real, ω::Real)
    n²_O,n²_X = appleton_hartree_rad(n_e, ν, B, θ, ω)

    αO = 2 * imag(sqrt(n²_O)) * ω / C_LIGHT
    αX = 2 * imag(sqrt(n²_X)) * ω / C_LIGHT 

    Absorp_O = NEPER_TO_DB * αO # dB/m
    Absorp_X = NEPER_TO_DB * αX # dB/m
    return (O=Absorp_O, X=Absorp_X)
end

function radio_absorption_hz(n_e, ν, B_vec::AbstractVector, k::AbstractVector, f::Real)
    ω = 2π * f
    return radio_absorption_rad(n_e, ν, B_vec, k, ω)
end

function radio_absorption_hz(n_e, ν, B_vec::AbstractVector, f::Real)
    ω = 2π * f
    return radio_absorption_rad(n_e, ν, B_vec, ω)
end

function radio_absorption_hz(n_e, ν, B::Real, θ::Real, f::Real)
    ω = 2π * f
    return radio_absorption_rad(n_e, ν, B, θ, ω)
end