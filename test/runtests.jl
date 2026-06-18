using MagnetoIonic
using Test
using LinearAlgebra

# Alias for reaching non-exported helpers and module constants.
const MI = MagnetoIonic

@testset "MagnetoIonic.jl" begin

    # ------------------------------------------------------------------
    # utils.jl — plasma frequency
    # ------------------------------------------------------------------
    @testset "Plasma frequency" begin
        n_e = 1.0e11   # ~F-region electron density, m⁻³

        # Well-known constant: fₚ ≈ 8.98 √(n_e) Hz  (n_e in m⁻³)
        @test plasma_frequency_hz(n_e) ≈ 8.98 * sqrt(n_e) rtol = 1e-3

        # The four forms must be mutually consistent.
        @test plasma_frequency_rad(n_e)   ≈ sqrt(plasma_frequency_rad_sq(n_e))
        @test plasma_frequency_hz(n_e)    ≈ plasma_frequency_rad(n_e) / (2π)
        @test plasma_frequency_hz_sq(n_e) ≈ plasma_frequency_rad_sq(n_e) / (4π^2)
        @test plasma_frequency_hz_sq(n_e) ≈ plasma_frequency_hz(n_e)^2

        # ωp² scales linearly with electron density.
        @test plasma_frequency_rad_sq(2n_e) ≈ 2 * plasma_frequency_rad_sq(n_e)

        # Zero density ⇒ zero frequency.
        @test plasma_frequency_rad_sq(0.0) == 0.0
        @test plasma_frequency_hz(0.0)     == 0.0

        # Negative density is unphysical.
        @test_throws DomainError plasma_frequency_rad_sq(-1.0)
        @test_throws DomainError plasma_frequency_rad(-1.0)
        @test_throws DomainError plasma_frequency_hz(-1.0)
        @test_throws DomainError plasma_frequency_hz_sq(-1.0)
    end

    # ------------------------------------------------------------------
    # utils.jl — gyrofrequency
    # ------------------------------------------------------------------
    @testset "Gyrofrequency" begin
        B = 5.0e-5   # ~Earth surface field, T

        # Well-known constant: f_H ≈ 2.799e10 · B Hz  (≈ 28 GHz/T)
        @test MI.gyrofrequency_hz(B) ≈ 2.799e10 * B rtol = 1e-3

        # Consistency between forms.
        @test MI.gyrofrequency_hz(B)     ≈ MI.gyrofrequency_rad(B) / (2π)
        @test MI.gyrofrequency_rad_sq(B) ≈ MI.gyrofrequency_rad(B)^2
        @test MI.gyrofrequency_hz_sq(B)  ≈ MI.gyrofrequency_rad_sq(B) / (4π^2)
        @test MI.gyrofrequency_hz_sq(B)  ≈ MI.gyrofrequency_hz(B)^2

        # Linear in |B|.
        @test MI.gyrofrequency_rad(2B) ≈ 2 * MI.gyrofrequency_rad(B)

        # Vector form depends only on the magnitude.
        Bvec = [3.0e-5, 0.0, 4.0e-5]   # |B| = 5e-5
        @test MI.gyrofrequency_rad(Bvec) ≈ MI.gyrofrequency_rad(5.0e-5)
        @test MI.gyrofrequency_rad(Bvec) ≈ MI.gyrofrequency_rad(norm(Bvec))

        # Zero / negative field.
        @test MI.gyrofrequency_rad(0.0) == 0.0
        @test_throws DomainError MI.gyrofrequency_rad(-1.0)
    end

    # ------------------------------------------------------------------
    # utils.jl — wave geometry helpers
    # ------------------------------------------------------------------
    @testset "Wave geometry helpers" begin
        ẑ = [0.0, 0.0, 1.0]
        x̂ = [1.0, 0.0, 0.0]

        # Parallel / perpendicular / antiparallel.
        @test MI.wave_normal_angle(ẑ, ẑ)  ≈ 0.0 atol = 1e-12
        @test MI.wave_normal_angle(x̂, ẑ)  ≈ π / 2
        @test MI.wave_normal_angle(ẑ, -ẑ) ≈ π
        @test MI.wave_normal_angle([1.0, 0.0, 1.0], ẑ) ≈ π / 4

        # Direction only — insensitive to vector length.
        @test MI.wave_normal_angle(2x̂, 3ẑ) ≈ π / 2

        # clamp guards against NaN for (anti)parallel inputs.
        @test !isnan(MI.wave_normal_angle(ẑ, ẑ))
        @test !isnan(MI.wave_normal_angle(ẑ, -ẑ))

        # field_magnitude_and_angle: |B| = 5e-5, 36.87° from vertical.
        Bvec = [0.0, 3.0e-5, 4.0e-5]
        B, θ = MI.field_magnitude_and_angle(Bvec, ẑ)
        @test B ≈ 5.0e-5
        @test θ ≈ acos(4.0e-5 / 5.0e-5)
        @test θ ≈ MI.wave_normal_angle(ẑ, Bvec)

        # Default wave normal is vertical.
        Bd, θd = MI.field_magnitude_and_angle(Bvec)
        @test Bd ≈ B
        @test θd ≈ θ

        # vertical_wave_normal_angle == angle between vertical and B.
        @test MI.vertical_wave_normal_angle(Bvec) ≈ MI.wave_normal_angle(ẑ, Bvec)

        # field_magnitude_and_dircos returns (|B|, sinθ, cosθ) consistent with θ.
        Bc, s, c = MI.field_magnitude_and_dircos(Bvec, ẑ)
        @test Bc ≈ B
        @test c ≈ cos(θ)
        @test s ≈ sin(θ)
        @test s ≥ 0                    # θ ∈ [0,π] ⇒ sinθ ≥ 0
        @test s^2 + c^2 ≈ 1

        # Default wave normal is vertical here too.
        Bc2, s2, c2 = MI.field_magnitude_and_dircos(Bvec)
        @test Bc2 ≈ Bc && s2 ≈ s && c2 ≈ c
    end

    # ------------------------------------------------------------------
    # appletonhartree.jl
    # ------------------------------------------------------------------
    @testset "Appleton-Hartree" begin
        n_e = 1.0e11
        f   = 5.0e6          # 5 MHz (above fₚ ≈ 2.84 MHz ⇒ X < 1)
        ω   = 2π * f
        B   = 5.0e-5
        ν   = 1.0e5

        X = plasma_frequency_rad_sq(n_e) / ω^2
        Y = MI.gyrofrequency_rad(B) / ω

        @testset "rad / Hz agreement" begin
            for θ in (0.0, π / 6, π / 3, π / 2)
                r = appleton_hartree_rad(n_e, ν, B, θ, ω)
                h = appleton_hartree_hz(n_e, ν, B, θ, f)
                @test r.O ≈ h.O
                @test r.X ≈ h.X
            end
        end

        @testset "Unmagnetized, collisionless ⇒ n² = 1 - X" begin
            for θ in (0.0, π / 4, π / 2)
                res = appleton_hartree_rad(n_e, 0.0, 0.0, θ, ω)
                @test res.O ≈ 1 - X
                @test res.X ≈ 1 - X
            end
        end

        @testset "Longitudinal (θ=0), collisionless" begin
            res = appleton_hartree_rad(n_e, 0.0, B, 0.0, ω)
            @test res.O ≈ 1 - X / (1 + Y)
            @test res.X ≈ 1 - X / (1 - Y)
        end

        @testset "Transverse (θ=π/2), collisionless" begin
            res = appleton_hartree_rad(n_e, 0.0, B, π / 2, ω)
            # Ordinary mode is unaffected by the field.
            @test res.O ≈ 1 - X
            # Extraordinary mode.
            @test res.X ≈ 1 - X * (1 - X) / (1 - X - Y^2)
        end

        @testset "Vector form matches scalar form" begin
            Bvec = [0.0, 3.0e-5, 4.0e-5]   # |B| = 5e-5
            k    = [0.0, 0.0, 1.0]
            Bmag = norm(Bvec)
            θ    = MI.wave_normal_angle(k, Bvec)

            res_vec = appleton_hartree_rad(n_e, ν, Bvec, k, ω)
            res_sca = appleton_hartree_rad(n_e, ν, Bmag, θ, ω)
            @test res_vec.O ≈ res_sca.O
            @test res_vec.X ≈ res_sca.X

            # Vertical-wave-normal convenience form (k assumed vertical).
            res_vert = appleton_hartree_rad(n_e, ν, Bvec, ω)
            @test res_vert.O ≈ res_vec.O
            @test res_vert.X ≈ res_vec.X

            # Hz vector forms.
            @test appleton_hartree_hz(n_e, ν, Bvec, k, f).O ≈ res_vec.O
            @test appleton_hartree_hz(n_e, ν, Bvec, f).X    ≈ res_vert.X
        end

        @testset "Internal helpers" begin
            Xp, Yp, Zp = MI._magnetoionic_parameters(n_e, ν, B, ω)
            @test Xp ≈ X
            @test Yp ≈ Y
            @test Zp ≈ ν / ω

            # Kernel directly: collisionless transverse ⇒ O mode = 1 - X.
            res = MI._appleton_hartree_kernel(X, Y, 0.0, 1.0, 0.0)  # sinθ=1, cosθ=0
            @test res.O ≈ 1 - X
        end
    end

    # ------------------------------------------------------------------
    # radioabsorption.jl
    # ------------------------------------------------------------------
    @testset "Radio absorption" begin
        n_e = 1.0e11
        f   = 5.0e6
        ω   = 2π * f
        B   = 5.0e-5
        ν   = 1.0e6

        @testset "NEPER_TO_DB constant" begin
            @test MI.NEPER_TO_DB ≈ 10 * log10(ℯ)
        end

        @testset "rad / Hz agreement" begin
            for θ in (0.0, π / 4, π / 2)
                r = radio_absorption_rad(n_e, ν, B, θ, ω)
                h = radio_absorption_hz(n_e, ν, B, θ, f)
                @test r.O ≈ h.O
                @test r.X ≈ h.X
            end
        end

        @testset "Matches the Appleton-Hartree definition" begin
            θ   = π / 4
            nsq = appleton_hartree_rad(n_e, ν, B, θ, ω)
            αO  = MI.NEPER_TO_DB * 2 * imag(sqrt(nsq.O)) * ω / MI.C_LIGHT
            αX  = MI.NEPER_TO_DB * 2 * imag(sqrt(nsq.X)) * ω / MI.C_LIGHT
            res = radio_absorption_rad(n_e, ν, B, θ, ω)
            @test res.O ≈ αO
            @test res.X ≈ αX
        end

        @testset "Collisions ⇒ positive absorption" begin
            res = radio_absorption_rad(n_e, ν, B, π / 4, ω)
            @test res.O > 0
            @test res.X > 0
        end

        @testset "Lossless propagating wave ⇒ zero absorption" begin
            # B=0, ν=0, X<1 ⇒ n² = 1-X is real & positive ⇒ no imaginary part.
            res = radio_absorption_rad(n_e, 0.0, 0.0, π / 4, ω)
            @test res.O ≈ 0.0 atol = 1e-9
            @test res.X ≈ 0.0 atol = 1e-9
        end

        @testset "Vector and scalar forms agree" begin
            Bvec = [0.0, 3.0e-5, 4.0e-5]
            k    = [0.0, 0.0, 1.0]
            Bmag = norm(Bvec)
            θ    = MI.wave_normal_angle(k, Bvec)

            rv = radio_absorption_rad(n_e, ν, Bvec, k, ω)
            rs = radio_absorption_rad(n_e, ν, Bmag, θ, ω)
            @test rv.O ≈ rs.O
            @test rv.X ≈ rs.X

            rvert = radio_absorption_rad(n_e, ν, Bvec, ω)
            @test rvert.O ≈ rv.O

            @test radio_absorption_hz(n_e, ν, Bvec, k, f).O ≈ rv.O
            @test radio_absorption_hz(n_e, ν, Bvec, f).X    ≈ rvert.X
        end
    end

end