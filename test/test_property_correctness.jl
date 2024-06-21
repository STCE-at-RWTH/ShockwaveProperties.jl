using LinearAlgebra
using ShockwaveProperties
using Test
using Unitful

@testset "Dimensional Analysis" begin
    # giving different units shouldn't mess with the actual results
    # these are the same state
    s1 = PrimitiveProps(1.225, [2.0, 0.0], 300.0)
    s2 = PrimitiveProps(0.002376892407u"slug/ft^3", [2.0, 0.0], 540.0u"Ra")
    u1 = ConservedProps(s1, DRY_AIR)
    u2 = ConservedProps(s2, DRY_AIR)
    #check dimensions
    for v ∈ (s1, s2, u1, u2)
        @test pressure(v, DRY_AIR) isa Unitful.Pressure
        @test temperature(v, DRY_AIR) isa Unitful.Temperature
        @test (
            specific_static_internal_energy(v, DRY_AIR) isa
            ShockwaveProperties.SpecificEnergy
        )
        @test (
            total_internal_energy_density(v, DRY_AIR) isa
            ShockwaveProperties.EnergyDensity
        )
        @test speed_of_sound(v, DRY_AIR) isa Unitful.Velocity
    end
    # check equivalency between primitive/conserved
    # cross test with unit change as well
    test_items = ((s1, u1), (s2, u2), (s1, u2), (s2, u1))
    @testset "Pairwise Equivalency $i" for i ∈ eachindex(test_items)
        (v1, v2) = test_items[i]
        @test pressure(v1, DRY_AIR) ≈ pressure(v2, DRY_AIR)
        @test temperature(v1, DRY_AIR) ≈ temperature(v2, DRY_AIR)
        @test all(
            momentum_density(v1, DRY_AIR) .≈ momentum_density(v2, DRY_AIR),
        )
        @test all(velocity(v1, DRY_AIR) .≈ velocity(v2, DRY_AIR))
        @test (
            specific_static_internal_energy(v1, DRY_AIR) ≈
            specific_static_internal_energy(v2, DRY_AIR)
        )
        @test speed_of_sound(v1, DRY_AIR) ≈ speed_of_sound(v2, DRY_AIR)
    end
end

@testset "Convert Primitve ↔ Conserved" begin
    s1 = PrimitiveProps(1.225, [2.0, 0.0], 300.0)
    u = ConservedProps(s1, DRY_AIR)
    s2 = PrimitiveProps(u, DRY_AIR)
    @test s1.ρ ≈ s2.ρ
    @test all(s1.M .≈ s2.M)
    @test s1.T ≈ s2.T
end

@testset "Equivalency Between Unitful and Unitless" begin
    n = [-1.0, 0]
    t = [0.0, 1.0]
    sL = PrimitiveProps(1.225, [2.0, 0.0], 300.0)
    sR = state_behind(sL, n, t, DRY_AIR)
    sL_nounits = state_to_vector(sL)
    sR_nounits = primitive_state_behind(sL_nounits, n, t, DRY_AIR)
    @testset "Primitive State Quantities" begin
        @test sR_nounits[1] ≈ ustrip(sR.ρ)
        @test all(sR_nounits[2:end-1] .≈ sR.M)
        @test sR_nounits[end] ≈ ustrip(sR.T)
        @test pressure(sL_nounits[1], sL_nounits[end], DRY_AIR) ≈
              pressure(sL, DRY_AIR)
        @test pressure(sR_nounits[1], sR_nounits[end], DRY_AIR) ≈
              pressure(sR, DRY_AIR)
    end

    uL = ConservedProps(sL, DRY_AIR)
    uL_nounits = state_to_vector(uL)
    uR = state_behind(uL, n, t, DRY_AIR)
    uR_nounits = conserved_state_behind(uL_nounits, n, t, DRY_AIR)

    @testset "Conserved State Quantities" begin
        @test uR_nounits[1] ≈ ustrip(uR.ρ)
        @test all(uR_nounits[2:end-1] .≈ ustrip.(uR.ρv))
        @test uR_nounits[end] ≈ ustrip(uR.ρE)
        ρe_nounits = static_internal_energy_density(
            uR_nounits[1],
            uR_nounits[2:end-1],
            uR_nounits[end],
        )
        @test ρe_nounits ≈ ustrip(static_internal_energy_density(uR))
        @test pressure(ρe_nounits, DRY_AIR) ≈ pressure(uR, DRY_AIR)
    end

    @testset "Conversion" begin
        uR_nounits2 = conserved_state_vector(sR_nounits, DRY_AIR)
        sR_nounits2 = primitive_state_vector(uR_nounits, DRY_AIR)
        @test uR_nounits[1] ≈ uR_nounits2[1]
        @test all(uR_nounits .≈ uR_nounits2)
        @test all(sR_nounits .≈ sR_nounits2)
    end
end

# flux for the euler equations
function F(u::ConservedProps)
    v = u.ρv / u.ρ # velocity
    P = pressure(u, DRY_AIR)
    return vcat(u.ρv', (u.ρv .* v' + I * P), (v .* (u.ρE + P))') # stack row vectors
end

@testset "Rankine-Hugoniot Condition" begin
    free_stream = PrimitiveProps(1.225, [2.0, 0.0], 300.0)
    u_L = ConservedProps(free_stream, DRY_AIR)
    # restricted domain here because formula from Anderson&Anderson
    # breaks down near β = 0
    # TODO investigate this.
    for θ ∈ range(0, π / 3; length = 40)
        n = [-cos(θ), sin(θ)]
        t = [0 1; -1 0] * n
        u_R = state_behind(u_L, n, t, DRY_AIR)
        # F(u_l)⋅n̂ - F(u_r)⋅n̂ = 0 ⟹ F(u_l)⋅n̂ = F(u_r)⋅n̂ 
        @test all(F(u_L) * n .≈ F(u_R) * n)
    end
end

@testset "Speed of Sound" begin
    gas = DRY_AIR
    s = PrimitiveProps(1.225, [2.0, 0.0], 300.0)
    u = ConservedProps(s, gas)

    a_s = speed_of_sound(s, gas)
    a_u = speed_of_sound(u, gas)

    a_s_pressure = speed_of_sound(density(s), pressure(s, gas), gas)
    a_u_pressure = speed_of_sound(density(u), pressure(u, gas), gas)

    a_s_ie = speed_of_sound(density(s), momentum_density(s, gas), total_internal_energy_density(s, gas), gas)
    a_u_ie = speed_of_sound(density(u), momentum_density(u), total_internal_energy_density(u), gas)

    @test a_s ≈ a_u
    @test a_s_pressure ≈ a_u_pressure
    @test a_s_ie ≈ a_u_ie
end