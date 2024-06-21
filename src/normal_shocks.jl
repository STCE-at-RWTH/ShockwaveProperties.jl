"""
**Equation 4.8** from Anderson&Anderson

Computes the density change across a _stationary_ shock wave. 
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
function shock_density_ratio(M, n̂, gas::CaloricallyPerfectGas)
    Mn2 = (M ⋅ n̂)^2
    return ((gas.γ + 1) * Mn2) / ((gas.γ - 1) * Mn2 + 2)
end

"""
**Equation 4.9** from Anderson&Anderson

Computes the pressure ratio across a _stationary_ shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
function shock_pressure_ratio(M, n̂, gas::CaloricallyPerfectGas)
    Mn2 = (M ⋅ n̂)^2
    return 1 + (2 * gas.γ) / (gas.γ + 1) * (Mn2 - 1)
end

"""
**Equation 4.10** from Anderson&Anderson

Computes the normal Mach number ratio across a _stationary_ shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
function shock_normal_mach_ratio(M, n̂, gas::CaloricallyPerfectGas)
    Mnsqr = (M ⋅ n̂)^2
    Mn2sqr = (Mnsqr + (2 / (gas.γ - 1))) / ((2 * gas.γ / (gas.γ - 1)) * Mnsqr - 1)
    mach_ratio = sqrt(Mn2sqr / Mnsqr)
    return mach_ratio
end

"""
Computes the normal velocity ratio across a _stationary_ shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.

Derived from speed of sound proportional to the square root of temperature.
Useful for computing momentum ratio.
"""
@inline function shock_normal_velocity_ratio(M, n̂, gas::CaloricallyPerfectGas)
    return sqrt(shock_temperature_ratio(M, n̂, gas)) *
           shock_normal_mach_ratio(M, n̂, gas)
end

"""
Computes the tangential mach number ratio across a _stationary_ shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.

Derived from speed of sound proportional to the square root of temperature.
"""
@inline function shock_tangent_mach_ratio(M, n̂, gas::CaloricallyPerfectGas)
    return sqrt(1 / shock_temperature_ratio(M, n̂, gas))
end

@inline function shock_normal_momentum_ratio(M, n̂, gas::CaloricallyPerfectGas)
    return shock_density_ratio(M, n̂, gas) *
           shock_normal_velocity_ratio(M, n̂, gas)
end

"""
**Equation 4.11** from Anderson&Anderson

Computes the temperature ratio across a _stationary_ shock wave.
The incoming flow has Mach number(s) ``M=[M_x, M_y]`` and the outward (away from body) normal is ``n̂``.
"""
@inline function shock_temperature_ratio(M, n̂, gas::CaloricallyPerfectGas)
    return shock_pressure_ratio(M, n̂, gas) / shock_density_ratio(M, n̂, gas)
end

"""
    state_behind(state_L, n̂, t̂; gas::CaloricallyPerfectGas)
Computes the gas state behind a _stationary_ shock wave.
The outward (away from body) normal to the shock wave is ``n̂`` and the tangent to the shock wave is ``t̂``.
"""
function state_behind(uL::ConservedProps, n̂, t̂, gas::CaloricallyPerfectGas)
    @assert ≈(t̂ ⋅ n̂, 0.0, atol = eps(Float64)) "tangent and normal vectors should be normal to each other."
    M_L = mach_number(uL, gas)
    ρ_R = uL.ρ * shock_density_ratio(M_L, n̂, gas)
    # momentum change
    ρv_n_R = (momentum_density(uL, gas) ⋅ n̂) * shock_normal_momentum_ratio(M_L, n̂, gas)
    # mach number component constant in tangential direction
    ρv_t_R = (momentum_density(uL, gas) ⋅ t̂) * shock_density_ratio(M_L, n̂, gas)
    ρv_R = ρv_n_R * n̂ + ρv_t_R * t̂
    # total internal energy change
    ρe_L = static_internal_energy_density(uL)
    # e = cT
    ρe_R =
        ρe_L *
        shock_density_ratio(M_L, n̂, gas) *
        shock_temperature_ratio(M_L, n̂, gas)
    ρE_R = ρe_R + (ρv_R ⋅ ρv_R) / (2 * ρ_R)
    return ConservedProps(ρ_R, ρv_R, ρE_R)
end

function state_behind(sL::PrimitiveProps, n̂, t̂, gas::CaloricallyPerfectGas)
    @assert ≈(t̂ ⋅ n̂, 0.0, atol = eps(Float64)) "tangent and normal vectors should be normal to each other."
    # mach number change
    M_n_R = (mach_number(sL, gas) ⋅ n̂) * shock_normal_mach_ratio(sL.M, n̂, gas)
    M_t_R =
        (mach_number(sL, gas) ⋅ t̂) * shock_tangent_mach_ratio(sL.M, n̂, gas)
    M_R = M_n_R * n̂ + M_t_R * t̂
    # density and temperature change
    ρ_R = density(sL) * shock_density_ratio(sL.M, n̂, gas)
    T_R = temperature(sL) * shock_temperature_ratio(sL.M, n̂, gas)
    return PrimitiveProps(ρ_R, M_R, T_R)
end

### COMPUTE STATES WITHOUT RESPECTING UNITS ###

function primitive_state_behind(state_L, n̂, t̂, gas::CaloricallyPerfectGas)
    @assert ≈(t̂ ⋅ n̂, 0.0, atol = eps(Float64)) "tangent and normal vectors should be normal to each other."
    M_L = state_L[2:end-1]
    # mach number change
    Mn_R = (M_L ⋅ n̂) * shock_normal_mach_ratio(M_L, n̂, gas)
    Mt_R = (M_L ⋅ t̂) * shock_tangent_mach_ratio(M_L, n̂, gas)
    M_R = Mn_R * n̂ + Mt_R * t̂
    # density and temperature change
    ρ_R = state_L[1] * shock_density_ratio(M_L, n̂, gas)
    T_R = state_L[end] * shock_temperature_ratio(M_L, n̂, gas)
    return vcat(ρ_R, M_R, T_R)
end

function conserved_state_behind(state_L, n̂, t̂, gas::CaloricallyPerfectGas)
    @assert ≈(t̂ ⋅ n̂, 0.0, atol = eps(Float64)) "tangent and normal vectors should be normal to each other."
    ρv_L = state_L[2:end-1]
    ρe_L = static_internal_energy_density(state_L[1], ρv_L, state_L[end])
    # find the speed of sound w/o units :)
    # I don't actually know how well this plays with ad tools e.g. zygote
    T_L = ρe_L / (state_L[1] * ustrip(_units_cvcp, gas.c_v))
    a_L = ustrip(u"m/s", speed_of_sound(T_L, gas))
    M_L = ρv_L / (state_L[1] * a_L)
    # density change
    ρ_R = state_L[1] * shock_density_ratio(M_L, n̂, gas)
    # momentum change
    ρv_n_R = (ρv_L ⋅ n̂) * shock_normal_momentum_ratio(M_L, n̂, gas)
    ρv_t_R = (ρv_L ⋅ t̂) * shock_density_ratio(M_L, n̂, gas)
    ρv_R = ρv_n_R * n̂ + ρv_t_R * t̂
    # total internal energy change
    ρe_R =
        ρe_L *
        shock_density_ratio(M_L, n̂, gas) *
        shock_temperature_ratio(M_L, n̂, gas)
    ρE_R = ρe_R + (ρv_R ⋅ ρv_R) / (2 * ρ_R)
    return vcat(ρ_R, ρv_R, ρE_R)
end
