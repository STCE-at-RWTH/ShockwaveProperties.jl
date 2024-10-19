using Unitful: Length, Acceleration

## DENSITY

"""
    density(state)
Compute the density at a given state in a gas.
"""
density(state) = state.ρ

## MOMENTUM
"""
    momentum_density(state, gas::CaloricallyPerfectGas)

Momentum in a gas at a given state.
"""
momentum_density(u::ConservedProps, gas = nothing) = u.ρv
function momentum_density(s::PrimitiveProps, gas::CaloricallyPerfectGas)
    v = velocity(s, gas)
    return density(s) .* v
end

## VELOCITY

"""
    velocity(state, gas::CaloricallyPerfectGas)

Velocity in a gas at a given state.
"""
function velocity(u::ConservedProps, gas = nothing)
    return momentum_density(u) ./ density(u)
end

function velocity(s::PrimitiveProps, gas::CaloricallyPerfectGas)
    return mach_number(s, gas) .* speed_of_sound(s, gas)
end

## MACH NUMBER

"""
    mach_number(state, gas::CaloricallyPerfectGas)

Mach number in a gas at a given state.
"""
function mach_number(u::ConservedProps, gas::CaloricallyPerfectGas)
    return uconvert.(NoUnits, velocity(u) ./ speed_of_sound(u, gas))
end
mach_number(s::PrimitiveProps, gas = nothing) = s.M

## SPECIFIC STATIC ENTHALPY AND INTERNAL ENERGY

"""
    specific_static_enthalpy(T, gas::CaloricallyPerfectGas)
Computes the static enthalpy of a calorically perfect gas at a temperature ``T`` (default in Kelvin).
"""
function specific_static_enthalpy(T::Real, gas::CaloricallyPerfectGas)
    return gas.c_p * Quantity(T, _units_T)
end
specific_static_enthalpy(T::Temperature, gas::CaloricallyPerfectGas) = gas.c_p * T

# there might be a better way for different property sets...
"""
    specific_static_enthalpy(state, gas::CaloricallyPerfectGas)

Compute the specific static enthalpy.
"""
function specific_static_enthalpy(s::PrimitiveProps, gas::CaloricallyPerfectGas)
    return specific_static_enthalpy(temperature(s, gas), gas)
end

function specific_static_enthalpy(u::ConservedProps, gas::CaloricallyPerfectGas)
    return gas.γ * specific_static_internal_energy(u)
end

"""
    specific_internal_energy(T, gas::CaloricallyPerfectGas)
Computes the specific internal energy of a calorically perfect gas.
"""
function specific_static_internal_energy(T::Real, gas::CaloricallyPerfectGas)
    gas.c_v * Quantity(T, _units_T)
end
specific_static_internal_energy(T::Temperature, gas::CaloricallyPerfectGas) = gas.c_v * T

"""
    specific_internal_energy(state, gas::CaloricallyPerfectGas)
Compute the internal energy (``e``) of a gas at a given state.
"""
function specific_static_internal_energy(u::ConservedProps, gas = nothing)
    return static_internal_energy_density(u) / density(u)
end

function specific_static_internal_energy(s::PrimitiveProps, gas::CaloricallyPerfectGas)
    return gas.c_v * temperature(s)
end

## STATIC INTERNAL ENERGY DENSITY

"""
    static_internal_energy_density(ρ, ρv, ρE)
    static_internal_energy_density(state, gas)

Compute the static internal energy density (``ρe``) from conserved state quantities.
"""
static_internal_energy_density(ρ, ρv, ρE) = ρE - (ρv ⋅ ρv) / (2 * ρ)

function static_internal_energy_density(u::ConservedProps, gas = nothing)
    return static_internal_energy_density(
        density(u),
        momentum_density(u),
        total_internal_energy_density(u),
    )
end

function static_internal_energy_density(s::PrimitiveProps, gas::CaloricallyPerfectGas)
    return density(s) * gas.c_v * temperature(s)
end

## TOTAL INTERNAL ENERGY DENSITY

"""
    total_internal_energy_density(state, gas::CaloricallyPerfectGas)

Compute the total internal energy density at a given state in a gas.
"""
function total_internal_energy_density(u::ConservedProps, gas = nothing)
    return u.ρE
end

function total_internal_energy_density(s::PrimitiveProps, gas::CaloricallyPerfectGas)
    v = velocity(s, gas)
    specific_KE = v ⋅ v / 2
    return static_internal_energy_density(s, gas) + density(s) * specific_KE
end

## TOTAL ENTHALPY
"""
    specific_total_enthalpy(state, gas::CaloricallyPerfectGas)

Compute the specific total enthalpy at a state in a gas. This quantity is constant along the stagnation line.
"""
function specific_total_enthalpy(state, gas::CaloricallyPerfectGas)
    v = velocity(state, gas)
    return specific_static_enthalpy(state, gas) + (v ⋅ v) / 2
end

""" 
    total_enthalpy_density(ρ, ρv, ρE, gas::CaloricallyPerfectGas)

Compute the total enthalpy at ``u=[ρ, ρv, ρE]``.
"""
function total_enthalpy_density(
    ρ::Density,
    ρv::AbstractVector{<:MomentumDensity},
    ρE::EnergyDensity,
    gas::CaloricallyPerfectGas,
)
    return ρE + pressure(static_internal_energy_density(ρ, ρv, ρE), gas)
end

function total_enthalpy_density(
    ρ::Real,
    ρv::AbstractVector{<:Real},
    ρE::Real,
    gas::CaloricallyPerfectGas,
)
    return ρE + ustrip(pressure(static_internal_energy_density(ρ, ρv, ρE), gas))
end

function total_enthalpy_density(state, gas::CaloricallyPerfectGas)
    return total_internal_energy_density(state, gas) + pressure(state, gas)
end

## TEMPERATURE

"""
    temperature(state, gas::CaloricallyPerfectGas)
Compute the temperature at a given state in a gas.
"""
function temperature(u::ConservedProps, gas::CaloricallyPerfectGas)
    return specific_static_internal_energy(u, gas) / gas.c_v
end
temperature(s::PrimitiveProps, gas = nothing) = s.T

## PRESSURE

"""
    pressure(ρ, T, gas::CaloricallyPerfectGas)
Compute the pressure in a calorically perfect gas from its density and temperature.
"""
function pressure(ρ::Real, T::Real, gas::CaloricallyPerfectGas)
    return Quantity(ρ, _units_ρ) * gas.R * Quantity(T, _units_T)
end

pressure(ρ::Density, T::Temperature, gas::CaloricallyPerfectGas) = ρ * gas.R * T

"""
    pressure(ρe, gas::CaloricallyPerfectGas)

Compute the pressure in a calorically perfect gas from its static internal energy density.
"""
pressure(ρe::Real, gas::CaloricallyPerfectGas) = (gas.γ - 1) * Quantity(ρe, _units_ρE)
pressure(ρe::EnergyDensity, gas::CaloricallyPerfectGas) = (gas.γ - 1) * ρe

# we don't want to restrict this to avoid quantities
"""
    pressure(ρ, ρv, ρE, gas::CaloricallyPerfectGas)

Compute the pressure at `u=[ρ, ρv, ρE]`. 
"""
function pressure(ρ, ρv, ρE, gas::CaloricallyPerfectGas)
    return pressure(static_internal_energy_density(ρ, ρv, ρE), gas)
end

"""
    pressure(state, gas::CaloricallyPerfectGas)
Compute the pressure at a given state in a gas.
"""
function pressure(u::ConservedProps, gas::CaloricallyPerfectGas)
    return pressure(static_internal_energy_density(u), gas)
end

function pressure(s::PrimitiveProps, gas::CaloricallyPerfectGas)
    return pressure(density(s), temperature(s), gas)
end

## SPEED OF SOUND

"""
    speed_of_sound(T, gas::CaloricallyPerfectGas)
Computes the speed of sound in an ideal gas at a temperature ``T``. 

*We assume that the gas is a non-dispersive medium.*
"""
function speed_of_sound(T::Real, gas::CaloricallyPerfectGas)
    return sqrt(gas.γ * gas.R * Quantity(T, _units_T))
end
speed_of_sound(T::Temperature, gas::CaloricallyPerfectGas) = sqrt(gas.γ * gas.R * T)

"""
    speed_of_sound(ρ, P, gas::CaloricallyPerfectGas)
Comptue the speed of sound in an ideal gas at density ``ρ`` and pressure ``P``.

*We assume that the gas is a non-dispersive medium.*
"""
function speed_of_sound(ρ::Real, P::Real, gas::CaloricallyPerfectGas)
    return sqrt(gas.γ * Quantity(P, _units_P) / Quantity(ρ, _units_ρ))
end
speed_of_sound(ρ::Density, P::Pressure, gas::CaloricallyPerfectGas) = sqrt(gas.γ * P / ρ)

"""
    speed_of_sound(ρ, ρv, ρE, gas::CaloricallyPerfectGas)

Compute the speed of sound from conserved state properties.
"""
function speed_of_sound(ρ::Real, ρv, ρE::Real, gas::CaloricallyPerfectGas)
    return speed_of_sound(ρ, ustrip(pressure(ρ, ρv, ρE, gas)), gas)
end

function speed_of_sound(ρ::Density, ρv, ρE::EnergyDensity, gas::CaloricallyPerfectGas)
    return speed_of_sound(ρ, pressure(ρ, ρv, ρE, gas), gas)
end

"""
    speed_of_sound(state, gas::CaloricallyPerfectGas)
Compute the speed of sound in a gas at a given state. 

*We assume that the gas is a non-dispersive medium.*
"""
function speed_of_sound(s::PrimitiveProps, gas::CaloricallyPerfectGas)
    return speed_of_sound(temperature(s, gas), gas)
end

function speed_of_sound(u::ConservedProps, gas::CaloricallyPerfectGas)
    return speed_of_sound(density(u), pressure(u, gas), gas)
end

## KINEMATIC VISCOSITY

function kinematic_viscosity(s, gas::CaloricallyPerfectGas)
    return gas.μ / density(s)
end

## THERMAL DIFFUSIVITY

function thermal_diffusivity(s, gas::CaloricallyPerfectGas)
    return gas.k / (gas.c_p * density(s))
end

## DIMENSIONLESS NUMBERS THAT DEPEND ON FLUID STATES

function reynolds_number(s, gas::CaloricallyPerfectGas, D::T) where {T<:Length}
    return uconvert(Unitful.NoUnits, density(s) * norm(velocity(s, gas)) * D / gas.μ)
end

function froude_number(
    s,
    gas,
    external_force_field::U1,
    L::U2,
) where {U1<:Acceleration,U2<:Length}
    return uconvert(
        Unitful.NoUnits,
        norm(velocity(s, gas)) / sqrt(norm(external_force_field) * L),
    )
end

function fourier_number(
    s,
    gas::CaloricallyPerfectGas,
    t::U1,
    L::U2,
) where {U1<:Time,U2<:Length}
    return uconvert(Unitful.NoUnits, thermal_diffusivity(s, gas) * t / L^2)
end


