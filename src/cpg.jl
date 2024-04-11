using Unitful: Temperature, Density, Velocity, @derived_dimension
using Unitful: 𝐋, 𝐓, 𝐌, 𝚯, 𝐍

@derived_dimension HeatCapacity 𝐋^2 * 𝐓^-2 * 𝚯^-1
@derived_dimension MolarMass 𝐌 * 𝐍^-1
@derived_dimension MomentumDensity 𝐌 * 𝐋^-2 * 𝐓^-1
@derived_dimension SpecificEnergy 𝐋^2 * 𝐓^-2
@derived_dimension EnergyDensity 𝐌 * 𝐋^-1 * 𝐓^-2

# default units for certain quantities.
const _units_cvcp = u"J/kg/K"
const _units_ℳ = u"kg/mol"
const _units_ρ = u"kg/m^3"
const _units_v = u"m/s"
const _units_T = u"K"

const _units_ρv = _units_ρ * _units_v
const _units_int_e = _units_cvcp * _units_T
const _units_ρE = _units_ρ * _units_int_e

"""
    CaloricallyPerfectGas
Provides the properties of a calorically perfect gas (or mixture of gases).
 - ``c_p``: Heat capacity at constant pressure
 - ``c_v``: Heat capacity at constant volume
 - ``ℳ``: Molar mass
 - ``γ``: Heat capacity ratio
 - ``R``: Gas constant
"""
struct CaloricallyPerfectGas{U1<:HeatCapacity,U2<:MolarMass}
    c_p::U1
    c_v::U1
    ℳ::U2

    γ::Float64
    R::U1
end

function CaloricallyPerfectGas(c_p, c_v, ℳ)
    q_cp = Quantity(c_p, _units_cvcp)
    q_cv = Quantity(c_v, _units_cvcp)
    q_ℳ = Quantity(ℳ, _units_ℳ)
    q_R = q_cp - q_cv
    return CaloricallyPerfectGas{typeof(q_cp),typeof(q_ℳ)}(q_cp, q_cv, q_ℳ, c_p / c_v, q_R)
end

function CaloricallyPerfectGas(
    c_p::T,
    c_v::T,
    ℳ::U,
) where {T<:HeatCapacity,U<:Unitful.MolarMass}
    R = c_p - c_v
    γ = c_p / c_v
    return CaloricallyPerfectGas{T,U}(c_p, c_v, ℳ, γ, R)
end

"""
    enthalpy(T; gas::CaloricallyPerfectGas)
Computes the static enthalpy of a calorically perfect gas at a temperature ``T``.
"""
static_enthalpy(T; gas::CaloricallyPerfectGas) = gas.c_p * Quantity(T, _units_T)
static_enthalpy(T::Temperature; gas::CaloricallyPerfectGas) = gas.c_p * T

enthalpy = static_enthalpy

"""
    internal_energy(T; gas::CaloricallyPerfectGas)
Computes the internal energy of a calorically perfect gas.
"""
internal_energy(T; gas::CaloricallyPerfectGas) = gas.c_v * Quantity(T, _units_T)
internal_energy(T::Temperature; gas::CaloricallyPerfectGas) = gas.c_v * T

"""
    speed_of_sound(T; gas::CaloricallyPerfectGas)
Computes the speed of sound in an ideal gas at a temperature ``T``. 

*We assume that the gas is a non-dispersive medium.*
"""
speed_of_sound(T; gas::CaloricallyPerfectGas) = sqrt(gas.γ * gas.R * Quantity(T, _units_T))
speed_of_sound(T::Temperature; gas::CaloricallyPerfectGas) = sqrt(gas.γ * gas.R * T)

"""
    pressure(ρ, T; gas::CaloricallyPerfectGas)
Compute the pressure in a calorically perfect gas from its density and temperature.

    pressure(ρe; gas::CaloricallyPerfectGas)
Compute the pressure in a calorically perfect gas from its internal energy density.
"""
function pressure(ρ, T; gas::CaloricallyPerfectGas)
    return (gas.γ - 1) * Quantity(ρ, _units_ρ) * internal_energy(T; gas = gas)
end

function pressure(ρ::Density, T::Temperature; gas::CaloricallyPerfectGas)
    return (gas.γ - 1) * ρ * internal_energy(T; gas = gas)
end

pressure(ρe; gas::CaloricallyPerfectGas) = (gas.γ - 1) * Quantity(ρe, _units_ρE)
pressure(ρe::EnergyDensity; gas::CaloricallyPerfectGas) = (gas.γ - 1) * ρe

"""
Properties that are easier to reason about than those in a `ConservationProps`. 
These completely determine the state of a calorically perfect gas.

 - ``ρ``: Density of the gas.
 - ``M``: The mach number, represented as a vector quantity.
 - ``T``: The absolute temperature of the gas.
"""
struct PrimitiveProps{U1<:Density,U2<:Temperature}
    ρ::U1
    M::Vector{Float64}
    T::U2
end
PrimitiveState = PrimitiveProps

"""
    PrimitiveProps(ρ::Float64, M::Vector{Float64}, T::Float64)
Construct a PrimitiveState and assign the default units.
"""
function PrimitiveProps(ρ::U, M::Vector{U}, T::U) where {U}
    return PrimitiveProps(Quantity(ρ, _units_ρ), M, Quantity(T, _units_T))
end

"""
    ConservationProps{<:Density, <:MomentumDensity, <:EnergyDensity}
The conserved quantities in the Euler equations.

 - ``ρ``: Density of the gas.
 - ``ρv``: Momentum density of the gas, represented as a vector quantity.
 - ``ρE``: The **sum** of the internal energy density of the gas and kinetic energy density of the moving gas.
"""
struct ConservationProps{U1<:Density,U2<:MomentumDensity,U3<:EnergyDensity}
    ρ::U1
    ρv::Vector{U2}
    ρE::U3
end
ConservedState = ConservationProps

"""
    PrimitiveState(state::ConservedState; gas::CaloricallyPerfectGas)
Compute the primitive state quantities from the conserved state quantities and 
the thermodynamic properies of the gas.
"""
function PrimitiveProps(state::ConservationProps; gas::CaloricallyPerfectGas)
    a = speed_of_sound(state; gas = gas)
    # mach number should be dimensionless and unitless (we avoid e.g. 2.0u"mm/m")
    return PrimitiveProps(
        state.ρ,
        uconvert.(NoUnits, state.ρv / (state.ρ * a)),
        temperature(state; gas = gas),
    )
end

"""
    ConservedState(ρ::Float64, ρv::Vector{Float64}, ρE::Float64)
Construct a ConservedState and assign the default units.
"""
function ConservationProps(ρ::T, ρv::Vector{T}, ρE::T) where {T}
    return ConservationProps(
        Quantity(ρ, _units_ρ),
        Quantity.(ρv, _units_ρv),
        Quantity(ρE, _units_ρE),
    )
end

"""
    ConservedState(state::PrimitiveProps, gas::CaloricallyPerfectGas)
Compute the conserved state quantities from primitive state quantities and
the thermodynamic properties of a gas.
"""
function ConservationProps(state::PrimitiveProps, gas::CaloricallyPerfectGas)
    v = state.M * speed_of_sound(state; gas = gas)
    e = gas.c_v * state.T
    return ConservationProps(state.ρ, state.ρ * v, state.ρ * (e + v ⋅ v / 2))
end

"""
    internal_energy_density(ρ, ρv, ρE)
    internal_energy_density(state::ConservedState)
Compute the internal energy volume density (``ρe``) from conserved state quantities.
"""
internal_energy_density(ρ, ρv, ρE) = ρE - (ρv ⋅ ρv) / (2 * ρ)
internal_energy_density(state::ConservationProps) =
    internal_energy_density(state.ρ, state.ρv, state.ρE)

"""
    internal_energy(state; gas::CaloricallyPerfectGas)
Compute the internal energy (``e``) of a gas at a given state.
"""
internal_energy(state::ConservationProps; gas::CaloricallyPerfectGas) =
    internal_energy_density(state) / state.ρ
internal_energy(state::PrimitiveProps; gas::CaloricallyPerfectGas) = gas.c_v * state.T

"""
    pressure(state; gas::CaloricallyPerfectGas)
Compute the pressure at a given state in a gas.
"""
pressure(state::ConservationProps; gas::CaloricallyPerfectGas) =
    pressure(internal_energy_density(state); gas = gas)
pressure(state::PrimitiveProps; gas::CaloricallyPerfectGas) =
    pressure(state.ρ, state.T; gas = gas)

"""
    density(state)
Compute the density at a given state in a gas.
"""
density(state::Union{ConservationProps,PrimitiveProps}) = state.ρ

"""
    temperature(state; gas::CaloricallyPerfectGas)
Compute the temperature at a given state in a gas.
"""
function temperature(state::ConservationProps; gas::CaloricallyPerfectGas)
    return internal_energy(state; gas = gas) / gas.c_v
end

temperature(state::PrimitiveProps; gas::CaloricallyPerfectGas) = state.T

"""
    speed_of_sound(state::Union{ConservedState, PrimitiveState}; gas::CaloricallyPerfectGas)
Compute the speed of sound in a gas at a given state. 

*We assume that the gas is a non-dispersive medium.*
"""
function speed_of_sound(
    state::Union{ConservationProps,PrimitiveProps};
    gas::CaloricallyPerfectGas,
)
    return speed_of_sound(temperature(state; gas = gas); gas = gas)
end

### DISRESPECT UNITS AND WORK WITH STATES AS COLLECTIONS ###

function state_to_vector(state::T) where {T<:Union{PrimitiveProps,ConservationProps}}
    return vcat(map(sym -> ustrip.(getfield(state, sym)), fieldnames(T))...)
end

"""
    primitive_state_vector(u; gas::CaloricallyPerfectGas)

Takes a vector of conserved quantities ``u=[ρ, ρv, ρE]`` and converts it into
``s=[ρ, M, T]``. 

**Assumes that everything is given in metric base units!**
"""
function primitive_state_vector(u; gas::CaloricallyPerfectGas)
    ρv = u[2:end-1]
    ρe = internal_energy_density(u[1], ρv, u[end])
    T = ρe / (u[1] * ustrip(_units_cvcp, gas.c_v))
    a = ustrip(u"m/s", speed_of_sound(T; gas = gas))
    return vcat(u[1], ρv / (u[1] * a), T)
end

function conserved_state_vector(s; gas::CaloricallyPerfectGas)
    a = ustrip(u"m/s", speed_of_sound(s[end]; gas = gas))
    ρv = s[1] * s[2:end-1] * a
    ρe = s[1] * ustrip(_units_cvcp, gas.c_v) * s[end]
    ρE = ρe + ρv ⋅ ρv / (2 * s[1])
    return vcat(s[1], ρv, ρE)
end