using Unitful: Temperature, Pressure, Density, Velocity, @derived_dimension
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
const _units_P = u"Pa"

const _units_ρv = _units_ρ * _units_v
const _units_int_e = _units_cvcp * _units_T
const _units_ρE = _units_ρ * _units_int_e

"""
    CaloricallyPerfectGas
Provides the properties of a calorically perfect gas (or mixture of gases).
 - ``c_p``: Specific heat capacity at constant pressure
 - ``c_v``: Specific heat capacity at constant volume
 - ``ℳ``: Molar mass
 - ``γ``: Heat capacity ratio
 - ``R``: Specific gas constant
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
Properties that are easier to reason about than those in a `ConservedProps`. 
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
Construct a PrimitiveProps and assign the default units.
"""
function PrimitiveProps(ρ::U, M::Vector{U}, T::U) where {U}
    return PrimitiveProps(Quantity(ρ, _units_ρ), M, Quantity(T, _units_T))
end

"""
    PrimitiveProps(s::AbstractVector)
Construct a PrimitiveProps from a vector and assign default units.
"""
function PrimitiveProps(s::AbstractVector)
    return PrimitiveProps(Quantity(s[1], _units_ρ), s[2:end-1], Quantity(s[end], _units_T))
end

"""
    ConservedProps{<:Density, <:MomentumDensity, <:EnergyDensity}
The conserved quantities in the Euler equations.

 - ``ρ``: Density of the gas.
 - ``ρv``: Momentum density of the gas, represented as a vector quantity.
 - ``ρE``: The **sum** of the internal energy density of the gas and kinetic energy density of the moving gas.
"""
struct ConservedProps{U1<:Density,U2<:MomentumDensity,U3<:EnergyDensity}
    ρ::U1
    ρv::Vector{U2}
    ρE::U3
end
ConservedState = ConservedProps

"""
    ConservedProps(ρ::Float64, ρv::Vector{Float64}, ρE::Float64)
Construct a ConservedState and assign the default units.
"""
function ConservedProps(ρ::T, ρv::Vector{T}, ρE::T) where {T}
    return ConservedProps(
        Quantity(ρ, _units_ρ),
        Quantity.(ρv, _units_ρv),
        Quantity(ρE, _units_ρE),
    )
end

"""
    ConservedProps(p::AbstractVector)
Construct a ConservedProps from a vector and assign default units.
"""
function ConservedProps(u::AbstractVector)
    return ConservedProps(
        Quantity(u[1], _units_ρ),
        Quantity.(u[2:end-1], _units_ρv),
        Quantity(u[end], _units_ρE),
    )
end

## DENSITY

"""
    density(state)
Compute the density at a given state in a gas.
"""
density(state) = state.ρ

## MOMENTUM
"""
    momentum_density(state; gas::CaloricallyPerfectGas)

Momentum in a gas at a given state.
"""
momentum_density(u::ConservedProps; gas = nothing) = u.ρv
function momentum_density(s::PrimitiveProps; gas::CaloricallyPerfectGas)
    v = velocity(s; gas)
    return density(s) * v
end

## VELOCITY

"""
    velocity(state; gas::CaloricallyPerfectGas)

Velocity in a gas at a given state.
"""
function velocity(u::ConservedProps; gas = nothing)
    return momentum_density(u) / density(u)
end

function velocity(s::PrimitiveProps; gas::CaloricallyPerfectGas)
    return mach_number(s; gas) * speed_of_sound(s; gas)
end

## MACH NUMBER

"""
    mach_number(state; gas::CaloricallyPerfectGas)

Mach number in a gas at a given state.
"""
function mach_number(u::ConservedProps; gas::CaloricallyPerfectGas)
    return map(v->ustrip(NoUnits, v), velocity(u) / speed_of_sound(u; gas))
end
mach_number(s::PrimitiveState; gas = nothing) = s.M

## SPECIFIC STATIC ENTHALPY AND INTERNAL ENERGY

"""
    specific_static_enthalpy(T; gas::CaloricallyPerfectGas)
Computes the static enthalpy of a calorically perfect gas at a temperature ``T`` (default in Kelvin).
"""
function specific_static_enthalpy(T::Number; gas::CaloricallyPerfectGas)
    return gas.c_p * Quantity(T, _units_T)
end
specific_static_enthalpy(T::Temperature; gas::CaloricallyPerfectGas) = gas.c_p * T

# there might be a better way for different property sets...
"""
    specific_static_enthalpy(state; gas::CaloricallyPerfectGas)

Compute the specific static enthalpy.
"""
function specific_static_enthalpy(state; gas::CaloricallyPerfectGas)
    return specific_static_enthalpy(temperature(state; gas); gas)
end

enthalpy = specific_static_enthalpy

"""
    specific_internal_energy(T; gas::CaloricallyPerfectGas)
Computes the specific internal energy of a calorically perfect gas.
"""
specific_internal_energy(T; gas::CaloricallyPerfectGas) = gas.c_v * Quantity(T, _units_T)
specific_internal_energy(T::Temperature; gas::CaloricallyPerfectGas) = gas.c_v * T

"""
    specific_internal_energy(state; gas::CaloricallyPerfectGas)
Compute the internal energy (``e``) of a gas at a given state.
"""
function specific_internal_energy(u::ConservedProps; gas = nothing)
    return static_internal_energy_density(u) / density(u)
end

function specific_internal_energy(s::PrimitiveProps; gas::CaloricallyPerfectGas)
    return gas.c_v * temperature(s)
end

## STATIC INTERNAL ENERGY DENSITY

"""
    static_internal_energy_density(ρ, ρv, ρE)
    static_internal_energy_density(state)

Compute the static internal energy density (``ρe``) from conserved state quantities.
"""
static_internal_energy_density(ρ, ρv, ρE) = ρE - (ρv ⋅ ρv) / (2 * ρ)

function static_internal_energy_density(u::ConservedProps; gas = nothing)
    return static_internal_energy_density(
        density(u),
        momentum_density(u),
        total_internal_energy_density(u),
    )
end

function static_internal_energy_density(s::PrimitiveProps; gas::CaloricallyPerfectGas)
    return density(s) * gas.c_v * temperature(s)
end

## TOTAL INTERNAL ENERGY DENSITY

"""
    total_internal_energy_density(state; gas::CaloricallyPerfectGas)

Compute the total internal energy density at a given state in a gas.
"""
function total_internal_energy_density(u::ConservedProps; gas = nothing)
    return u.ρE
end

function total_internal_energy_density(s::PrimitiveProps; gas::CaloricallyPerfectGas)
    v = velocity(s; gas)
    specific_KE = v ⋅ v / 2
    return static_internal_energy_density(s; gas) + density(s) * specific_KE
end

## TOTAL ENTHALPY DENSITY
""" 
    total_enthalpy_density(ρ, ρv, ρE; gas::CaloricallyPerfectGas)

Compute the total enthalpy at ``u=[ρ, ρv, ρE]``.
"""
function total_enthalpy_density(ρ, ρv, ρE; gas::CaloricallyPerfectGas)
    return ρE + pressure(static_internal_energy_density(ρ, ρv, ρE); gas)
end

function total_enthalpy_density(state; gas::CaloricallyPerfectGas)
    return total_internal_energy_density(state; gas) + pressure(state; gas)
end

## TEMPERATURE

"""
    temperature(state; gas::CaloricallyPerfectGas)
Compute the temperature at a given state in a gas.
"""
function temperature(u::ConservedProps; gas::CaloricallyPerfectGas)
    return specific_internal_energy(u; gas) / gas.c_v
end
temperature(s::PrimitiveProps; gas = nothing) = s.T

## PRESSURE

"""
    pressure(ρ, T; gas::CaloricallyPerfectGas)
Compute the pressure in a calorically perfect gas from its density and temperature.
"""
function pressure(ρ, T; gas::CaloricallyPerfectGas)
    return Quantity(ρ, _units_ρ) * gas.R * Quantity(T, _units_T)
end

pressure(ρ::Density, T::Temperature; gas::CaloricallyPerfectGas) = ρ * gas.R * T

"""
    pressure(ρe; gas::CaloricallyPerfectGas)
Compute the pressure in a calorically perfect gas from its static internal energy density.
"""
pressure(ρe; gas::CaloricallyPerfectGas) = (gas.γ - 1) * Quantity(ρe, _units_ρE)
pressure(ρe::EnergyDensity; gas::CaloricallyPerfectGas) = (gas.γ - 1) * ρe

"""
    pressure(state; gas::CaloricallyPerfectGas)
Compute the pressure at a given state in a gas.
"""
function pressure(u::ConservedProps; gas::CaloricallyPerfectGas)
    return pressure(static_internal_energy_density(u); gas)
end

function pressure(s::PrimitiveProps; gas::CaloricallyPerfectGas)
    return pressure(density(s), temperature(s); gas)
end

## SPEED OF SOUND

"""
    speed_of_sound(T; gas::CaloricallyPerfectGas)
Computes the speed of sound in an ideal gas at a temperature ``T``. 

*We assume that the gas is a non-dispersive medium.*
"""
function speed_of_sound(T::Number; gas::CaloricallyPerfectGas)
    return sqrt(gas.γ * gas.R * Quantity(T, _units_T))
end
speed_of_sound(T::Temperature; gas::CaloricallyPerfectGas) = sqrt(gas.γ * gas.R * T)

"""
    speed_of_sound(ρ, P; gas::CaloricallyPerfectGas)
Comptue the speed of sound in an ideal gas at density ``ρ`` and pressure ``P``.

*We assume that the gas is a non-dispersive medium.*
"""
function speed_of_sound(ρ, P; gas::CaloricallyPerfectGas)
    return sqrt(gas.γ * Quantity(P, _units_P) / Quantity(ρ, _units_ρ))
end
speed_of_sound(ρ::Density, P::Pressure; gas::CaloricallyPerfectGas) = sqrt(gas.γ * P / ρ)

"""
    speed_of_sound(state::Union{ConservedState, PrimitiveState}; gas::CaloricallyPerfectGas)
Compute the speed of sound in a gas at a given state. 

*We assume that the gas is a non-dispersive medium.*
"""
function speed_of_sound(state; gas::CaloricallyPerfectGas)
    return speed_of_sound(temperature(state; gas); gas)
end

## CONVERT PROPERTY REPRESENTATIONS

"""
    PrimitiveProps(state::ConservedProps; gas::CaloricallyPerfectGas)
Compute the primitive state quantities from the conserved state quantities and 
the thermodynamic properies of the gas.
"""
function PrimitiveProps(u::ConservedProps; gas::CaloricallyPerfectGas)
    a = speed_of_sound(u; gas)
    # mach number should be dimensionless and unitless (we avoid e.g. 2.0u"mm/m")
    return PrimitiveProps(u.ρ, uconvert.(NoUnits, u.ρv / (u.ρ * a)), temperature(u; gas))
end

"""
    ConservedState(state::PrimitiveProps, gas::CaloricallyPerfectGas)
Compute the conserved state quantities from primitive state quantities and
the thermodynamic properties of a gas.
"""
function ConservedProps(s::PrimitiveProps; gas::CaloricallyPerfectGas)
    v = s.M * speed_of_sound(s; gas)
    e = gas.c_v * temperature(s; gas)
    return ConservedProps(s.ρ, s.ρ * v, s.ρ * (e + v ⋅ v / 2))
end

### DISRESPECT UNITS AND WORK WITH STATES AS COLLECTIONS ###

function state_to_vector(state::T) where {T}
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
    ρe = static_internal_energy_density(u[1], ρv, u[end])
    T = ρe / (u[1] * ustrip(_units_cvcp, gas.c_v))
    a = ustrip(u"m/s", speed_of_sound(T; gas))
    return vcat(u[1], ρv / (u[1] * a), T)
end

"""
    conserved_state_vector(u; gas::CaloricallyPerfectGas)

Takes a vector of primitive quantities ``s=[ρ, M, T]`` 
    and converts it into ``u=[ρ, ρv, ρE]``. 

**Assumes that everything is given in metric base units!**
"""
function conserved_state_vector(s; gas::CaloricallyPerfectGas)
    a = ustrip(u"m/s", speed_of_sound(s[end]; gas))
    ρv = s[1] * s[2:end-1] * a
    ρe = s[1] * ustrip(_units_cvcp, gas.c_v) * s[end]
    ρE = ρe + ρv ⋅ ρv / (2 * s[1])
    return vcat(s[1], ρv, ρE)
end