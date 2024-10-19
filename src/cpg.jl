using Unitful: Length, Time, MolarMass, DynamicViscosity
@derived_dimension HeatCapacity 𝐋^2 * 𝐓^-2 * 𝚯^-1 true
@derived_dimension ThermalConductivity 𝐌 * 𝐋 * 𝐓^-3 * 𝚯^-1 true
# @derived_dimension MolarMass 𝐌 * 𝐍^-1 true
@derived_dimension MomentumDensity 𝐌 * 𝐋^-2 * 𝐓^-1 true
@derived_dimension SpecificEnergy 𝐋^2 * 𝐓^-2 true
@derived_dimension EnergyDensity 𝐌 * 𝐋^-1 * 𝐓^-2 true

# default units for certain quantities.
const _units_cvcp = u"J/kg/K"
const _units_ℳ = u"kg/mol"
const _units_ρ = u"kg/m^3"
const _units_v = u"m/s"
const _units_T = u"K"
const _units_P = u"Pa"
const _units_k = u"W/m/K"
const _units_μ = u"Pa*s"

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
struct CaloricallyPerfectGas{
    U1<:HeatCapacity,
    U2<:MolarMass,
    U3<:ThermalConductivity,
    U4<:DynamicViscosity,
}
    c_p::U1
    c_v::U1
    ℳ::U2

    γ::Float64
    R::U1

    k::U3
    μ::U4
end

function CaloricallyPerfectGas(c_p, c_v, ℳ, k, μ)
    q_cp = Quantity(c_p, _units_cvcp)
    q_cv = Quantity(c_v, _units_cvcp)
    q_ℳ = Quantity(ℳ, _units_ℳ)
    q_R = q_cp - q_cv
    q_k = Quantity(k, _units_k)
    q_μ = Quantity(k, _units_μ)
    return CaloricallyPerfectGas{typeof(q_cp),typeof(q_ℳ),typeof(q_k),typeof(q_μ)}(
        q_cp,
        q_cv,
        q_ℳ,
        c_p / c_v,
        q_R,
        q_k,
        q_μ,
    )
end

function CaloricallyPerfectGas(
    c_p::T,
    c_v::T,
    ℳ::U,
    k::V,
    μ::W,
) where {T<:HeatCapacity,U<:MolarMass,V<:ThermalConductivity,W<:DynamicViscosity}
    R = c_p - c_v
    γ = c_p / c_v
    return CaloricallyPerfectGas{T,U,V,W}(c_p, c_v, ℳ, γ, R, k, μ)
end

## DIMENSIONLESS NUMBERS THAT DEPEND ON GAS PROPERTIES

function prandtl_number(gas::CaloricallyPerfectGas)
    return uconvert(Unitful.NoUnits, gas.c_p * gas.μ / gas.k)
end

## STATES

"""
    PrimitiveProps{N, DTYPE, Q1<:Density, Q2<:Temperature}

Properties that are easier to reason about than those in a `ConservedProps`. 
These completely determine the state of a calorically perfect gas.

 - ``ρ``: Density of the gas.
 - ``M``: The mach number, represented as a statically sized-vector quantity.
 - ``T``: The absolute temperature of the gas.
"""
struct PrimitiveProps{N,DTYPE,Q1<:Density{DTYPE},Q2<:Temperature{DTYPE}}
    ρ::Q1
    M::SVector{N,DTYPE}
    T::Q2
end

"""
    PrimitiveProps(ρ::Real, M, T::Real)
Construct a PrimitiveProps and assign the default units.
"""
function PrimitiveProps(ρ, M, T)
    return PrimitiveProps(
        Quantity(ρ, _units_ρ),
        SVector{length(M)}(M),
        Quantity(T, _units_T),
    )
end

"""
    PrimitiveProps(ρ::Density, M, T::Temperature)

Construct a `PrimitiveProps` and strip any vestigal units from the mach number.
"""
function PrimitiveProps(ρ::Density, M, T::Temperature)
    return PrimitiveProps(ρ, SVector{length(M)}(uconvert.(NoUnits, M)), T)
end

"""
    PrimitiveProps(ρ::Density, M, P::Pressure, gas::CaloricallyPerfectGas)

Construct a PrimitiveProps from `[ρ, M, P]`. 
"""
function PrimitiveProps(ρ::Density, M, P::Pressure, gas::CaloricallyPerfectGas)
    T = P / (ρ * gas.R)
    return PrimitiveProps(ρ, SVector{length(M)}(uconvert.(NoUnits, M)), T)
end

"""
    PrimitiveProps(s::AbstractVector)
Construct a PrimitiveProps from a vector and assign default units.
"""
function PrimitiveProps(s::AbstractVector)
    N = length(u)
    idcs = SVector{N - 2}(ntuple(i -> i + 1, N - 2))
    return PrimitiveProps(Quantity(s[1], _units_ρ), s[idcs], Quantity(s[end], _units_T))
end

# convert to different units

function Base.convert(
    ::Type{PrimitiveProps{N,T,A1,A2}},
    x::PrimitiveProps{N,T,B1,B2},
) where {N,T,A1<:Density,A2<:Temperature,B1<:Density,B2<:Temperature}
    PrimitiveProps(convert(A1, density(x)), mach_number(x), convert(A2, temperature(x)))
end

"""
    ConservedProps{N, DTYPE, <:Density, <:MomentumDensity, <:EnergyDensity}
The conserved quantities in the Euler equations.

 - ``ρ``: Density of the gas.
 - ``ρv``: Momentum density of the gas, represented as a statically-sized quantity.
 - ``ρE``: The **sum** of the internal energy density of the gas and kinetic energy density of the moving gas.
"""
struct ConservedProps{
    N,
    DTYPE,
    U1<:Density{DTYPE},
    U2<:MomentumDensity{DTYPE},
    U3<:EnergyDensity{DTYPE},
}
    ρ::U1
    ρv::SVector{N,U2}
    ρE::U3

    function ConservedProps(
        ρ::U1,
        ρv::Union{StaticVector{N,U2},NTuple{N,U2}},
        ρE::U3,
    ) where {N,DTYPE,U1<:Density{DTYPE},U2<:MomentumDensity{DTYPE},U3<:EnergyDensity{DTYPE}}
        return new{N,DTYPE,U1,U2,U3}(ρ, SVector{N}(ρv), ρE)
    end
end

"""
    ConservedProps(u::AbstractVector{<:Real})
Construct a ConservedProps from a vector and assign default units.
"""
function ConservedProps(u::AbstractVector{<:Real})
    N = length(u)
    idcs = SVector{N - 2}(ntuple(i -> i + 1, N - 2))
    return ConservedProps(
        Quantity(u[1], _units_ρ),
        Quantity.(u[idcs], _units_ρv),
        Quantity(u[end], _units_ρE),
    )
end

# exclude the unitful constructor
"""
    ConservedProps(ρ, ρv, ρE)

Construct a ConservedProps from the individual components and assign default units.
"""
function ConservedProps(ρ::Real, ρv, ρE::Real)
    return ConservedProps(
        Quantity(ρ, _units_ρ),
        Quantity.(SVector{length(ρv)}(ρv), _units_ρv),
        Quantity(ρE, _units_ρE),
    )
end

function Base.convert(
    ::Type{ConservedProps{N,T,A1,A2,A3}},
    x::ConservedProps{N,T,B1,B2,B3},
) where {
    N,
    T,
    A1<:Density,
    A2<:MomentumDensity,
    A3<:EnergyDensity,
    B1<:Density,
    B2<:MomentumDensity,
    B3<:EnergyDensity,
}
    return ConservedProps(
        convert(A1, density(x)),
        convert.(A2, momentum_density(x)),
        convert(A3, total_internal_energy_density(x)),
    )
end

## CONVERT PROPERTY REPRESENTATIONS

"""
    PrimitiveProps(state::ConservedProps, gas::CaloricallyPerfectGas)
Compute the primitive state quantities from the conserved state quantities and 
the thermodynamic properies of the gas.
"""
function PrimitiveProps(u::ConservedProps, gas::CaloricallyPerfectGas)
    a = speed_of_sound(u, gas)
    # mach number should be dimensionless and unitless (we avoid e.g. 2.0u"mm/m")
    return PrimitiveProps(u.ρ, uconvert.(NoUnits, u.ρv ./ (u.ρ .* a)), temperature(u, gas))
end

"""
    ConservedState(state::PrimitiveProps, gas::CaloricallyPerfectGas)
Compute the conserved state quantities from primitive state quantities and
the thermodynamic properties of a gas.
"""
function ConservedProps(s::PrimitiveProps, gas::CaloricallyPerfectGas)
    v = s.M .* speed_of_sound(s, gas)
    e = gas.c_v .* temperature(s, gas)
    return ConservedProps(s.ρ, s.ρ .* v, s.ρ * (e + v ⋅ v / 2))
end

## TEST EQUALITY

function Base.isapprox(a::ConservedProps, b::ConservedProps; kwargs...)
    return (
        isapprox(a.ρ, b.ρ; kwargs...) &&
        all(isapprox.(a.ρv, b.ρv; kwargs...)) &&
        isapprox(a.ρE, b.ρE; kwargs...)
    )
end

function Base.:(==)(a::ConservedProps, b::ConservedProps)
    return (a.ρ == b.ρ && all(a.ρv .== b.ρv) && a.ρE == b.ρE)
end

function Base.isapprox(a::PrimitiveProps, b::PrimitiveProps; kwargs...)
    return (
        isapprox(a.ρ, b.ρ; kwargs...) &&
        all(isapprox.(a.M, b.M; kwargs...)) &&
        isapprox(a.T, b.T; kwargs...)
    )
end

function Base.:(==)(a::PrimitiveProps, b::PrimitiveProps)
    return (a.ρ == b.ρ && all(a.M .== b.M) && a.T == b.T)
end

### DISRESPECT UNITS AND WORK WITH STATES AS COLLECTIONS ###

function state_to_vector(state::ConservedProps{N,T,U1,U2,U3}) where {N,T,U1,U2,U3}
    SVector{N + 2}(
        ustrip(_units_ρ, density(state)),
        ustrip.(_units_ρv, momentum_density(state))...,
        ustrip(_units_ρE, total_internal_energy_density(state)),
    )
end

function state_to_vector(state::PrimitiveProps{N,T,U1,U2}) where {N,T,U1,U2}
    SVector{N + 2}(
        ustrip(_units_ρ, density(state)),
        mach_number(state)...,
        ustrip(_units_T, temperature(state)),
    )
end

"""
    primitive_state_vector(u, gas::CaloricallyPerfectGas)

Takes a vector of conserved quantities ``u=[ρ, ρv, ρE]`` and converts it into
``s=[ρ, M, T]``. 

**Assumes that everything is given in metric base units!**
"""
function primitive_state_vector(u, gas::CaloricallyPerfectGas)
    ρv = u[2:end-1]
    ρe = static_internal_energy_density(u[1], ρv, u[end])
    T = ρe / (u[1] * ustrip(_units_cvcp, gas.c_v))
    a = ustrip(u"m/s", speed_of_sound(T, gas))
    return SVector{length(u)}(u[1], (ρv / (u[1] * a))..., T)
end

"""
    conserved_state_vector(u, gas::CaloricallyPerfectGas)

Takes a vector of primitive quantities ``s=[ρ, M, T]`` 
    and converts it into ``u=[ρ, ρv, ρE]``. 

**Assumes that everything is given in metric base units!**
"""
function conserved_state_vector(s, gas::CaloricallyPerfectGas)
    a = ustrip(u"m/s", speed_of_sound(s[end], gas))
    ρv = s[1] * s[2:end-1] * a
    ρe = s[1] * ustrip(_units_cvcp, gas.c_v) * s[end]
    ρE = ρe + ρv ⋅ ρv / (2 * s[1])
    return SVector{length(s)}(s[1], ρv..., ρE)
end