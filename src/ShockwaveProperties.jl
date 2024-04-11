module ShockwaveProperties

using LinearAlgebra
using Unitful
using UnitfulChainRules

# constants
export DRY_AIR
# submodules
export BilligShockParametrization
# gas properties
export CaloricallyPerfectGas, ConservationProps, PrimitiveProps
export density, temperature, pressure, speed_of_sound
export static_enthalpy, internal_energy, internal_energy_density

# shock wave jumps
export shock_density_ratio, shock_pressure_ratio, shock_temperature_ratio
export state_behind

# what if we ignored units?
export state_to_vector
export conserved_state_vector,
    conserved_state_behind, primitive_state_vector, primitive_state_behind

include("cpg.jl")

# At 300K
const DRY_AIR = CaloricallyPerfectGas(1004.9u"J/kg/K", 717.8u"J/kg/K", 0.0289647u"kg/mol")

include("billig.jl")
include("normal_shocks.jl")
include("rules.jl")

end
