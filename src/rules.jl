# rrules for AD
using Unitful
using ChainRulesCore
using UnitfulChainRules

# ChainRulesCore.debug_mode() = true

function ChainRulesCore.rrule(::RuleConfig{>:HasReverseMode}, ::Type{ConservationProps{U1, U2, U3}}, ρ::U1, ρv::Vector{U2}, ρE::U3) where {U1<:MomentumDensity, U2, U3}
    ConservedState_back(Δs) = begin
        return NoTangent(), Δs.ρ, Δs.ρv, Δs.ρE
    end
    return ConservedState{U1, U2, U3}(ρ, ρv, ρE), ConservedState_back
end

function ChainRulesCore.rrule(::RuleConfig{>:HasReverseMode}, ::Type{ConservationProps{U1, U2, U3}}, s::PrimitiveProps{S1, S2}, gas::CaloricallyPerfectGas) 



# function ChainRulesCore.rrule(::typeof(speed_of_sound), )

# ShockwaveProperties.ConservedState{
#     Unitful.Quantity{Float64, 𝐌 𝐋^-3, Unitful.FreeUnits{(kg, m^-3), 𝐌 𝐋^-3, nothing}}, 
#     Unitful.Quantity{Float64, 𝐌 𝐋^-2 𝐓^-1, Unitful.FreeUnits{(kg^1/2, J^1/2, m^-3), 𝐌 𝐋^-2 𝐓^-1, nothing}}, 
#     Unitful.Quantity{Float64, 𝐌 𝐋^-1 𝐓^-2, Unitful.FreeUnits{(J, m^-3), 𝐌 𝐋^-1 𝐓^-2, nothing}}}

# ShockwaveProperties.ConservedState{
#     Unitful.Quantity{Float64, 𝐌 𝐋^-3, Unitful.FreeUnits{(kg, m^-3), 𝐌 𝐋^-3, nothing}}, 
#     Unitful.Quantity{Float64, 𝐌 𝐋^-2 𝐓^-1, Unitful.FreeUnits{(kg, m^-2, s^-1), 𝐌 𝐋^-2 𝐓^-1, nothing}}, 
#     Unitful.Quantity{Float64, 𝐌 𝐋^-1 𝐓^-2, Unitful.FreeUnits{(J, m^-3), 𝐌 𝐋^-1 𝐓^-2, nothing}}}