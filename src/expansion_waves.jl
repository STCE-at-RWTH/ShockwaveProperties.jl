function forward_mach_angle()
end

function rear_mach_angle()
end

function expansion_wave_temperature_raio()
end

function expansion_wave_pressure_ratio()
end

function prandtl_meyer_func(M, γ)
    γ_r = (γ + 1) / (γ - 1)
    return sqrt(γ_r) * atan(sqrt((M^2 - 1) / γ_r)) - atan(sqrt(M^2 - 1))
end