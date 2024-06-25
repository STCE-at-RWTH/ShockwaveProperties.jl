# ShockwaveProperties.jl

---

This package provides routines for computing property changes across shock waves in perfect gases. The ultimate goal is to support a solver that can be used to analyze complex interactions of shock waves in dynamic systems.

## Usage

The `ConservedProps` and `PrimitiveProps` types allow for conversion between specifying gas properties as $\rho, \vec M, T$ and $\rho, \rho \vec v,\rho E$. These types will take on the SI standard units by default.

```julia
free_stream = PrimitiveProps(1.225u"kg/m^3", [2.0, 0.0], 300u"K")
u = ConservedProps(free_stream, DRY_AIR)
```

Property computation and conversion can be done via a myriad of functions. Any arguments given as `Float64` are assumed to be in SI base units (e.g. `Pa`, not `atm` or `kPa`), even if the numerical value is unfriendly!

Interrogating `ConservedProps` and `PrimitiveProps` is easily done via:

- `density(u)`: Density of state `u`.
- `momentum(u, gas)`: Momentum of state `u`
- `velocity(u, gas)`: Velocity of state `u`
- `mach_number(u, gas)`: Mach number of state `u`
- `temperature(u, gas)`: Temperature of state `u`
- `total_internal_energy_density(u, gas)`: Total (internal + kinetic) energy density of state `u`

Other properties can be directly computed:

- `pressure(u, gas)`: also offers overloads for certain special cases.
- `static_enthalpy_density(u, gas)`: Computes the static enthalpy density of `u`. Does **NOT** include kinetic energy!
- `total_enthalpy_density(u, gas)`: Computes total enthalpy density of `u`. **DOES** include kinetic energy!

"Specific" properties are related to the mass of a state, rather than its volume.

- `specific_internal_energy(u, gas)`
- `specific_static_enthalpy(u, gas)`
- `specific_total_enthalpy(u, gas)`

We can compute the change in properties across a shock wave from the shock normal $\hat n$ and shock tangent $\hat t$ via

```julia
stream_R = state_behind(free_stream, n̂, t̂)
u_R = state_behind(u_L, n̂, t̂)
```

## Tests

The test suite does some basic dimensional analysis -- the speed of sound should be a velocity, pressure should be a pressure, etc. Additionally, it tests the Rankine-Hugoniot conditions for a stable shock and verifies some algebraic properties of the Billig shockwave parametrization.

## Development Goals

- Extend this module with routines to compute interactions between shock waves
- Integrate this module into a larger project that can handle time-dependent situations and simulations.

## Related Packages

- [Euler2D.jl](https://github.com/STCE-at-RWTH/Euler2D.jl)
- [Unitful.jl](https://github.com/PainterQubits/Unitful.jl)
- [ChainRules.jl](https://github.com/JuliaDiff/ChainRules.jl)

## Contributing

Simply fork this repository and submit a pull request, should you wish.
