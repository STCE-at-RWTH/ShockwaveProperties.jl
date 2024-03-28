# ShockwaveProperties.jl

---

This package provides routines for computing property changes across shock waves in perfect gases. The ultimate goal is to support a solver that can be used to analyze complex interactions of shock waves in dynamic systems.

# Usage

The `PrimitiveState` and `ConservedState` types allow for conversion between specifying gas properties as $\rho, \vec M, T$ and $\rho, \rho \vec v,\rho E$. These types will take on the SI standard units by default. 

```julia
free_stream = PrimitiveState(1.225u"kg/m^3", [2.0, 0.0], 300u"K")
u = ConservedState(free_stream, DRY_AIR)
```

We can compute the change in properties across a shock wave from the shock normal $\hat n$ and shock tangent $\hat t$ via

```julia
stream_R = state_behind(free_stream, n̂, t̂)
u_R = state_behind(u_L, n̂, t̂)
```

## Tests
The test suite does some basic dimensional analysis -- the speed of sound should be a velocity, pressure should be a pressure, etc. Additionally, it tests the Rankine-Hugoniot conditions for a stable shock and verifies some algebraic properties of the Billig shockwave parametrization.

# Development Goals

- Extend this module with a solver for the Euler equations to accurately compute flow properties behind shock waves.
- Extend this module with routines to compute interactions between shock waves
- Integrate this module into a larger project that can handle time-dependent situations and simulations.


