# RayTraceHeatTransfer.jl Documentation

Radiative heat transfer in participating media

## Description

This Julia package can be used for radiation heat transfer calculations in an enclosure with partly or fully reflecting walls containing an absorbing-emitting-scattering participating medium.
This phenomenon is governed by the **Radiative Transfer Equation** (RTE), which is a integro-differential equation:

```math
\frac{\partial I_{\lambda}(S,\Omega)}{\partial S} = \kappa_{\lambda} I_{\lambda \mathrm{b}}(S) - \kappa_{\lambda} I_{\lambda}(S,\Omega) - \sigma_{\mathrm{s},\lambda} I_{\lambda}(S,\Omega) + \frac{\sigma_{\mathrm{s},\lambda}}{4 \pi} \int_{\Omega_i=4 \pi} I_{\lambda}(S,\Omega_i) \Phi_{\lambda}(\Omega_i,\Omega) d\Omega_i
```

Analytical solution of the RTE is only possible in simple cases.
This repository contains a number of types and functions which collectively enables the user to build, mesh and ray trace a user defined 2D geometry to solve the RTE.
The result of the ray tracing is four 'Exchange Factor' matrices which together describe how the enclosure is radiatively connected.
Using the exchange factor matrices it is possible to quickly perform a heat transfer calculation on the entire enclosure, which would otherwise be computationally expensive to ray trace.

### Features

- Define a custom 2D geometry in a modular and interactive way.
- Mesh the user defined geometry.
- Ray trace the meshed geometry in parallel.
- Save and load ray tracing results as Exchange Factor matrices in CSV-format.
- Quickly calculate heat transfer in the geometry using the Exchange Factors.
- Avoid ray tracing the same geometry multiple times when changing boundary conditions.

### Limitations

This package is limited to:
- A uniformly distributed gray participating medium.
- Diffuse surface emission/reflection.