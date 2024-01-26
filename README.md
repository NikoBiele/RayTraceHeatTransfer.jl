# RayTraceHeatTransfer.jl

[![Build Status](https://travis-ci.com/NikoBiele/RayTraceHeatTransfer.jl.svg?branch=main)](https://travis-ci.com/NikoBiele/RayTraceHeatTransfer.jl)
[![Coverage](https://codecov.io/gh/NikoBiele/RayTraceHeatTransfer.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/NikoBiele/RayTraceHeatTransfer.jl)

## Description

This repository can be used for radiation heat transfer calculations in an enclosure including a participating medium.
It contains a number of functions which collectively enables the user to ray trace a user defined geometry.
The result of the ray tracing are four 'exchange factor' matrices which together describe how the enclosure is radiatively connected.
Using the exchange factor matrices it is possible to quickly perform a heat transfer calculation on the entire enclosure, which would otherwise be computationally expensive to ray trace.
This package is limited to a uniformly distributed participating medium.
It is also limited to 2D enclosures (or more accurately: specular/mirrorlike front and back, since the sampled distributions are 3D).

## Features

- Define a custom geometry
- Ray trace the geometry rapidly
- 'Save' ray tracing result as Exchange Factor matrices (also as CSV-files)
- Quickly calculate heat transfer in the geometry using the Exchange Factors
- Avoid ray tracing the same geometry multiple times

## Installation

For now this repository must be downloaded from GitHub since it is not yet registered.

## Usage

### Generate geometry

'''julia
yLayersHeight = [0.0, 1.0];
    N_subs = length(yLayersHeight)-1; # number of sub-enclosures
    xLayersWidth = zeros(2, length(yLayersHeight));
    xLayersWidth[:,1] = [0.0, 1.0];
    xLayersWidth[:,2] = [0.0, 1.0];
'''
