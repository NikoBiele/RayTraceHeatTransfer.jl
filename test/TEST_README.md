# RayTraceHeatTransfer.jl Test Suite

## Philosophy

This test suite follows the principle of **validation against reference solutions** rather than exact reproduction of specific numerical values. Since the package relies on Monte Carlo ray tracing, results will have statistical noise. The tests are designed to:

1. **Compare with analytical/reference solutions** - Test against established benchmarks (Narayanaswamy view factors, EES solutions, Crosbie & Schrenker)
2. **Test invariance properties** - Solutions should be invariant under rotations and similar for different mesh refinements
3. **Verify physical constraints** - Energy conservation, reciprocity, temperature bounds
4. **Check consistency** - Grey vs spectral with black surfaces should match

## Test Structure

```
test/
├── runtests.jl                    # Main test runner
├── test_3d_viewfactors.jl         # 3D view factor calculations
├── test_3d_heat_transfer.jl       # 3D transparent surface heat transfer
├── test_2d_grey.jl                # 2D grey participating media
├── test_2d_spectral.jl            # 2D spectral participating media
├── test_spectral_consistency.jl   # Grey vs spectral consistency
└── TEST_README.md                 # This file
```

## Test Categories

### 1. 3D View Factors (`test_3d_viewfactors.jl`)

**Purpose**: Validate view factor calculations against analytical solutions

**Tests**:
- Narayanaswamy examples (7 test cases with various geometries)
- Cube aligned with coordinate system (comparison with EES reference)
- Rotated cube (invariance under rotation)
- Mesh refinement independence

**Tolerances**:
- View factors: `1e-5` (absolute)
- Reciprocity: `1e-10` (relative)

### 2. 3D Heat Transfer (`test_3d_heat_transfer.jl`)

**Purpose**: Validate heat transfer solutions for transparent surfaces

**Tests**:
- Isothermal enclosure (all walls same temperature)
- Two hot walls at different temperatures
- Energy conservation with specified heat flux
- Rotational invariance of solutions
- Grey surface properties (non-black surfaces)

**Tolerances**:
- Temperature: `5.0` K (absolute)
- Energy balance: `1e-4` W (absolute)

### 3. 2D Grey Participating Media (`test_2d_grey.jl`)

**Purpose**: Validate 2D ray tracing with absorbing/emitting media

**Tests**:
- Aligned square with single hot wall (4 configurations)
- Rotated square at 45°
- Multiple rotation angles (invariance testing)
- Mesh refinement convergence
- Energy conservation

**Reference**: Crosbie & Schrenker (1984) analytical solution

**Tolerances**:
- vs analytical: `0.05` (5% relative error, to avoid excessive sampling)

### 4. 2D Spectral Participating Media (`test_2d_spectral.jl`)

**Purpose**: Validate spectral radiation modeling

**Tests**:
- Spectral uniform vs grey (should match with black walls)
- Exchange factor vs direct method consistency
- Spectral variable (wavelength-dependent properties)
- Spectral energy conservation
- Selective surface properties
- Varying number of spectral bins

**Tolerances**:
- Spectral vs grey: `0.05` (5% relative)

### 5. Spectral Consistency (`test_spectral_consistency.jl`)

**Purpose**: Verify grey and spectral implementations are consistent

**Tests**:
- 3D spectral vs grey comparison
- Spectral bin integration (sum over bins = total)
- Selective vs black surface comparison
- Planck function integration check
- Spectral mode auto-detection
- Energy conservation per bin and total

**Tolerances**:
- Consistency: `0.05` (5% relative)

## Running the Tests

### Run all tests (prints progress to REPL):
```julia
using Pkg
Pkg.test("RayTraceHeatTransfer")
```

### Run specific test file:
```julia
using Test, RayTraceHeatTransfer
include("test/test_3d_viewfactors.jl")
```

## Understanding Test Failures

### Statistical Noise
Since ray tracing is Monte Carlo-based, individual runs may occasionally fail due to statistical fluctuations. If a test fails:
1. Run it again to see if it's reproducible
2. Check if the error is small (near the tolerance boundary)
3. Consider increasing `N_rays` for more accuracy

### Tolerance Philosophy
- **Tight tolerances** (`1e-10`): Mathematical properties (reciprocity, self-view = 0)
- **Moderate tolerances** (`1e-5` to `1e-4`): Numerical accuracy (view factors, energy conservation)
- **Loose tolerances** (`0.05` and `5` K): Statistical comparisons (MC methods, rotational invariance, temperature)

### Common Issues

**High energy error**: 
- May need more rays
- Check boundary conditions are consistent
- Verify spectral bins integrate properly

**Rotational invariance failure**:
- Solution statistics should be invariant, not individual values
- Check that mesh refinement is sufficient
- Geometric transformations may introduce small numerical differences

**Grey vs spectral mismatch**:
- Verify black surface properties (ε = 1 in all bins)
- Check wavelength range covers significant Planck function range
- Ensure sufficient spectral bins for integration

## Adding New Tests

When adding tests, follow these guidelines:

1. **Name descriptively**: Test names should clearly indicate what's being tested
2. **Use tolerances**: Always use appropriate tolerance for `@test` (avoid exact equality)
3. **Check physical bounds**: Verify temperatures, heat fluxes are physically reasonable
4. **Test edge cases**: Include degenerate or extreme geometries
5. **Document references**: Cite analytical solutions or benchmarks
6. **Keep tests fast**: Use smaller meshes/fewer rays for routine testing, add slow tests separately

## Performance Benchmarks

Typical run times on a modern workstation (approximate):

- `test_3d_viewfactors.jl`: ~10 seconds
- `test_3d_heat_transfer.jl`: ~15 seconds  
- `test_2d_grey.jl`: ~30 seconds
- `test_2d_spectral.jl`: ~45 seconds
- `test_spectral_consistency.jl`: ~20 seconds

**Total**: ~2 minutes

For faster testing during development, reduce:
- `N_rays_total` (fewer rays, more noise)
- `Ndim` (coarser mesh, less accuracy)
- `n_bins` (fewer spectral bins)

## References

1. **Narayanaswamy, A.** (2015). "An Analytic Expression for Radiation View Factor Between Two Arbitrarily Oriented Planar Polygons." *International Journal of Heat and Mass Transfer*.

2. **Crosbie, A. L., & Schrenker, R. G.** (1984). "Radiative Transfer in a Two-Dimensional Rectangular Medium Exposed to Diffuse Radiation" *Journal of Quantitative Spectroscopy and Radiative Transfer*.

3. **EES, Klein, S. & Nellis G.** (Engineering Equation Solver) - Reference view factor calculations for cube geometry.

## Continuous Integration

The test suite is designed to run in CI environments:
- All tests should complete in < 5 minutes
- Tests use reasonable but not excessive resources
- Deterministic tests (view factors, energy conservation) should always pass
- Statistical tests may rarely fail due to random fluctuations (acceptable)

## Questions or Issues?

If you encounter persistent test failures or have questions about the test suite, please open an issue on the GitHub repository with:
1. Which test(s) failed
2. Error messages
3. Julia version and system information
4. Whether the failure is reproducible