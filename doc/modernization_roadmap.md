# XNet Modernization Roadmap

These phases are sequential — each assumes the prior phase is substantially
complete. Within a phase, work items can be parallelized.

## Phase 0: Establish Test Infrastructure (prerequisite for everything)
- [ ] Scrap existing CI (`makefile_tests.yml`, `test_xnet.sh` as CI driver)
- [ ] Add test-drive as git submodule
- [ ] Create first unit tests (Layer 0-1: types, constants, util, fd)
- [ ] Build new GitHub Actions workflow with Tier 1 tests
- [ ] Add `make test_unit` target to Makefile
- [ ] Move `xnet_integrate_bd.F90` to `archive/`
- [ ] Fix `make clean` to remove executables

## Phase 1: Expand Test Coverage
- [ ] Unit tests for Layer 2-3 (data, conditions, screening, flux, ffn, nnu)
- [ ] Tier 2 physics-informed tests (conservation laws, NSE consistency,
      Jacobian accuracy, convergence order)
- [ ] Create new Tier 3 regression tests with proper tolerance checking
- [ ] Retire old regression tests as replacements are established
- [ ] Establish tolerance standards

## Phase 2: Code Restructuring
- [ ] Begin directory restructuring (`source/` -> `src/` with subdirectories)
- [ ] Extract `xnet_jacobian_common.F90` (requires Jacobian unit tests first)
- [ ] Separate `extern/` from `tools/`
- [ ] Documentation: bring module headers to target standard (FORD markup)

## Phase 3: Build Modernization
- [ ] Introduce CMake (alongside Make, not replacing it yet)
- [ ] Composable compiler/site configs (`cmake/compilers/`, `cmake/sites/`)
- [ ] Out-of-source build support
- [ ] Expand CI compiler matrix (Intel ifx, etc.)

## Phase 4: State Object Migration (deferred)
- [ ] Design state types (start from leaves: conditions, abundances)
- [ ] Incremental migration, module by module
- [ ] Requires comprehensive test coverage from Phases 0-1 as safety net
