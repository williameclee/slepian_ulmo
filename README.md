# slepian_ulmo

A fork of `csdms-contrib` repositories that are dedicated to constructing Slepian functions over ocean basins.

## Functionalities

Functions in this repository may call functions from the `slepian_alpha`, `slepian_bravo`, or `slepian_delta` packages. Please ensure they are installed before running the functions in this repository.

### Geographic domains

The following ocean basins are supported:

- `oceans`: All oceans, excluding the Arctic Ocean at the moment
- `pacific`: The Pacific Ocean
- `npacific`: The North Pacific Ocean
- `spacific`: The South Pacific Ocean
- `atlantic`: The Atlantic Ocean
- `natlantic`: The North Atlantic Ocean
- `satlantic`: The South Atlantic Ocean
- `indian`: The Indian Ocean
- `arctic`: The Arctic Ocean, which is by default rotated to the equator

The boundaries of these ocean basins are given by the International Hydrographic Organisation (IHO)'s *Limits of Oceans and Seas*.
Differing from the `slepian_alpha` package, the coastline data in this package is from GSHHG (A Global Self-consistent, Hierarchical, High-resolution Geography Database).

Along with the ocean basins, the following geographic domains are also supported:

- `earthquakes`: A mask of coastal megathrust earthquakes.

A new class is introduced as an interface to geographic domains:

- `GeoDomain`: A class that supports a few new methods that make fetching vertices and defining file names easier.

### Modifications to functions from other packages

To support the class `GeoDomain`, the following functions have been modified:

- `kernelcp` (from `slepian_alpha`) -> `kernelcp_new`
- `glmalpha` (from `slepian_alpha`) -> `glmalpha_new`

The ultimate goal is to drop the `_new` suffixes and replace the original functions with the modified ones. Note that these functions are only tested for geographic domains, and may not work as expected for, e.g. circular caps.

---
Last modified by:
- [En-Chi Lee](mailto:williameclee@arizona.edu), 2024/08/09