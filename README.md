# Phyzo

Local PDE computation tool written in C++20. Windows x64, CLI and Win32 GUI.

Scalar FEM, a mixed P2/P1 formulation for incompressible flow in 2D/3D, and a cosine Galerkin method. Steady and transient solves, mesh, quadrature and time-step sensitivity studies. Run archives store the model, parameters, mesh, diagnostics, source code and executable for replay.

## Build

Requires MSVC with C++20 support, Windows SDK and CMake 3.24+. Install dependencies separately; the builder does not download them.

Place the following dependencies under `%PHYZO_SDK%\deps`:

- `json.hpp`: nlohmann/json 3.12.0.
- `eigen-3.4.0\Eigen` and `boost_1_86_0\boost`: Eigen 3.4.0 and Boost 1.86.0.
- `sundials\include` and `sundials\lib`: SUNDIALS 7.4.0, static x64 Release build with CVODE, serial NVector, dense SUNMatrix and dense SUNLinearSolver.
- `miniz`: miniz 3.0.2 sources; `tinyxml2`: tinyxml2 11.0.0 sources.

From the repository root in an x64 Native Tools Command Prompt:

```bat
set PHYZO_SDK=C:\Phyzo-native
cl /std:c++20 /EHsc /utf-8 build.cpp /Fe:BuildPhyzo.exe
BuildPhyzo.exe
```

Output: `%PHYZO_SDK%\bin\Phyzo.exe`. Build files are generated inside the SDK directory.

## Usage

```bat
Phyzo.exe examples
Phyzo.exe run heat3
Phyzo.exe stat nseq2
Phyzo.exe study heat3 --kind mesh --levels 3
Phyzo.exe selftest --count 1000 --seed 20261115 -o tests.json
Phyzo.exe --help
```

Pass a JSON model path instead of a built-in example name to use your own model. Running without arguments opens the GUI. `PHYZO_DATA` sets the output directory; `verify` checks an archive and `replay` repeats the computation.
