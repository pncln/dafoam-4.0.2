DAFoam with DASonicPimpleFoam
======================================
Super-sonic sovler for DAFoam

Mathematical Summary
-------------

  | Aspect     | DAPimpleFoam                 | DASonicPimpleFoam                               |
  |------------|------------------------------|-------------------------------------------------|
  | Momentum   | ∂U/∂t + ∇·(UU) = -∇p/ρ + ∇·τ | ∂(ρU)/∂t + ∇·(ρUU) = -∇p + ∇·τ                  |
  | Continuity | ∇·U = 0                      | ∂ρ/∂t + ∇·(ρU) = 0                              |
  | Energy     | ∂T/∂t + ∇·(UT) = ∇·(α∇T)     | ∂(ρh)/∂t + ∇·(ρUh) + ∇·(pU) = ∇·(k∇T) + viscous |
  | Pressure   | ∇²p = ∇·(∇·UU)               | ∂(ψp)/∂t + ∇·(ρU) = ∇·(ρ∇p/A)                   |
  | State      | ρ = const                    | p = ρRT, ρ = f(p,T)                             |

Documentation
-------------

cd $DAFOAM_ROOT_PATH/repos && \
wget https://github.com/pncln/dafoam-4.0.2/archive/refs/tags/test.tar.gz -O dafoam.tar.gz && \
tar -xvf dafoam.tar.gz && cd dafoam-* && \
. $DAFOAM_ROOT_PATH/loadDAFoam.sh && \
./Allmake

Refer to https://dafoam.github.io for installation, documentation, and tutorials.
