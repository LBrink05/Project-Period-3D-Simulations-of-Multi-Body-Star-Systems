# 🌟 3D Simulations of Multi-Body Star Systems

A high-performance computational framework for simulating and analyzing gravitational dynamics in multi-body stellar systems using advanced numerical methods and Lyapunov stability analysis.

---

## 📋 Overview

This project provides a comprehensive toolkit for investigating stellar dynamics through numerical simulation of N-body gravitational systems. By combining numerical integration methods with real-time visualization capabilities, researchers can explore chaotic orbital mechanics, close encounters, and long-term stability patterns in multi-star configurations.

### Key Features

- 🚀 **Numerical Methods**: Compare symplectic, adaptive, arbitrary-precision, and series-based integrators
- 📊 **Lyapunov Stability Analysis**: Full spectrum computation for orbital predictability assessment
- 🎨 **Multi-Platform Visualization**: Switch between Matplotlib, Unity, and Manim renderers
- ⚡ **C++ Core**: C++ engine optimized for computational efficiency
- 🖥️ **Python Interface**: Python-based GUI built with CustomTkinter
- 📈 **Comprehensive Error Metrics**: Phase-space trajectory accuracy and Hamiltonian conservation tracking

---

## 🔬 Scientific Background

### Stellar Configurations

The simulation framework supports various test configurations for benchmarking and analysis:

#### Pythagorean Three-Body Problem
Three bodies with mass ratios **3:4:5** arranged in a right triangle configuration. This highly chaotic system serves as an ideal benchmark for comparing numerical methods due to its sensitivity to initial conditions.

```
         Mass 5 (●)
           ╱│
          ╱ │
         ╱  │ 4 units
        ╱   │
       ╱    │
      ╱_____│
  Mass 3   Mass 4
  (●)  3 units  (●)
```

#### Lagrange-Euler Configuration
Three equal-mass bodies positioned equidistantly on a circle with tangential velocities, producing stable circular orbits around the center of mass—perfect for validating numerical accuracy against analytical solutions.

### Numerical Methods Compared

| Category | Method | Characteristics |
|----------|--------|-----------------|
| **Arbitrary-Precision** | BRUTUS | Minimized round-off error for reference trajectories |
| **Symplectic** | Leapfrog/Verlet | Preserves Hamiltonian structure and phase-space volume |
| **Adaptive** | Adaptive RK/Verlet | Dynamic timestep adjustment based on error estimates |
| **Semi-Analytical** | Adomian Decomposition | Series-based solution approximations |
| **Fixed-Step** | Standard RK4 | Non-symplectic classical method |

---

## 📐 Mathematical Framework

### Trajectory Error Metrics

**Phase-Space Error** measures deviation from reference trajectories:

$$E_{\%}(t) = \frac{\lVert\vec{x}_{\text{reference}}(t) - \vec{x}_{\text{approximate}}(t)\rVert}{\lVert\vec{x}_{\text{reference}}(t)\rVert} \times 100\\%$$

Where phase-space vector: $\vec{x}(t) = (q_1, q_2, q_3, \ldots, p_1, p_2, p_3, \ldots)$ for all bodies

**Root-Mean-Square Error** aggregates trajectory accuracy:

$$E_{\%\text{RMS}} = \sqrt{\frac{1}{T} \sum_{t=0}^{T} E_{\%}(t)^2 \Delta t}$$

**Hamiltonian Conservation Error** quantifies energy drift:

$$E_{H\%} = \frac{\lvert H_0 - H(t)\rvert}{H(t)} \times 100\\%$$

The simulation computes the full Lyapunov spectrum `{λ₁, λ₂, ..., λ₆ₙ}` to characterize orbital stability:

1. **Initialize** perturbation vectors as identity matrix in 6N-dimensional phase space
2. **Evolve** perturbations alongside trajectories using the Jacobian matrix
3. **Renormalize** via QR decomposition at regular intervals to prevent alignment
4. **Accumulate** logarithmic growth rates across all renormalization steps

**Maximum Lyapunov Exponent** (λₘₐₓ):
- λₘₐₓ > 0: Chaotic, exponentially diverging trajectories
- λₘₐₓ ≈ 0: Regular, quasi-periodic orbits
- λₘₐₓ < 0: Stable, converging dynamics

**Lyapunov Time**: `t_L = 1/λₘₐₓ` quantifies predictability horizon

---

## 🛠️ Technology Stack

### Core Engine
- **C++**: High-performance numerical integration and physics calculations
- Optimized for computational efficiency in long-term simulations

### User Interface & Visualization
- **Python**: Main interface logic and data analysis
- **CustomTkinter**: Modern, responsive GUI framework
- **Matplotlib**: 2D trajectory plotting and error visualization
- **Unity**: Real-time 3D interactive visualization
- **Manim**: Publication-quality animation generation

---

## 🚀 Getting Started

### Prerequisites

```bash
# C++ Compiler (GCC 9.0+ or Clang 10.0+)
# Python 3.8+
# Required Python packages:
pip install customtkinter matplotlib numpy scipy
# or use included env with after Installation:
source env/bin/activate
```

### Installation

```bash
git clone https://github.com/UM-MSP/3D-Simulations-of-Multi-Body-Star-Systems
cd 3D-Simulations-of-Multi-Body-Star-Systems
```

### Quick Start

```bash
# Compile C++ engine
make build

# Launch GUI
python main.py
```

---

## 📊 Evaluation Criteria

Numerical methods are benchmarked across multiple dimensions:

| Criterion | Description |
|-----------|-------------|
| **Computation Time** | Time under fixed resource constraints |
| **E_%RMS (Distance)** | Accuracy for close encounters vs. distant interactions |
| **E_%RMS (Mass)** | Performance across different mass ratios |
| **E_%RMS (Velocity)** | Precision for fast-moving vs. slow-moving bodies |
| **E_H%** | Hamiltonian conservation (symplectic structure preservation) |

---

## 📖 Usage Examples

### Running a Pythagorean Configuration

```python
from stellar_sim import Simulation, Configuration

# Initialize Pythagorean three-body problem
config = Configuration.pythagorean(masses=[3, 4, 5])
sim = Simulation(config, method='adaptive_verlet')

# Run simulation
results = sim.run(t_max=100, collision_threshold=0.1)

# Analyze stability
lyapunov = sim.compute_lyapunov_spectrum()
print(f"Max Lyapunov Exponent: {lyapunov.max():.4f}")
print(f"Lyapunov Time: {1/lyapunov.max():.2f} time units")
```

### Comparing Methods

```python
methods = ['brutus', 'verlet', 'adaptive_rk', 'rk4']
results = {}

for method in methods:
    sim = Simulation(config, method=method)
    results[method] = sim.benchmark()

# Plot comparison
plot_method_comparison(results)
```

---

## 🎨 Visualization Options

Switch between rendering engines based on your needs:

- **Matplotlib**: Quick 2D plots, error graphs, phase-space diagrams
- **Unity**: Interactive 3D exploration with real-time camera control
- **Manim**: Smooth, publication-ready animations for presentations

---

## 📚 Project Paper

- Read our Document on https://www.overleaf.com/project/68fc924425dbd906ca8ddfc4
