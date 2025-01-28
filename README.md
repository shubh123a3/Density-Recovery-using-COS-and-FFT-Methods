# Density Recovery Methods: COS and FFT Implementations

A comprehensive implementation of numerical methods for recovering probability density functions from characteristic functions, featuring both COS (Cosine expansion) and FFT (Fast Fourier Transform) methods with applications to Normal and LogNormal distributions.

## Overview

This project provides tools and implementations for recovering probability density functions using two distinct numerical approaches. It includes an interactive web application for visualization and comparison of different methods and parameters.

## Table of Contents
1. [Mathematical Foundation](#mathematical-foundation)
2. [Implementation Details](#implementation-details)
3. [Interactive Application](#interactive-application)
4. [Technical Details](#technical-details)
5. [Usage Examples](#usage-examples)
6. [Performance Analysis](#performance-analysis)
7. [Installation and Setup](#installation-and-setup)

## Mathematical Foundation

### Characteristic Functions
The characteristic function is the foundation of both recovery methods:

Ctrl+K to generate a command
@LogNormal density recovery using the COS method.ipynb @Normal density recovery using FFT.ipynb @Normal density recovery using the COS method.ipynb  create a 3-0 word repo discription of projeact and then create very very detailed readme file of my project everything about it with theory

Here's a 3-word description of your project:
Density Recovery Methods
And here's a detailed README.md for your project:
math
φ_X(u) = E[e^{iuX}] = ∫{-∞}^{∞} e^{iux} f_X(x) dx
math
f_X(x) ≈ ∑{k=0}^{N-1}' A_k cos(kπ (x-a)/(b-a))
math
A_k ≈ \frac{2}{b-a} Re{φ_X(\frac{kπ}{b-a}) e^{-ikaπ/(b-a)}}
math
f_X(x) = \frac{1}{2π} ∫{-∞}^{∞} e^{-iux} φ_X(u) du
python
# COS Method
def COSDensity(cf, x, N, a, b):
"""
Recovers density using COS method
Args:
cf: Characteristic function
x: Points for density evaluation
N: Number of expansion terms
a,b: Truncation range
"""
FFT Method
def FFTDensity(cf, x, N):
"""
Recovers density using FFT method
Args:
cf: Characteristic function
x: Points for density evaluation
N: Number of FFT points
"""
bash
pip install numpy scipy matplotlib streamlit plotly
bash
streamlit run app.py
This README provides a comprehensive overview of your project, including the theoretical background, implementation details, usage instructions, and performance characteristics. It's structured to be both informative for users and helpful for developers who might want to understand or contribute to the project.
give me this in one readme file in toltal  in vbery detailed .

Here's a complete, detailed README.md file for your project:
math
φ_X(u) = E[e^{iuX}] = ∫{-∞}^{∞} e^{iux} f_X(x) dx

#### Normal Distribution
- Characteristic function: `φ_X(u) = exp(iuμ - ½σ²u²)`
- PDF: `f_X(x) = (1/σ√(2π)) exp(-(x-μ)²/2σ²)`

#### LogNormal Distribution
- Related to normal distribution through: Y = exp(X)
- PDF: `f_Y(y) = (1/yσ√(2π)) exp(-(ln y-μ)²/2σ²)`

### COS Method
The COS method uses a Fourier-cosine series expansion:

Ctrl+K to generate a command
@LogNormal density recovery using the COS method.ipynb @Normal density recovery using FFT.ipynb @Normal density recovery using the COS method.ipynb  create a 3-0 word repo discription of projeact and then create very very detailed readme file of my project everything about it with theory

Here's a 3-word description of your project:
Density Recovery Methods
And here's a detailed README.md for your project:
math
φ_X(u) = E[e^{iuX}] = ∫{-∞}^{∞} e^{iux} f_X(x) dx
math
f_X(x) ≈ ∑{k=0}^{N-1}' A_k cos(kπ (x-a)/(b-a))
math
A_k ≈ \frac{2}{b-a} Re{φ_X(\frac{kπ}{b-a}) e^{-ikaπ/(b-a)}}
math
f_X(x) = \frac{1}{2π} ∫{-∞}^{∞} e^{-iux} φ_X(u) du
python
# COS Method
def COSDensity(cf, x, N, a, b):
"""
Recovers density using COS method
Args:
cf: Characteristic function
x: Points for density evaluation
N: Number of expansion terms
a,b: Truncation range
"""
FFT Method
def FFTDensity(cf, x, N):
"""
Recovers density using FFT method
Args:
cf: Characteristic function
x: Points for density evaluation
N: Number of FFT points
"""
bash
pip install numpy scipy matplotlib streamlit plotly
bash
streamlit run app.py
This README provides a comprehensive overview of your project, including the theoretical background, implementation details, usage instructions, and performance characteristics. It's structured to be both informative for users and helpful for developers who might want to understand or contribute to the project.
give me this in one readme file in toltal  in vbery detailed .

Here's a complete, detailed README.md file for your project:
math
φ_X(u) = E[e^{iuX}] = ∫{-∞}^{∞} e^{iux} f_X(x) dx
math
f_X(x) ≈ ∑{k=0}^{N-1}' A_k cos(kπ (x-a)/(b-a))

Key components:
- Truncation range [a,b]
- Expansion terms N
- Fourier coefficients:
  ```math
  A_k ≈ (2/(b-a)) Re{φ_X(kπ/(b-a)) e^{-ikaπ/(b-a)}}
  ```

### FFT Method
Based on the inverse Fourier transform:
math
f_X(x) = (1/2π) ∫{-∞}^{∞} e^{-iux} φ_X(u) du
## Implementation Details

### Core Functions

python
def COSDensity(cf, x, N, a, b):
"""
Recovers density using COS method.
Parameters:
cf (function): Characteristic function
x (array): Points for density evaluation
N (int): Number of expansion terms
a, b (float): Truncation range bounds
Returns:
array: Recovered density values
"""
i = 1j
k = np.arange(N)
u = k np.pi / (b - a)
f_k = 2.0 / (b - a) np.real(cf(u) np.exp(-i u a))
f_k[0] = 0.5
f_x = np.dot(f_k, np.cos(np.outer(u, x - a)))
return f_x
def FFTDensity(cf, x, N):
"""
Recovers density using FFT method.
Parameters:
cf (function): Characteristic function
x (array): Points for density evaluation
N (int): Number of FFT points
Returns:
array: Recovered density values
"""
i = 1j
u_max = 20.0
du = u_max / N
u = np.linspace(0, N-1, N) du
b = np.min(x)
dx = 2.0 np.pi / (N du)
x_i = b + np.linspace(0, N-1, N) dx
phi = np.exp(-i b u) cf(u)
gamma_1 = np.exp(-i x_i u[0]) cf(u[0])
gamma_2 = np.exp(-i x_i u[-1]) cf(u[-1])
phi_boundary = 0.5 (gamma_1 + gamma_2)
f_xi = du / np.pi np.real(np.fft.fft(phi) - phi_boundary)
return np.interp(x, x_i, f_xi)


## Interactive Application

The project includes a Streamlit web application (`app.py`) with features:

### User Interface Components
- Distribution selection (Normal/LogNormal)
- Method selection (COS/FFT)
- Parameter controls:
  - Mean (μ)
  - Standard deviation (σ)
  - Truncation range [a,b]
  - Number of terms (N)

### Visualization Features
- Real-time density plotting
- Comparison with exact density
- Multiple N-value comparisons
- Interactive parameter adjustment

## Technical Details

### Dependencies
- Python 3.x
- NumPy: Numerical computations
- SciPy: Statistical functions
- Matplotlib: Basic plotting
- Plotly: Interactive visualization
- Streamlit: Web interface

### Key Features
1. **Adaptive Computation**
   - Dynamic N selection
   - Automatic range adjustment
   - Error handling

2. **Visualization Options**
   - Multiple plot types
   - Error analysis
   - Convergence visualization

3. **Performance Optimizations**
   - Vectorized operations
   - Efficient FFT implementation
   - Memory management

## Usage Examples

### Basic Usage

Normal distribution example
mu, sigma = 0.0, 1.0
cf = lambda u: np.exp(1j mu u - 0.5 sigma2 u2)
x = np.linspace(-5, 5, 1000)
COS method
density_cos = COSDensity(cf, x, N=64, a=-10, b=10)
FFT method
density_fft = FFTDensity(cf, x, N=1024)
### Web Application

Ctrl+K to generate a command
@LogNormal density recovery using the COS method.ipynb @Normal density recovery using FFT.ipynb @Normal density recovery using the COS method.ipynb  create a 3-0 word repo discription of projeact and then create very very detailed readme file of my project everything about it with theory

Here's a 3-word description of your project:
Density Recovery Methods
And here's a detailed README.md for your project:
math
φ_X(u) = E[e^{iuX}] = ∫{-∞}^{∞} e^{iux} f_X(x) dx
math
f_X(x) ≈ ∑{k=0}^{N-1}' A_k cos(kπ (x-a)/(b-a))
math
A_k ≈ \frac{2}{b-a} Re{φ_X(\frac{kπ}{b-a}) e^{-ikaπ/(b-a)}}
math
f_X(x) = \frac{1}{2π} ∫{-∞}^{∞} e^{-iux} φ_X(u) du
python
# COS Method
def COSDensity(cf, x, N, a, b):
"""
Recovers density using COS method
Args:
cf: Characteristic function
x: Points for density evaluation
N: Number of expansion terms
a,b: Truncation range
"""
FFT Method
def FFTDensity(cf, x, N):
"""
Recovers density using FFT method
Args:
cf: Characteristic function
x: Points for density evaluation
N: Number of FFT points
"""
bash
pip install numpy scipy matplotlib streamlit plotly
bash
streamlit run app.py
This README provides a comprehensive overview of your project, including the theoretical background, implementation details, usage instructions, and performance characteristics. It's structured to be both informative for users and helpful for developers who might want to understand or contribute to the project.
give me this in one readme file in toltal  in vbery detailed .

Here's a complete, detailed README.md file for your project:
math
φ_X(u) = E[e^{iuX}] = ∫{-∞}^{∞} e^{iux} f_X(x) dx
math
f_X(x) ≈ ∑{k=0}^{N-1}' A_k cos(kπ (x-a)/(b-a))
math
f_X(x) = (1/2π) ∫{-∞}^{∞} e^{-iux} φ_X(u) du
python
def COSDensity(cf, x, N, a, b):
"""
Recovers density using COS method.
Parameters:
cf (function): Characteristic function
x (array): Points for density evaluation
N (int): Number of expansion terms
a, b (float): Truncation range bounds
Returns:
array: Recovered density values
"""
i = 1j
k = np.arange(N)
u = k np.pi / (b - a)
f_k = 2.0 / (b - a) np.real(cf(u) np.exp(-i u a))
f_k[0] = 0.5
f_x = np.dot(f_k, np.cos(np.outer(u, x - a)))
return f_x
def FFTDensity(cf, x, N):
"""
Recovers density using FFT method.
Parameters:
cf (function): Characteristic function
x (array): Points for density evaluation
N (int): Number of FFT points
Returns:
array: Recovered density values
"""
i = 1j
u_max = 20.0
du = u_max / N
u = np.linspace(0, N-1, N) du
b = np.min(x)
dx = 2.0 np.pi / (N du)
x_i = b + np.linspace(0, N-1, N) dx
phi = np.exp(-i b u) cf(u)
gamma_1 = np.exp(-i x_i u[0]) cf(u[0])
gamma_2 = np.exp(-i x_i u[-1]) cf(u[-1])
phi_boundary = 0.5 (gamma_1 + gamma_2)
f_xi = du / np.pi np.real(np.fft.fft(phi) - phi_boundary)
return np.interp(x, x_i, f_xi)
python
Normal distribution example
mu, sigma = 0.0, 1.0
cf = lambda u: np.exp(1j mu u - 0.5 sigma2 u2)
x = np.linspace(-5, 5, 1000)
COS method
density_cos = COSDensity(cf, x, N=64, a=-10, b=10)
FFT method
density_fft = FFTDensity(cf, x, N=1024)
bash
streamlit run app.py


## Performance Analysis

### COS Method
Advantages:
- Exponential convergence for smooth densities
- Efficient with low N
- Excellent tail accuracy

Limitations:
- Requires truncation range selection
- Potential oscillations with very low N

### FFT Method
Advantages:
- No truncation range needed
- Natural for periodic functions
- Stable implementation

Limitations:
- Higher N requirements
- Boundary artifacts
- Numerical sensitivity

## Installation and Setup

1. Clone the repository:
bash
git clone https://github.com/yourusername/density-recovery-methods.git

2. Install dependencies:
   bash
pip install -r requirements.txt

3. Run tests:
   python -m pytest tests/
4. Launch web application:
5. bash
streamlit run app.py


## Contributing

1. Fork the repository
2. Create feature branch
3. Commit changes
4. Push to branch
5. Create Pull Request

## License

This project is licensed under the MIT License.

## Citation

If you use this code in your research, please cite:
****
@software{density_recovery_methods,
author = shubh shrishrimal,
title = {Density Recovery Methods: COS and FFT Implementations},
year = {2024},


## Contact

For questions and support, please open an issue or contact [your-email@example.com].
