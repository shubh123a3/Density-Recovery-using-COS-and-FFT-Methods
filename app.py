import streamlit as st
import numpy as np
import plotly.graph_objects as go
import scipy.stats as stats

def COSDensity(cf, x, N, a, b):
    i = 1j
    k = np.arange(N)
    u = k * np.pi / (b - a)
    f_k = 2.0 / (b - a) * np.real(cf(u) * np.exp(-i * u * a))
    f_k[0] *= 0.5
    f_x = np.dot(f_k, np.cos(np.outer(u, x - a)))
    return f_x

def FFTDensity(cf, x, N):
    i = 1j
    u_max = 20.0
    du = u_max / N
    u = np.linspace(0, N-1, N) * du
    b = np.min(x)
    dx = 2.0 * np.pi / (N * du)
    x_i = b + np.linspace(0, N-1, N) * dx
    phi = np.exp(-i * b * u) * cf(u)
    gamma_1 = np.exp(-i * x_i * u[0]) * cf(u[0])
    gamma_2 = np.exp(-i * x_i * u[-1]) * cf(u[-1])
    phi_boundary = 0.5 * (gamma_1 + gamma_2)
    f_xi = du / np.pi * np.real(np.fft.fft(phi) - phi_boundary)
    return np.interp(x, x_i, f_xi)

def plot_density_comparison(mu, sigma, a, b, method, distribution, N_values):
    if distribution == 'Normal':
        x = np.linspace(a, b, 1000)
        f_XExact = stats.norm.pdf(x, mu, sigma)
        cF = lambda u: np.exp(1j * mu * u - 0.5 * sigma**2 * u**2)
    else:  # LogNormal
        x = np.linspace(0.05, 5, 1000)
        f_XExact = stats.lognorm.pdf(x, s=sigma, scale=np.exp(mu))
        cF = lambda u: np.exp(1j * mu * u - 0.5 * sigma**2 * u**2)
    
    fig = go.Figure()
    
    colors = ['blue', 'orange', 'green', 'red', 'purple']
    
    for N, color in zip(N_values, colors):
        if method == 'COS':
            if distribution == 'Normal':
                f_X = COSDensity(cF, x, N, a, b)
            else:
                f_X = 1/x * COSDensity(cF, np.log(x), N, a, b)
        else:
            f_X = FFTDensity(cF, x, N)
        fig.add_trace(go.Scatter(x=x, y=f_X, mode='lines', name=f'N = {N}', line=dict(color=color)))
    
    fig.add_trace(go.Scatter(x=x, y=f_XExact, mode='lines', name='Exact', line=dict(color='black', dash='dash')))
    
    fig.update_layout(
        title=f"{distribution} Density Recovery using {method} Method<br>(μ={mu:.2f}, σ={sigma:.2f}, a={a:.1f}, b={b:.1f})",
        xaxis_title="x" if distribution == 'Normal' else "y",
        yaxis_title="f_X(x)" if distribution == 'Normal' else "f_Y(y)",
        legend_title="Legend",
        font=dict(size=14),
        xaxis=dict(range=[a, b] if distribution == 'Normal' else [0, 5]),
        yaxis=dict(range=[0, max(f_XExact)*1.1])
    )
    
    return fig

st.title("Density Recovery using COS and FFT Methods")

st.sidebar.header("Parameters")
distribution = st.sidebar.radio("Distribution", ['Normal', 'LogNormal'])
mu = st.sidebar.slider("Mean (μ)", -5.0, 5.0, 0.5, 0.1)
sigma = st.sidebar.slider("Standard Deviation (σ)", 0.1, 5.0, 0.2, 0.1)
a = st.sidebar.slider("Lower bound (a)", -20.0, 0.0, -10.0, 0.5)
b = st.sidebar.slider("Upper bound (b)", 0.0, 20.0, 10.0, 0.5)
method = st.sidebar.radio("Method", ['COS', 'FFT'])

if method == 'COS':
    N_values = st.sidebar.multiselect("N values", [4, 8, 16, 32, 64, 128], default=[16, 64, 128])
else:
    N_values = st.sidebar.multiselect("N values", [64, 128, 256, 512, 1024], default=[64, 256, 1024])

st.plotly_chart(plot_density_comparison(mu, sigma, a, b, method, distribution, N_values))

if distribution == 'Normal':
    st.markdown("""
    ## Normal Distribution
    
    The probability density function (PDF) of a Normal distribution is:
    """)
    st.latex(r"f_X(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(x - \mu)^2}{2\sigma^2}\right)")
else:
    st.markdown("""
    ## LogNormal Distribution
    
    A random variable Y follows a LogNormal distribution if X = ln(Y) follows a Normal distribution. 
    If X ~ N(μ, σ²), then Y ~ LogNormal(μ, σ²).
    
    The probability density function (PDF) of a LogNormal distribution is:
    """)
    st.latex(r"f_Y(y) = \frac{1}{y\sigma\sqrt{2\pi}} \exp\left(-\frac{(\ln y - \mu)^2}{2\sigma^2}\right)")

st.markdown("""
### Characteristic Function

The characteristic function of the distribution is:
""")
st.latex(r"\phi_X(u) = E[e^{iuX}] = \exp\left(iu\mu - \frac{1}{2}\sigma^2u^2\right)")

if method == 'COS':
    st.markdown("""
    ## COS Method Theory

    The COS method approximates the density function as:
    """)
    st.latex(r"f_X(x) \approx \sum_{k=0}^{N-1}{}' A_k \cos(k\pi \frac{x-a}{b-a})")
    st.markdown("""
    where:
    - $[a,b]$ is the truncation range
    - $N$ is the number of expansion terms
    - $A_k$ are the Fourier-cosine series coefficients

    The coefficients $A_k$ are approximated using the characteristic function:
    """)
    st.latex(r"A_k \approx \frac{2}{b-a} \text{Re}\left\{\phi_X\left(\frac{k\pi}{b-a}\right) e^{-ika\pi/(b-a)}\right\}")
else:
    st.markdown("""
    ## FFT Method Theory

    The FFT method recovers the density function using the inverse Fourier transform:
    """)
    st.latex(r"f_X(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} e^{-iux} \phi_X(u) du")
    st.markdown("""
    where $\phi_X(u)$ is the characteristic function.

    The FFT algorithm is used to efficiently compute the discrete version of this integral:
    """)
    st.latex(r"f_X(x_j) \approx \frac{\Delta u}{2\pi} \sum_{k=0}^{N-1} e^{-iu_k x_j} \phi_X(u_k)")

st.markdown("""
## Plot Insights

1. **Convergence**: As N increases, both methods converge to the exact PDF (black dashed line). However, the rate of convergence may differ between COS and FFT.

2. **Oscillations**: For lower N values, you may observe oscillations in the recovered PDF, especially at the tails of the distribution. These oscillations are more pronounced in the COS method for very low N.

3. **Accuracy at Tails**: The accuracy of the recovered PDF at the tails of the distribution is often a good indicator of the method's performance. Higher N values generally provide better accuracy at the tails.

4. **COS Method Specifics**:
   - The truncation range [a, b] significantly affects the accuracy. If the range is too narrow, it may cut off important parts of the distribution.
   - Even with low N, the COS method often captures the general shape of the PDF well.

5. **FFT Method Specifics**:
   - The FFT method typically requires higher N values to achieve the same accuracy as the COS method.
   - It may show artifacts near the boundaries of the plotted range, especially for lower N values.

6. **Computational Efficiency**: While not visible in the plot, the COS method is generally more computationally efficient for smooth functions like the normal and lognormal distributions.

7. **LogNormal Specifics**:
   - The LogNormal distribution is skewed, which can make accurate recovery more challenging, especially at the tails.
   - The COS method for LogNormal uses a change of variables, which can affect the interpretation of the truncation range [a, b].

Experiment with different parameters to see how they affect the recovery accuracy and observe these insights in action!
""")import streamlit as st
import numpy as np
import plotly.graph_objects as go
import scipy.stats as stats

def COSDensity(cf, x, N, a, b):
    i = 1j
    k = np.arange(N)
    u = k * np.pi / (b - a)
    f_k = 2.0 / (b - a) * np.real(cf(u) * np.exp(-i * u * a))
    f_k[0] *= 0.5
    f_x = np.dot(f_k, np.cos(np.outer(u, x - a)))
    return f_x

def FFTDensity(cf, x, N):
    i = 1j
    u_max = 20.0
    du = u_max / N
    u = np.linspace(0, N-1, N) * du
    b = np.min(x)
    dx = 2.0 * np.pi / (N * du)
    x_i = b + np.linspace(0, N-1, N) * dx
    phi = np.exp(-i * b * u) * cf(u)
    gamma_1 = np.exp(-i * x_i * u[0]) * cf(u[0])
    gamma_2 = np.exp(-i * x_i * u[-1]) * cf(u[-1])
    phi_boundary = 0.5 * (gamma_1 + gamma_2)
    f_xi = du / np.pi * np.real(np.fft.fft(phi) - phi_boundary)
    return np.interp(x, x_i, f_xi)

def plot_density_comparison(mu, sigma, a, b, method, distribution, N_values):
    if distribution == 'Normal':
        x = np.linspace(a, b, 1000)
        f_XExact = stats.norm.pdf(x, mu, sigma)
        cF = lambda u: np.exp(1j * mu * u - 0.5 * sigma**2 * u**2)
    else:  # LogNormal
        x = np.linspace(0.05, 5, 1000)
        f_XExact = stats.lognorm.pdf(x, s=sigma, scale=np.exp(mu))
        cF = lambda u: np.exp(1j * mu * u - 0.5 * sigma**2 * u**2)
    
    fig = go.Figure()
    
    colors = ['blue', 'orange', 'green', 'red', 'purple']
    
    for N, color in zip(N_values, colors):
        if method == 'COS':
            if distribution == 'Normal':
                f_X = COSDensity(cF, x, N, a, b)
            else:
                f_X = 1/x * COSDensity(cF, np.log(x), N, a, b)
        else:
            f_X = FFTDensity(cF, x, N)
        fig.add_trace(go.Scatter(x=x, y=f_X, mode='lines', name=f'N = {N}', line=dict(color=color)))
    
    fig.add_trace(go.Scatter(x=x, y=f_XExact, mode='lines', name='Exact', line=dict(color='black', dash='dash')))
    
    fig.update_layout(
        title=f"{distribution} Density Recovery using {method} Method<br>(μ={mu:.2f}, σ={sigma:.2f}, a={a:.1f}, b={b:.1f})",
        xaxis_title="x" if distribution == 'Normal' else "y",
        yaxis_title="f_X(x)" if distribution == 'Normal' else "f_Y(y)",
        legend_title="Legend",
        font=dict(size=14),
        xaxis=dict(range=[a, b] if distribution == 'Normal' else [0, 5]),
        yaxis=dict(range=[0, max(f_XExact)*1.1])
    )
    
    return fig

st.title("Density Recovery using COS and FFT Methods")

st.sidebar.header("Parameters")
distribution = st.sidebar.radio("Distribution", ['Normal', 'LogNormal'])
mu = st.sidebar.slider("Mean (μ)", -5.0, 5.0, 0.5, 0.1)
sigma = st.sidebar.slider("Standard Deviation (σ)", 0.1, 5.0, 0.2, 0.1)
a = st.sidebar.slider("Lower bound (a)", -20.0, 0.0, -10.0, 0.5)
b = st.sidebar.slider("Upper bound (b)", 0.0, 20.0, 10.0, 0.5)
method = st.sidebar.radio("Method", ['COS', 'FFT'])

if method == 'COS':
    N_values = st.sidebar.multiselect("N values", [4, 8, 16, 32, 64, 128], default=[16, 64, 128])
else:
    N_values = st.sidebar.multiselect("N values", [64, 128, 256, 512, 1024], default=[64, 256, 1024])

st.plotly_chart(plot_density_comparison(mu, sigma, a, b, method, distribution, N_values))

if distribution == 'Normal':
    st.markdown("""
    ## Normal Distribution
    
    The probability density function (PDF) of a Normal distribution is:
    """)
    st.latex(r"f_X(x) = \frac{1}{\sigma\sqrt{2\pi}} \exp\left(-\frac{(x - \mu)^2}{2\sigma^2}\right)")
else:
    st.markdown("""
    ## LogNormal Distribution
    
    A random variable Y follows a LogNormal distribution if X = ln(Y) follows a Normal distribution. 
    If X ~ N(μ, σ²), then Y ~ LogNormal(μ, σ²).
    
    The probability density function (PDF) of a LogNormal distribution is:
    """)
    st.latex(r"f_Y(y) = \frac{1}{y\sigma\sqrt{2\pi}} \exp\left(-\frac{(\ln y - \mu)^2}{2\sigma^2}\right)")

st.markdown("""
### Characteristic Function

The characteristic function of the distribution is:
""")
st.latex(r"\phi_X(u) = E[e^{iuX}] = \exp\left(iu\mu - \frac{1}{2}\sigma^2u^2\right)")

if method == 'COS':
    st.markdown("""
    ## COS Method Theory

    The COS method approximates the density function as:
    """)
    st.latex(r"f_X(x) \approx \sum_{k=0}^{N-1}{}' A_k \cos(k\pi \frac{x-a}{b-a})")
    st.markdown("""
    where:
    - $[a,b]$ is the truncation range
    - $N$ is the number of expansion terms
    - $A_k$ are the Fourier-cosine series coefficients

    The coefficients $A_k$ are approximated using the characteristic function:
    """)
    st.latex(r"A_k \approx \frac{2}{b-a} \text{Re}\left\{\phi_X\left(\frac{k\pi}{b-a}\right) e^{-ika\pi/(b-a)}\right\}")
else:
    st.markdown("""
    ## FFT Method Theory

    The FFT method recovers the density function using the inverse Fourier transform:
    """)
    st.latex(r"f_X(x) = \frac{1}{2\pi} \int_{-\infty}^{\infty} e^{-iux} \phi_X(u) du")
    st.markdown("""
    where $\phi_X(u)$ is the characteristic function.

    The FFT algorithm is used to efficiently compute the discrete version of this integral:
    """)
    st.latex(r"f_X(x_j) \approx \frac{\Delta u}{2\pi} \sum_{k=0}^{N-1} e^{-iu_k x_j} \phi_X(u_k)")

st.markdown("""
## Plot Insights

1. **Convergence**: As N increases, both methods converge to the exact PDF (black dashed line). However, the rate of convergence may differ between COS and FFT.

2. **Oscillations**: For lower N values, you may observe oscillations in the recovered PDF, especially at the tails of the distribution. These oscillations are more pronounced in the COS method for very low N.

3. **Accuracy at Tails**: The accuracy of the recovered PDF at the tails of the distribution is often a good indicator of the method's performance. Higher N values generally provide better accuracy at the tails.

4. **COS Method Specifics**:
   - The truncation range [a, b] significantly affects the accuracy. If the range is too narrow, it may cut off important parts of the distribution.
   - Even with low N, the COS method often captures the general shape of the PDF well.

5. **FFT Method Specifics**:
   - The FFT method typically requires higher N values to achieve the same accuracy as the COS method.
   - It may show artifacts near the boundaries of the plotted range, especially for lower N values.

6. **Computational Efficiency**: While not visible in the plot, the COS method is generally more computationally efficient for smooth functions like the normal and lognormal distributions.

7. **LogNormal Specifics**:
   - The LogNormal distribution is skewed, which can make accurate recovery more challenging, especially at the tails.
   - The COS method for LogNormal uses a change of variables, which can affect the interpretation of the truncation range [a, b].

Experiment with different parameters to see how they affect the recovery accuracy and observe these insights in action!
""")