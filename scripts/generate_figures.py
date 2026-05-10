import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from numpy.linalg import eigvals
from pathlib import Path

plt.rcParams.update({
    'font.size': 16,
    'axes.titlesize': 17,
    'axes.labelsize': 16,
    'xtick.labelsize': 14,
    'ytick.labelsize': 14,
    'legend.fontsize': 13,
    'figure.titlesize': 17,
})

OUT = Path(__file__).resolve().parents[1] / 'figures'
OUT.mkdir(parents=True, exist_ok=True)

# ABE dynamics: dI/dt = gamma*S - lambda*I - beta*I^3, dS/dt = eta*I - xi*S

def f(t, y, gamma=1, lam=1, beta=1, eta=0.7, xi=1):
    I, S = y
    return [gamma*S - lam*I - beta*I**3, eta*I - xi*S]

# Figure 1: phase portraits below and above threshold
fig, axes = plt.subplots(1, 2, figsize=(13.0, 6.2))
settings = [(0.7, 'R < 1: quiescent branch stable'), (2.5, 'R > 1: active branch stable')]
for ax, (eta, title) in zip(axes, settings):
    gamma=lam=beta=xi=1.0
    x = np.linspace(-2, 2, 25)
    y = np.linspace(-2, 2, 25)
    X, Y = np.meshgrid(x, y)
    U = gamma*Y - lam*X - beta*X**3
    V = eta*X - xi*Y
    speed = np.sqrt(U**2 + V**2)
    U2 = U/(speed+1e-8)
    V2 = V/(speed+1e-8)
    ax.quiver(X, Y, U2, V2, angles='xy', alpha=0.35, width=0.0025)
    initials = [(-1.8, 1.8), (1.8, 1.8), (-1.8,-1.8), (1.8,-1.8), (-0.7,0.2),(0.7,-0.2)]
    for y0 in initials:
        sol = solve_ivp(lambda t,z: f(t,z,gamma,lam,beta,eta,xi), (0, 20), y0, max_step=0.05, rtol=1e-9, atol=1e-11)
        ax.plot(sol.y[0], sol.y[1], lw=1.0, alpha=0.85)
    ax.axhline(0, lw=0.5, alpha=0.5)
    ax.axvline(0, lw=0.5, alpha=0.5)
    ax.scatter([0],[0], marker='o', s=65, label='quiescent equilibrium' if eta<1 else 'unstable saddle', zorder=4)
    if eta > 1:
        Istar = np.sqrt((gamma*eta/xi - lam)/beta)
        Sstar = (eta/xi)*Istar
        ax.scatter([Istar, -Istar],[Sstar, -Sstar], marker='o', s=65, label='active equilibria', zorder=4)
    ax.set_title(title)
    ax.set_xlabel('I')
    ax.set_ylabel('S')
    ax.set_xlim(-2,2); ax.set_ylim(-2,2)
    ax.legend(fontsize=12, loc='upper left')
fig.tight_layout()
fig.savefig(OUT/'figure1_phase_portrait.pdf', bbox_inches='tight')
fig.savefig(OUT/'figure1_phase_portrait.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# Figure 2: amplitude scaling from ODE integration
etas = np.linspace(1.01, 4.1, 90)
gamma=lam=beta=xi=1.0
xvals=[]; yvals=[]
rng = np.random.default_rng(4)
for eta in etas:
    y0 = rng.normal(0.05,0.01,2)
    y0[0] = abs(y0[0]) + 0.01
    sol = solve_ivp(lambda t,z: f(t,z,gamma,lam,beta,eta,xi), (0, 80), y0, rtol=1e-9, atol=1e-11, max_step=0.1)
    Iend = sol.y[0,-1]
    xvals.append(eta - 1.0)
    yvals.append(Iend**2)
fig, ax = plt.subplots(figsize=(10.5,6.6))
ax.scatter(xvals, yvals, s=20, alpha=0.7, label='ODE integration endpoint')
xx=np.linspace(0, max(xvals), 200)
ax.plot(xx, xx, label=r'Theory: $(I^*)^2=\eta-\eta_c$')
ax.set_title('Active-branch amplitude from numerical integration')
ax.set_xlabel(r'$\eta-\eta_c$')
ax.set_ylabel('Measured $(I_{end})^2$')
ax.legend(fontsize=13)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(OUT/'figure2_amplitude_scaling.pdf', bbox_inches='tight')
fig.savefig(OUT/'figure2_amplitude_scaling.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# Figure 3: random stability check
rng=np.random.default_rng(12)
N=1000
Rs=[]; eigmax=[]; colors=[]
for _ in range(N):
    gamma=10**rng.uniform(-1,1)
    lam=10**rng.uniform(-1,1)
    beta=10**rng.uniform(-1,1)
    xi=10**rng.uniform(-1,1)
    R=10**rng.uniform(-4,5)
    eta=R*xi*lam/gamma
    if R <= 1:
        J=np.array([[-lam, gamma],[eta, -xi]])
        colors.append(0)
    else:
        J=np.array([[2*lam - 3*gamma*eta/xi, gamma],[eta, -xi]])
        colors.append(1)
    Rs.append(R)
    eigmax.append(np.max(eigvals(J).real))
fig, ax = plt.subplots(figsize=(10.5,6.6))
Rs=np.array(Rs); eigmax=np.array(eigmax); colors=np.array(colors)
ax.scatter(Rs[colors==0], eigmax[colors==0], s=14, alpha=0.45, label='quiescent branch used')
ax.scatter(Rs[colors==1], eigmax[colors==1], s=14, alpha=0.45, label='active branch used')
ax.axhline(0, ls='--', lw=1)
ax.axvline(1, ls='--', lw=1)
ax.set_xscale('log')
ax.set_title('Stability check on synthetic parameter sets')
ax.set_xlabel(r'$R=\gamma\eta/(\xi\lambda)$')
ax.set_ylabel('Largest real eigenvalue')
ax.legend(fontsize=13)
ax.grid(alpha=0.25)
fig.tight_layout()
fig.savefig(OUT/'figure3_stability_check.pdf', bbox_inches='tight')
fig.savefig(OUT/'figure3_stability_check.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# Figure 4: synthetic operational demonstration
rng=np.random.default_rng(22)
eta_vals=np.linspace(0.4,3.2,60)
noise_sigma=0.04
Rvals=[]; measured_abs=[]; above_x=[]; above_y=[]
for eta in eta_vals:
    R=eta
    if eta <= 1:
        Itrue=0.0
    else:
        Itrue=np.sqrt(eta-1.0)
    Imeas=abs(Itrue + rng.normal(0, noise_sigma))
    Rvals.append(R); measured_abs.append(Imeas)
    if eta > 1:
        above_x.append(eta-1.0)
        above_y.append(Imeas**2)
fig, axes = plt.subplots(1,2, figsize=(13.0,6.2))
axes[0].scatter(Rvals, measured_abs, s=24, alpha=0.75, label='measured equilibrium $|I|$')
rr=np.linspace(min(Rvals), max(Rvals), 300)
theory=np.where(rr>1, np.sqrt(rr-1), 0)
axes[0].plot(rr, theory, label=r'theory: $|I|=\sqrt{R-1}$')
axes[0].axvline(1, ls='--', lw=1)
axes[0].set_title('Threshold crossing in synthetic operational data')
axes[0].set_xlabel(r'$R=\gamma\eta/(\xi\lambda)$')
axes[0].set_ylabel('Measured $|I|$ at equilibrium')
axes[0].legend(fontsize=11)
axes[0].grid(alpha=0.25)
axes[1].scatter(above_x, above_y, s=24, alpha=0.75, label='measured $|I|^2$ above threshold')
xx=np.linspace(0, max(above_x), 200)
axes[1].plot(xx, xx, label=r'theory: $|I|^2=\eta-\eta_c$')
axes[1].set_title('Amplitude scaling above threshold')
axes[1].set_xlabel(r'$\eta-\eta_c$')
axes[1].set_ylabel('Measured $|I|^2$ at equilibrium')
axes[1].legend(fontsize=11)
axes[1].grid(alpha=0.25)
fig.tight_layout()
fig.savefig(OUT/'figure4_synthetic_demo.pdf', bbox_inches='tight')
fig.savefig(OUT/'figure4_synthetic_demo.png', dpi=300, bbox_inches='tight')
plt.close(fig)

# Figure 5: critical slowing near threshold
# Fixed normalized parameters: gamma = xi = lambda = beta = 1.
# The active-branch determinant satisfies D* = 2*xi*lambda*(R - 1), so the
# dominant relaxation rate tends to zero as R -> 1+.
gamma = 1.0
lam = 1.0
xi = 1.0
beta = 1.0
delta = np.logspace(-4, np.log10(2.0), 300)
R_values = 1.0 + delta
tau_slow_values = []
for R in R_values:
    eta = R
    Istar = np.sqrt((gamma * eta / xi - lam) / beta)
    Sstar = (eta / xi) * Istar
    J = np.array([[2.0 * lam - 3.0 * gamma * eta / xi, gamma], [eta, -xi]])
    ev = eigvals(J)
    slow_rate = np.min(np.abs(ev.real))
    tau_slow_values.append(1.0 / slow_rate)
    _ = (Istar, Sstar)
tau_slow_values = np.array(tau_slow_values)
reference_trend = tau_slow_values[0] * delta[0] / delta
fig, ax = plt.subplots(figsize=(10.5, 6.6))
ax.loglog(delta, tau_slow_values, label='Active-branch Jacobian relaxation time')
ax.loglog(delta, reference_trend, linestyle='--', label=r'Reference trend $\propto 1/(R-1)$')
ax.set_xlabel(r'$R - 1$')
ax.set_ylabel(r'Dominant relaxation time $\tau_{\mathrm{slow}}$')
ax.set_title('Critical slowing near viability threshold')
ax.grid(True, which='both', linestyle=':', linewidth=0.7)
ax.legend(fontsize=13)
fig.tight_layout()
fig.savefig(OUT/'figure5_critical_slowing.pdf', bbox_inches='tight')
fig.savefig(OUT/'figure5_critical_slowing.png', dpi=300, bbox_inches='tight')
plt.close(fig)
