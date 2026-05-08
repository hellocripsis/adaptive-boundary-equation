#!/usr/bin/env python3
"""Generate Figure 5: critical slowing near the ABE viability threshold."""

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np


def main() -> None:
    # Fixed normalized ABE parameters for the paper-validation sweep.
    gamma = 1.0
    lambda_ = 1.0
    xi = 1.0
    beta = 1.0

    # Sweep R just above threshold. With gamma=lambda=xi=1, eta=R.
    delta = np.logspace(-4, np.log10(2.0), 300)
    r_values = 1.0 + delta

    tau_slow_values = []

    for r_value in r_values:
        eta = r_value
        eta_c = xi * lambda_ / gamma

        if eta <= eta_c:
            raise ValueError("Active branch only exists for R > 1.")

        i_star = np.sqrt((gamma * eta / xi - lambda_) / beta)
        s_star = (eta / xi) * i_star

        # The active-branch Jacobian has determinant
        # D_star = 2*(gamma*eta - lambda*xi) = 2*xi*lambda*(R - 1).
        # Therefore D_star tends to zero as R -> 1+, driving the dominant
        # relaxation rate toward zero and making tau_slow large near threshold.
        j_star = np.array(
            [
                [2.0 * lambda_ - 3.0 * gamma * eta / xi, gamma],
                [eta, -xi],
            ]
        )

        eigenvalues = np.linalg.eigvals(j_star)
        real_parts = np.real(eigenvalues)
        slow_rate = np.min(np.abs(real_parts))
        tau_slow = 1.0 / slow_rate

        # Keep the active fixed point variables explicit for validation clarity.
        _ = (i_star, s_star)
        tau_slow_values.append(tau_slow)

    tau_slow_values = np.array(tau_slow_values)

    # Reference trend proportional to 1/(R - 1), normalized to the first point.
    reference_trend = tau_slow_values[0] * delta[0] / delta

    fig, ax = plt.subplots(figsize=(7.0, 4.8))
    ax.loglog(delta, tau_slow_values, label="Jacobian eigenvalue relaxation time")
    ax.loglog(
        delta,
        reference_trend,
        linestyle="--",
        label="reference trend proportional to 1/(R - 1)",
    )
    ax.set_xlabel("R - 1")
    ax.set_ylabel("dominant relaxation time tau_slow")
    ax.set_title("Critical slowing near ABE threshold")
    ax.grid(True, which="both", linestyle=":", linewidth=0.7)
    ax.legend()
    fig.tight_layout()

    output_path = (
        Path(__file__).resolve().parents[1]
        / "outputs"
        / "figure5_critical_slowing.png"
    )
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(output_path, dpi=200)
    plt.close(fig)

    print(f"sweep_points: {len(r_values)}")
    print(f"min_R: {np.min(r_values):.6f}")
    print(f"max_R: {np.max(r_values):.6f}")
    print(f"max_tau_slow: {np.max(tau_slow_values):.6f}")
    print(f"min_tau_slow: {np.min(tau_slow_values):.6f}")
    print(f"output_path: {output_path}")


if __name__ == "__main__":
    main()
