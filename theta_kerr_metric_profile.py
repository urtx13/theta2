#!/usr/bin/env python3
"""
theta_kerr_sigma_scan.py

Scan en la "sharpness" del gate, sigma_Pi, per a un forat negre
quasi extremal (a* ~ 0.99) i un valor fix d'χ.

Objectiu:
    - Comprovar com canvien <G>_H i Δr_ph/r_ph quan fem
      el gate més suau o més brusc en Π.
"""

from typing import List, Tuple

from theta_kerr_horizon_gate import KerrBH, horizon_gate_profile
from theta_kerr_diagnostics import (
    ThetaKerrParams,
    photon_sphere_shift_estimate,
    first_law_residual,
    triage_decision,
)


def sigma_scan(
    M: float = 10.0,
    a_star: float = 0.99,
    eps: float = 1.0e-3,   # χ ~ 1e-5 per Λ=100
    Lambda: float = 100.0,
    sigma_K: float = 0.1,
    sigma_Pi_list: List[float] | None = None,
) -> List[Tuple[float, float, float, float, float, str]]:
    """
    Retorna una llista de tuples:

      (sigma_Pi, chi, <G>_H, Δr_ph/r_ph_min, Δr_ph/r_ph_max, decisio)
    """

    if sigma_Pi_list is None:
        sigma_Pi_list = [0.05, 0.10, 0.20, 0.30, 0.50]

    results: List[Tuple[float, float, float, float, float, str]] = []

    for sigma_Pi in sigma_Pi_list:
        # Paràmetres EFT per aquest sigma_Pi
        params = ThetaKerrParams(
            M=M,
            eps=eps,
            Lambda=Lambda,
            a_star=a_star,
            sigma_K=sigma_K,
            sigma_Pi=sigma_Pi,
        )

        chi = params.chi

        # Forat negre de Kerr
        bh = KerrBH(M=params.M, a=params.a)

        # Càlcul de <G>_H amb aquest sigma_Pi
        _, _, _, G_avg = horizon_gate_profile(
            bh,
            sigma_K=params.sigma_K,
            sigma_Pi=params.sigma_Pi,
            n_theta=400,
        )

        # Estimació Δr_ph/r_ph
        dmin, dmax = photon_sphere_shift_estimate(
            chi, G_avg, A_min=0.5, A_max=2.0
        )

        # Residual de la First Law (mateix placeholder que a la resta)
        dM = 1.0
        TH_dS = 0.99
        OmegaH_dJ = 0.0
        FL_R = first_law_residual(
            dM, TH_dS, OmegaH_dJ, epsM=params.eps * params.M
        )

        # Triatge
        decision, _ = triage_decision(G_avg, chi, FL_R)

        results.append((sigma_Pi, chi, G_avg, dmin, dmax, decision))

    return results


def demo() -> None:
    """
    Demo amb "progress bar" i taula resum.
    """

    print("===========================================")
    print("  Θ–Kerr sigma_Pi scan demo (a* ~ 0.99)   ")
    print("===========================================\n")

    steps = [
        "Configurar paràmetres globals (M, a*, ε, Λ, σ_K)",
        "Definir valors de sigma_Pi",
        "Executar sigma-scan",
        "Imprimir taula resum",
    ]
    total = len(steps)

    for i, step in enumerate(steps, start=1):
        bar = "[" + "#" * i + "-" * (total - i) + "]"
        print(f"{bar} {step}...")

        if i == 1:
            M = 10.0
            a_star = 0.99
            eps = 1.0e-3
            Lambda = 100.0
            sigma_K = 0.1

        elif i == 2:
            sigma_Pi_list = [0.05, 0.10, 0.20, 0.30, 0.50]

        elif i == 3:
            results = sigma_scan(
                M=M,
                a_star=a_star,
                eps=eps,
                Lambda=Lambda,
                sigma_K=sigma_K,
                sigma_Pi_list=sigma_Pi_list,
            )

        elif i == 4:
            print("\n====== SIGMA_Pi SCAN (a* = 0.99, χ ~ 1e-5) ======")
            print(
                " sigma_Pi |   chi      |  <G>_H   |  Δr_ph/r_ph(min)  |  Δr_ph/r_ph(max)  | triage"
            )
            print(
                "----------+-----------+----------+--------------------+--------------------+--------"
            )
            for (sigma_Pi, chi, G_avg, dmin, dmax, decision) in results:
                print(
                    f"  {sigma_Pi:6.2f} | "
                    f"{chi:9.3e} | "
                    f"{G_avg:8.5f} | "
                    f"{dmin:18.3e} | "
                    f"{dmax:18.3e} | "
                    f"{decision}"
                )
            print("=================================================\n")

    print("Sigma_Pi scan Θ–Kerr completada.\n")


if __name__ == "__main__":
    demo()
