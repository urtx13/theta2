#!/usr/bin/env python3
"""
theta_kerr_param_map.py

Mapa 2D senzill en (a*, chi) per al framework Θ–Kerr:

- Es fixa (M, Λ, σ_K, σ_Π).
- Es defineix una llista de spins a* i una llista de ε (⇒ χ).
- Per a cada combinació:
    - es calcula <G>_H a partir de la geometria Kerr;
    - s'estima Δr_ph/r_ph;
    - s'avalua el residual de la First Law;
    - es fa triage (PASS / REJECT).
"""

from typing import List, Tuple

from theta_kerr_horizon_gate import KerrBH, horizon_gate_profile
from theta_kerr_diagnostics import (
    ThetaKerrParams,
    photon_sphere_shift_estimate,
    first_law_residual,
    triage_decision,
)


def param_map(
    M: float = 10.0,
    Lambda: float = 100.0,
    sigma_K: float = 0.1,
    sigma_Pi: float = 0.3,
    spins: List[float] | None = None,
    eps_values: List[float] | None = None,
) -> List[Tuple[float, float, float, float, float, str]]:
    """
    Retorna una llista de tuples:

      (a_star, chi, G_avg, dmin, dmax, decision)
    """

    if spins is None:
        spins = [0.5, 0.7, 0.9, 0.99]

    if eps_values is None:
        eps_values = [1.0e-4, 1.0e-3, 1.0e-2]  # chi ~ 1e-6, 1e-5, 1e-4

    results: List[Tuple[float, float, float, float, float, str]] = []

    # Bucle sobre spins
    for a_star in spins:
        # Per a cada spin, calculem <G>_H una sola vegada:
        params_ref = ThetaKerrParams(
            M=M,
            eps=eps_values[0],
            Lambda=Lambda,
            a_star=a_star,
            sigma_K=sigma_K,
            sigma_Pi=sigma_Pi,
        )
        bh_ref = KerrBH(M=params_ref.M, a=params_ref.a)

        _, _, _, G_avg = horizon_gate_profile(
            bh_ref,
            sigma_K=params_ref.sigma_K,
            sigma_Pi=params_ref.sigma_Pi,
            n_theta=400,
        )

        # Ara variem ε (⇒ chi), però <G>_H es manté
        for eps in eps_values:
            params = ThetaKerrParams(
                M=M,
                eps=eps,
                Lambda=Lambda,
                a_star=a_star,
                sigma_K=sigma_K,
                sigma_Pi=sigma_Pi,
            )

            chi = params.chi

            dmin, dmax = photon_sphere_shift_estimate(
                chi, G_avg, A_min=0.5, A_max=2.0
            )

            # Residual de First Law (placeholder coherent amb la resta)
            dM = 1.0
            TH_dS = 0.99
            OmegaH_dJ = 0.0
            FL_R = first_law_residual(
                dM, TH_dS, OmegaH_dJ, epsM=params.eps * params.M
            )

            decision, _ = triage_decision(G_avg, chi, FL_R)

            results.append((a_star, chi, G_avg, dmin, dmax, decision))

    return results


def demo() -> None:
    """
    Demo amb pseudo-progress bar i taula resum.
    """

    print("===========================================")
    print("  Θ–Kerr parameter map demo (a*, chi)      ")
    print("===========================================\n")

    steps = [
        "Configurar paràmetres globals (M, Λ, σ_K, σ_Π)",
        "Definir llistes de spins a* i ε",
        "Construir mapa en (a*, chi)",
        "Imprimir taula resum",
    ]
    total = len(steps)

    for i, step in enumerate(steps, start=1):
        bar = "[" + "#" * i + "-" * (total - i) + "]"
        print(f"{bar} {step}...")

        if i == 1:
            M = 10.0
            Lambda = 100.0
            sigma_K = 0.1
            sigma_Pi = 0.3

        elif i == 2:
            spins = [0.5, 0.7, 0.9, 0.99]
            eps_values = [1.0e-4, 1.0e-3, 1.0e-2]

        elif i == 3:
            results = param_map(
                M=M,
                Lambda=Lambda,
                sigma_K=sigma_K,
                sigma_Pi=sigma_Pi,
                spins=spins,
                eps_values=eps_values,
            )

        elif i == 4:
            print("\n=========== PARAMETER MAP (a*, chi) ===========")
            print(
                " a*   |   chi      |  <G>_H   |  Δr_ph/r_ph(min)  |  Δr_ph/r_ph(max)  | triage"
            )
            print(
                "------+-----------+----------+--------------------+--------------------+--------"
            )
            for (a_star, chi, G_avg, dmin, dmax, decision) in results:
                print(
                    f" {a_star:4.2f} | "
                    f"{chi:9.3e} | "
                    f"{G_avg:8.5f} | "
                    f"{dmin:18.3e} | "
                    f"{dmax:18.3e} | "
                    f"{decision}"
                )
            print("===============================================\n")

    print("Parameter map Θ–Kerr completada.\n")


if __name__ == "__main__":
    demo()
