#!/usr/bin/env python3
"""
theta_kerr_diagnostics.py

Mòdul de diagnòstic per al programa Θ–Kerr:

- Implementa el paràmetre EFT χ = ε M^2 / Λ^2.
- Crida theta_kerr_horizon_gate.horizon_gate_profile per obtenir <G>_H.
- Estima el desplaçament de l'esfera de fotons Δr_ph / r_ph ~ A χ <G>_H.
- Calcula el residual de la First Law.
- Aplica el protocol de triage (REJECT / PASS).
"""

from dataclasses import dataclass
import math

# Importem el codi del gate geomètric
from theta_kerr_horizon_gate import KerrBH, horizon_gate_profile


# ----------------------------
# 1. Paràmetres del model
# ----------------------------

@dataclass
class ThetaKerrParams:
    """
    Paràmetres bàsics de la teoria Θ–Kerr en unitats geomètriques (G = c = 1).

    M        : massa del forat negre
    eps      : paràmetre de couplatge ε
    Lambda   : escala EFT Λ
    a_star   : spin adimensional a* = a / M
    sigma_K  : amplada suau del gate en K
    sigma_Pi : amplada suau del gate en Π^2
    """
    M: float
    eps: float
    Lambda: float
    a_star: float = 0.5
    sigma_K: float = 0.1
    sigma_Pi: float = 0.3

    @property
    def chi(self) -> float:
        """
        χ = ε M^2 / Λ^2
        """
        return self.eps * (self.M ** 2) / (self.Lambda ** 2)

    @property
    def a(self) -> float:
        """
        Paràmetre de spin dimensional: a = a* M
        """
        return self.a_star * self.M


# ----------------------------
# 2. Gate geomètric (només wrappers numèrics)
# ----------------------------

def photon_sphere_shift_estimate(chi: float,
                                 G_avg: float,
                                 A_min: float = 0.5,
                                 A_max: float = 2.0) -> tuple[float, float]:
    """
    Estima Δr_ph / r_ph segons l'escalat del paper:

      Δr_ph / r_ph ≈ A(a*) * χ * <G>_H

    On A ∈ [A_min, A_max] ~ O(1).

    Retorna (Δr_min, Δr_max).
    """
    base = chi * G_avg
    return A_min * base, A_max * base


def first_law_residual(dM: float,
                       TH_dS: float,
                       OmegaH_dJ: float,
                       epsM: float) -> float:
    """
    R = |dM - TH dS - Ω_H dJ| / max(|dM|, |TH dS|, ε M)

    Residual normalitzat per fer el triage.
    """
    num = abs(dM - TH_dS - OmegaH_dJ)
    denom = max(abs(dM), abs(TH_dS), epsM)
    return num / denom if denom != 0.0 else 0.0


def triage_decision(G_avg: float,
                    chi: float,
                    FL_residual: float,
                    G_threshold: float = 0.10,
                    FL_threshold: float = 0.02) -> tuple[str, list[str]]:
    """
    Implementa el protocol de triage de la Taula 1:

      - Stage 1: gate coverage <G>_H
      - Stage 2: residual de la First Law R

    Retorna:
      (decisió, [llista de motius])
    """
    reasons: list[str] = []

    # Stage 1: gate coverage
    if G_avg < G_threshold:
        reasons.append(
            f"Gate coverage massa petit: <G>_H = {G_avg:.3e} < {G_threshold:.2f}"
        )

    # Stage 2: First Law residual
    if FL_residual > FL_threshold:
        reasons.append(
            f"Residual de la First Law massa gran: R = {FL_residual:.3e} > {FL_threshold:.2f}"
        )

    if reasons:
        return "REJECT", reasons
    else:
        return "PASS", ["Configuració apta per passar al solver de mètrica."]


# ----------------------------
# 3. Demo amb “progress bar”
# ----------------------------

def demo() -> None:
    """
    Petit exemple que reprodueix lògica de la Sec. 3:

    - Tria (M, ε, Λ, a*).
    - Calcula χ.
    - Calcula <G>_H dinàmicament amb horizon_gate_profile.
    - Estima Δr_ph / r_ph.
    - Avalua residual de la First Law.
    - Aplica triage.
    """

    print("===========================================")
    print("  Θ–Kerr diagnostic demo (slow/fast spin) ")
    print("===========================================\n")

    steps = [
        "Fixar paràmetres (M, ε, Λ, a*)",
        "Calcular χ i <G>_H dinàmicament",
        "Estimar Δr_ph / r_ph",
        "Avaluar residual de la First Law",
        "Aplicar triage (REJECT / PASS)",
        "Resum final"
    ]

    total = len(steps)

    # Variables que anirem omplint
    params: ThetaKerrParams | None = None
    chi: float | None = None
    G_avg: float | None = None
    dmin = dmax = None
    FL_R: float | None = None
    decision: str | None = None
    reasons: list[str] = []

    for i, step in enumerate(steps, start=1):
        bar = "[" + "#" * i + "-" * (total - i) + "]"
        print(f"{bar} {step}...")

        if i == 1:
            # Exemple: M = 10, ε = 1e-3, Λ = 100, a* = 0.99
            params = ThetaKerrParams(
                M=10.0,
                eps=1.0e-3,
                Lambda=100.0,
                a_star=0.99,
                sigma_K=0.1,
                sigma_Pi=0.3,
            )

        elif i == 2:
            chi = params.chi

            # Construir el BH de Kerr i calcular <G>_H
            bh = KerrBH(M=params.M, a=params.a)
            _, _, _, G_avg = horizon_gate_profile(
                bh,
                sigma_K=params.sigma_K,
                sigma_Pi=params.sigma_Pi,
                n_theta=400
            )

        elif i == 3:
            dmin, dmax = photon_sphere_shift_estimate(
                chi, G_avg, A_min=0.5, A_max=2.0
            )

        elif i == 4:
            # Valors ficticis per il·lustrar el residual de la First Law.
            # Quan tinguis el codi termodinàmic detallat, aquí hi posaràs
            # dM, TH dS i Ω_H dJ reals.
            dM = 1.0
            TH_dS = 0.99
            OmegaH_dJ = 0.0
            FL_R = first_law_residual(
                dM, TH_dS, OmegaH_dJ,
                epsM=params.eps * params.M
            )

        elif i == 5:
            decision, reasons = triage_decision(G_avg, chi, FL_R)

        elif i == 6:
            print("\n=========== SUMMARY ===========")
            print(f"M                 = {params.M:.3f}")
            print(f"a*                = {params.a_star:.3f}")
            print(f"a                 = {params.a:.3f}")
            print(f"ε                 = {params.eps:.3e}")
            print(f"Λ                 = {params.Lambda:.3f}")
            print(f"χ = ε M^2 / Λ^2   = {chi:.3e}")
            print(f"<G>_H            = {G_avg:.5f}")
            print(f"Δr_ph/r_ph (min) = {dmin:.3e}")
            print(f"Δr_ph/r_ph (max) = {dmax:.3e}")
            print(f"First Law R      = {FL_R:.3e}")
            print(f"Decisió triage   = {decision}")
            for r in reasons:
                print(f"  - {r}")
            print("================================\n")

    print("Demo Θ–Kerr completada.\n")


# ----------------------------
# 4. Punt d'entrada
# ----------------------------

if __name__ == "__main__":
    demo()
