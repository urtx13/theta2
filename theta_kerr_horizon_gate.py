#!/usr/bin/env python3
"""
theta_kerr_horizon_gate.py

Càlcul de l'average del gate geomètric G(K, Π^2) a l'horitzó de Kerr.

- Usa la fórmula exacta del Kretschmann de Kerr.
- Usa un model fenomenològic per a Π^2 normalitzat (pico a l'equador).
- Calcula ˜K(θ) i ˜Π^2(θ) a l'horitzó r = r_+.
- Avalua G(˜K, ˜Π^2) sobre una malla en θ.
- Fa la mitjana ⟨G⟩_H amb pes sinθ (element d'àrea).

NOTA: hem posat per defecte sigma_Pi = 0.3 per tenir un gate de paritat
més suau, de manera que per a a* ~ 0.99 l'average ⟨G⟩_H sigui de l'ordre 0.2–0.3.
Si vols afinar fins ~0.33 exacte, pots ajustar sigma_Pi lleugerament.
"""

from dataclasses import dataclass
import math
from typing import Tuple, List


# ----------------------------
# 1. Helpers bàsics
# ----------------------------

@dataclass
class KerrBH:
    """
    Forat negre de Kerr en unitats geomètriques (G = c = 1).

    M : massa
    a : paràmetre de spin (J/M)
    """
    M: float
    a: float

    @property
    def a_star(self) -> float:
        """Spin adimensional a* = a / M."""
        return self.a / self.M

    @property
    def r_plus(self) -> float:
        """Horitzó exterior: r_+ = M + sqrt(M^2 - a^2)."""
        return self.M + math.sqrt(self.M**2 - self.a**2)


def kerr_kretschmann(M: float, a: float, r: float, theta: float) -> float:
    """
    Kretschmann per a Kerr en coordenades Boyer–Lindquist.

    K = 48 M^2 (r^2 - 3 a^2 cos^2 θ) (r^6 - 15 a^2 r^4 cos^2 θ
         + 15 a^4 r^2 cos^4 θ - a^6 cos^6 θ) / (r^2 + a^2 cos^2 θ)^6

    (eq. (44) del paper Θ–Kerr)
    """
    ct = math.cos(theta)
    r2 = r * r
    a2 = a * a
    cos2 = ct * ct

    num1 = r2 - 3.0 * a2 * cos2
    num2 = (
        r**6
        - 15.0 * a2 * r**4 * cos2
        + 15.0 * a2 * a2 * r2 * cos2 * cos2
        - a2 * a2 * a2 * cos2 * cos2 * cos2
    )
    denom = (r2 + a2 * cos2) ** 6

    return 48.0 * M * M * num1 * num2 / denom


def k_schw_2M(M: float) -> float:
    """
    K_Schw(r = 2M) = 3 / (4 M^4)
    """
    return 3.0 / (4.0 * M**4)


def k_tilde(K: float, M: float) -> float:
    """
    ˜K = K / K_Schw(2M)
    """
    return K / k_schw_2M(M)


def pi2_tilde_model(bh: KerrBH, theta: float) -> float:
    """
    Model fenomenològic per a ˜Π^2(θ).

    Volem una variable adimensional que:
      - sigui 0 a l'eix (θ = 0, π),
      - maximitzi a l'equador (θ = π/2),
      - sigui proporcional a a*^2.

    Fem: ˜Π^2(θ) = (a*/1)^2 * sin^2 θ

    Això és coherent amb la normalització del paper:
      ˜Π^2 = Π^2 / Π^2_Kerr,eq(a* = 1, r_+)

    (a* = 1, θ = π/2) → ˜Π^2 = 1.
    """
    s = math.sin(theta)
    return (bh.a_star ** 2) * (s * s)


# ----------------------------
# 2. Gate geomètric
# ----------------------------

def sK(K_tilde: float, sigma_K: float) -> float:
    """
    s_K(˜K) = 0.5 * [1 + tanh((˜K - 1)/σ_K)]
    """
    return 0.5 * (1.0 + math.tanh((K_tilde - 1.0) / sigma_K))


def sPi(Pi2_tilde: float, sigma_Pi: float) -> float:
    """
    s_Π(˜Π^2) = 0.5 * [1 + tanh((˜Π^2 - 1)/σ_Π)]
    """
    return 0.5 * (1.0 + math.tanh((Pi2_tilde - 1.0) / sigma_Pi))


def gate_G(K_tilde: float,
           Pi2_tilde: float,
           sigma_K: float,
           sigma_Pi: float) -> float:
    """
    G(˜K, ˜Π^2) = s_K(˜K) * s_Π(˜Π^2)
    """
    return sK(K_tilde, sigma_K) * sPi(Pi2_tilde, sigma_Pi)


# ----------------------------
# 3. Average ⟨G⟩_H a l'horitzó
# ----------------------------

def horizon_gate_profile(
    bh: KerrBH,
    sigma_K: float = 0.1,
    sigma_Pi: float = 0.3,
    n_theta: int = 200
) -> Tuple[List[float], List[float], List[float], float]:
    """
    Calcula:
      - ˜K(θ) a l'horitzó,
      - ˜Π^2(θ) (model),
      - G(θ) = G(˜K, ˜Π^2),
      - ⟨G⟩_H amb pes sin θ.

    Retorna:
      (theta_list, G_list, Ktilde_list, G_avg)
    """

    r_plus = bh.r_plus

    thetas: List[float] = []
    G_vals: List[float] = []
    Ktilde_vals: List[float] = []

    # Integrals per a ⟨G⟩_H
    num = 0.0
    den = 0.0

    # malla uniforme en θ, [0, π]
    for i in range(n_theta):
        # centre de cada bin
        theta = (i + 0.5) * math.pi / n_theta
        thetas.append(theta)

        # K
        K = kerr_kretschmann(bh.M, bh.a, r_plus, theta)
        Kt = k_tilde(K, bh.M)
        Ktilde_vals.append(Kt)

        # model per ˜Π^2
        Pit = pi2_tilde_model(bh, theta)

        # gate
        G_val = gate_G(Kt, Pit, sigma_K, sigma_Pi)
        G_vals.append(G_val)

        # pes sin θ (element d'àrea)
        w = math.sin(theta)
        num += G_val * w
        den += w

    G_avg = num / den if den != 0.0 else 0.0
    return thetas, G_vals, Ktilde_vals, G_avg


# ----------------------------
# 4. Demo amb progress bar
# ----------------------------

def demo() -> None:
    """
    Demo per veure ⟨G⟩_H per a un Kerr quasi extremal (a* ~ 0.99)
    amb σ_K = 0.1, σ_Π = 0.3.
    """

    print("===========================================")
    print("  Θ–Kerr horizon gate demo                 ")
    print("===========================================\n")

    steps = [
        "Definir el forat negre de Kerr",
        "Construir malla angular i calcular K(θ)",
        "Aplicar el model de ˜Π^2(θ)",
        "Avaluar G(˜K, ˜Π^2) i ⟨G⟩_H",
        "Resum final"
    ]
    total = len(steps)

    # Configuració d'exemple: quasi extremal
    bh = KerrBH(M=10.0, a=9.9)  # a* ≈ 0.99
    sigma_K = 0.1
    sigma_Pi = 0.3

    thetas = []
    G_vals = []
    Ktilde_vals = []
    G_avg = None

    for i, step in enumerate(steps, start=1):
        bar = "[" + "#" * i + "-" * (total - i) + "]"
        print(f"{bar} {step}...")

        if i == 1:
            pass
        elif i == 2:
            pass
        elif i == 3:
            pass
        elif i == 4:
            thetas, G_vals, Ktilde_vals, G_avg = horizon_gate_profile(
                bh,
                sigma_K=sigma_K,
                sigma_Pi=sigma_Pi,
                n_theta=400
            )
        elif i == 5:
            print("\n=========== SUMMARY ===========")
            print(f"M               = {bh.M:.3f}")
            print(f"a               = {bh.a:.3f}")
            print(f"a*              = {bh.a_star:.3f}")
            print(f"r_+             = {bh.r_plus:.6f}")
            print(f"sigma_K         = {sigma_K:.3f}")
            print(f"sigma_Pi        = {sigma_Pi:.3f}")
            print(f"<G>_H           = {G_avg:.5f}")
            print("================================\n")

    print("Demo Θ–Kerr (horizon gate) completada.\n")


# ----------------------------
# 5. Punt d'entrada
# ----------------------------

if __name__ == "__main__":
    demo()
