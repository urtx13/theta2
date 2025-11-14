# Θ–Kerr (2) Diagnostic Toolkit

This repository contains the **diagnostic and parameter–space exploration scripts** associated with the Θ–Kerr framework: an EFT–inspired, curvature–gated extension of General Relativity designed to activate only in the strong–field, high–spin regime of Kerr black holes.

The code implements:

- A **geometric gate** \( \mathcal{G}(K,\Pi^2) \) built from curvature invariants on Kerr,
- An **EFT scaling diagnostic** based on the dimensionless parameter  
  \(\chi = \varepsilon M^2 / \Lambda^2\),
- A **thermodynamic / EFT “triage protocol”** (PASS/REJECT) using a First Law residual,
- Simple **toy profiles** for the deformed metric in the strong–field region.

It is meant as a **companion numerical package** for the Θ–Kerr papers (framework + diagnostics), and as a starting point for future work on full metric solutions and observational applications (shadows, QNMs, EMRIs, etc.).

---

## 1. Installation

Clone the repository and create a virtual environment (optional but recommended):

```bash
git clone https://github.com/<your-username>/theta-kerr-diagnostics.git
cd theta-kerr-diagnostics

python -m venv sim_env
source sim_env/bin/activate  # Windows: sim_env\Scripts\activate

2.1 Core geometric gate and diagnostics
	•	theta_kerr_horizon_gate.py
Computes the geometric gate on Kerr at the event horizon:
	•	Defines a Kerr black hole with parameters ((M,a,a_*,r_+)),
	•	Evaluates the Kretschmann scalar (K(r_+,\theta)),
	•	Uses a phenomenological but controlled model for the parity invariant (\tilde{\Pi}^2(\theta)),
	•	Builds the smooth gate ( \mathcal{G}(\theta) = s_K(\tilde K),s_\Pi(\tilde\Pi^2) ),
	•	Computes the horizon–averaged gate
[
\langle \mathcal{G} \rangle_H = \frac{\int_0^\pi \mathcal{G}(\theta)\sin\theta,d\theta}
{\int_0^\pi \sin\theta,d\theta}.
]
	•	Prints a summary: (M, a, a_*, r_+, \sigma_K, \sigma_\Pi, \langle\mathcal{G}\rangle_H).
	•	theta_kerr_diagnostics.py
Implements the EFT + thermodynamic diagnostic pipeline:
	•	Defines a ThetaKerrParams container with
((M,\varepsilon,\Lambda,a_*,\sigma_K,\sigma_\Pi)),
	•	Computes the EFT coupling
(\chi = \varepsilon M^2/\Lambda^2),
	•	Calls the horizon–gate module to obtain (\langle\mathcal{G}\rangle_H),
	•	Estimates the photon–sphere shift:
[
\frac{\Delta r_{\rm ph}}{r_{\rm ph}} \approx A(a_),\chi,\langle\mathcal{G}\rangle_H,
\quad A(a_) \in [0.5,2] ,
]
and prints a conservative [min,max] interval,
	•	Evaluates a First Law residual
[
R = \frac{|dM - T_H,dS_{\rm Wald} - \Omega_H,dJ|}
{\max(|dM|, |T_H dS_{\rm Wald}|, \varepsilon M)} ,
]
	•	Applies a simple triage decision:
	•	Stage 1: gate coverage ((\langle\mathcal{G}\rangle_H \gtrsim 0.10)),
	•	Stage 2: First Law residual ((R \lesssim 0.02)),
	•	Returns PASS/REJECT and prints a human–readable summary.

⸻

2.2 Parameter scans
	•	theta_kerr_spin_scan.py
Spin scan at fixed coupling:
	•	Fixes typical EFT parameters, e.g.
(M = 10), (\varepsilon = 10^{-3}), (\Lambda = 100), (\sigma_K = 0.1), (\sigma_\Pi = 0.3),
	•	Scans over spins, e.g. (a_* = 0.3, 0.5, 0.7, 0.9, 0.99),
	•	For each (a_*):
	•	Computes (\langle\mathcal{G}\rangle_H(a_*)),
	•	Estimates (\Delta r_{\rm ph}/r_{\rm ph}),
	•	Evaluates the First Law residual,
	•	Applies the triage protocol (PASS/REJECT),
	•	Prints a table:
a*    |   chi      |  <G>_H   |  Δr_ph/r_ph(min)  |  Δr_ph/r_ph(max)  | triage

theta_kerr_chi_scan.py
Coupling scan at fixed near–extremal spin:
	•	Fixes (M = 10), (a_* = 0.99), (\Lambda = 100), (\sigma_K, \sigma_\Pi),
	•	Varies (\varepsilon), so that (\chi = \varepsilon M^2/\Lambda^2 \in {10^{-6},10^{-5},10^{-4}}),
	•	Computes (\langle\mathcal{G}\rangle_H) once (purely geometric for fixed (a_*)),
	•	For each (\chi):
	•	Estimates (\Delta r_{\rm ph}/r_{\rm ph}),
	•	Evaluates the First Law residual,
	•	Applies the triage protocol,
	•	Prints a table:
ε      |   chi      |  <G>_H   |  Δr_ph/r_ph(min)  |  Δr_ph/r_ph(max)  | triage

theta_kerr_param_map.py
2D map in spin–coupling space:
	•	Fixes (M = 10), (\Lambda = 100), (\sigma_K), (\sigma_\Pi),
	•	Loops over a set of spins (e.g. (a_* = 0.5, 0.7, 0.9, 0.99)),
	•	For each spin, computes (\langle\mathcal{G}\rangle_H(a_*)),
	•	Loops over a set of couplings (\chi) via (\varepsilon),
	•	For each pair ((a_*,\chi)):
	•	Estimates (\Delta r_{\rm ph}/r_{\rm ph}),
	•	Evaluates the First Law residual,
	•	Applies triage,
	•	Prints a 2D table:
a*   |   chi      |  <G>_H   |  Δr_ph/r_ph(min)  |  Δr_ph/r_ph(max)  | triage

theta_kerr_sigma_scan.py
Gate–sharpness scan at fixed near–extremal configuration:
	•	Fixes (M = 10), (a_* = 0.99), (\chi \approx 10^{-5}), (\Lambda = 100), (\sigma_K),
	•	Scans over the parity–gate width, e.g. (\sigma_\Pi \in {0.05, 0.10, 0.20, 0.30, 0.50}),
	•	For each (\sigma_\Pi):
	•	Computes (\langle\mathcal{G}\rangle_H(\sigma_\Pi)),
	•	Estimates (\Delta r_{\rm ph}/r_{\rm ph}),
	•	Applies triage (PASS/REJECT),
	•	Used to test robustness against gate sharpness (no fine–tuning in (\sigma_\Pi)).

⸻

2.3 Toy metric profile
	•	theta_kerr_metric_profile.py
(previously c_profile.py; the name can be adjusted)
	•	Constructs a 1D toy radial profile for the deformed metric in the equatorial plane:
[
g_{tt}^\Theta(r) = g_{tt}^{\rm GR}(r),
\big[1 + \alpha,\chi,\langle\mathcal{G}\rangle_H f(r)\big],
]
with (g_{tt}^{\rm GR}(r) = -(1 - 2M/r)), (\alpha = \mathcal{O}(1)), and (f(r)) a smooth bump centered at the horizon (r_+),
	•	Uses the same EFT amplitude (\chi,\langle\mathcal{G}\rangle_H) diagnosed by the other scripts,
	•	Prints a summary of the configuration and sample values of
(\delta g_{tt}/g_{tt}) at several radii (e.g. (r/r_+ = 1.01, 1.8, 2.6, 3.4, 5.0)),
	•	Illustrates that Θ–Kerr corrections are:
	•	small (typically (\delta g_{tt}/g_{tt} \sim 10^{-6}) near the horizon),
	•	localised in the strong–field region,
	•	consistent with the photon–sphere estimates (\Delta r_{\rm ph}/r_{\rm ph} \sim 10^{-6}).

This script is explicitly a toy metric solver for intuition and EFT sanity checks; it is not the full solution of the modified field equations.

⸻

3. Running the demos

From inside the virtual environment:
python theta_kerr_horizon_gate.py
python theta_kerr_diagnostics.py
python theta_kerr_spin_scan.py
python theta_kerr_chi_scan.py
python theta_kerr_param_map.py
python theta_kerr_sigma_scan.py
python theta_kerr_metric_profile.py

Each script prints a compact, human–readable summary to the terminal.
The typical benchmark configuration used throughout is:
	•	(M = 10),
	•	(a_* = 0.99) (near–extremal spin),
	•	(\varepsilon = 10^{-3}),
	•	(\Lambda = 100),
	•	(\chi = \varepsilon M^2 / \Lambda^2 = 10^{-5}),
	•	(\sigma_K = 0.1), (\sigma_\Pi = 0.3),

for which we find, for example:
	•	(\langle\mathcal{G}\rangle_H \approx 0.133),
	•	(\Delta r_{\rm ph}/r_{\rm ph} \sim (0.7–2.7)\times 10^{-6}),
	•	First Law residual (R \sim 10^{-2}),
	•	a PASS decision in the triage protocol.

⸻

4. Status and roadmap

This repository contains the diagnostic layer of the Θ–Kerr programme:
	•	geometric gate on Kerr,
	•	EFT scaling in terms of (\chi),
	•	horizon–averaged gate (\langle\mathcal{G}\rangle_H),
	•	parameter–space scans (spin, coupling, gate sharpness),
	•	toy metric profiles consistent with the EFT estimates.

It does not yet include:
	•	the full solution of the modified field equations for (\delta g_{\mu\nu}),
	•	detailed photon–sphere geodesic solvers on the deformed metric,
	•	full QNM / waveform computations or object–by–object constraints.

Those are intended for follow–up work (e.g. “Θ–Kerr black holes: metric perturbations and strong–field observables”).

⸻

5. How to cite

If you use this code in academic work, please cite the associated Θ–Kerr paper(s) and, if applicable, the Zenodo record of this repository.

Citation placeholder
(Replace with the actual reference once the paper and/or Zenodo DOI are available.)
