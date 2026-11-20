"""
part6e_window_processes.py -- monthly shock processes for SVAR A's window, by filter.

The model side of SVAR A (ags_svar_sim.jl) needs, for each filter, monthly AR(1) processes whose
simulated quarterly averages, filtered as the data are, reproduce the data's AR(1) statistics on
SVAR A's window. part6b already does this for HP (section 9, columns bprime_*). This program
adds the Hamilton filter (h = 8, p = 4), with part6b's own functions:
    z, s   Hamilton cycles over 1951-2019 (one-sided, so the first 11 quarters are lost); target
           = AR(1) on 1995Q2-2019Q4, SVAR A's Hamilton window; the simulation filters the same
           276 quarters and fits on the last 99.
    delta  part6b's Hamilton-matched process (section 8), fitted on delta's own Hamilton window
           (1995Q2-2019Q4) from the (z, delta) VAR, as under HP.
Under Hamilton no stationary monthly AR(1) may be persistent enough for z; the cap (0.999) is
then reported as not attained.

Output: data/results/shock_calibration_window.csv, one row per series (z, delta, s), columns
hp_rho_m, hp_sigma_m (= part6b's bprime) and ham_rho_m, ham_sigma_m, ham_attained, and the
targets. Read by Programs AGS/ags_svar_sim.jl (shock_processes, prefix "hp" or "ham").
"""

import numpy as np
import pandas as pd

import part6b_shock_persistence_cyclical as p6b
import q3_data_svar as q3

OUT_CSV   = p6b.RESULTS_DIR / "shock_calibration_window.csv"
HAM_START = "1995Q2"


if __name__ == "__main__":
    calib = pd.read_csv(p6b.RESULTS_DIR / "shock_calibration.csv").set_index("series")
    cyc = q3.build_cycles("hamilton")
    rng = np.random.default_rng(p6b.MATCH_SEED)
    lost = p6b.HAMILTON_H + p6b.HAMILTON_P - 1
    rows = []
    for name in ("z", "delta", "s"):
        c = calib.loc[name]
        row = {"series": name, "hp_rho_m": c["bprime_rho_m"], "hp_sigma_m": c["bprime_sigma_m"]}
        if name == "delta":
            row.update({"ham_rho_m": c["hamilton_matched_rho_m"], "ham_sigma_m": c["hamilton_matched_sigma_m"],
                        "ham_attained": True, "ham_rho_Q": c["hamilton_rho_Q"],
                        "ham_sigma_Q": c["hamilton_sigma_Q"], "ham_target": "part6b section 8"})
        else:
            full = cyc[name].dropna()
            win = full[full.index >= HAM_START]
            fit = p6b.fit_ar1(win)
            nq, n_keep = len(full) + lost, len(win)
            rho_m, sigma_m, eps, C, attained = p6b.match_process(fit["rho_Q"], fit["sigma_Q"], nq, rng,
                                                                 method="hamilton", n_keep=n_keep)
            row.update({"ham_rho_m": rho_m, "ham_sigma_m": sigma_m, "ham_attained": attained,
                        "ham_rho_Q": fit["rho_Q"], "ham_sigma_Q": fit["sigma_Q"],
                        "ham_target": f"AR(1) on {win.index[0]}-{win.index[-1]} (n {n_keep})"})
        rows.append(row)
    out = pd.DataFrame(rows).set_index("series")
    print(out.round(4).to_string())
    out.to_csv(OUT_CSV)
    print(f"\nSaved: {OUT_CSV}")
