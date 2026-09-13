from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
import sys
import time
import warnings

HERE = Path(__file__).resolve().parent
PROJECT = HERE
sys.path.insert(0, str(PROJECT))
from solve import Grid, SEED, integrate, read_inputs
import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import eigh_tridiagonal
from scipy.special import j0, j1, jn_zeros

R0, TF, D, C0, AMP = 0.02, 3600.0, 1e-8, 1.0, 0.3
BETA = 1.0 / (2 * TF)
LAM = float(jn_zeros(1, 1)[0])
TIMES = np.linspace(0, TF, 129)
SPACE_LIMIT, TIME_LIMIT, MASS_LIMIT = 1e-5, 5e-6, 1e-10


def radius(t):
    return R0 * (1 - BETA * np.asarray(t))


def clock_integral(t):
    t = np.asarray(t)
    return t / (R0**2 * (1 - BETA * t))


def exact(x, times=TIMES):
    return C0 + AMP * j0(LAM * np.asarray(x)[:, None]) * np.exp(
        -LAM**2 * D * clock_integral(times)[None, :]
    )


def semidiscrete_reference(grid, initial, times=TIMES):
    lap = grid.jac_block(initial, 1.0, 0.0, 1.0, 0.0, same=True)
    sw = np.sqrt(grid.w)
    upper = lap.diagonal(1) * sw[:-1] / sw[1:]
    lower = lap.diagonal(-1) * sw[1:] / sw[:-1]
    assert np.allclose(upper, lower, rtol=1e-13, atol=1e-10)
    vals, vecs = eigh_tridiagonal(lap.diagonal(), upper)
    assert abs(vals[-1]) < 1e-7 and vals[-2] < 0
    vals[-1] = 0.0
    coeff = vecs.T @ (sw * initial)
    answer = (vecs @ (coeff[:, None] * np.exp(
        vals[:, None] * D * clock_integral(times)[None, :]
    ))) / sw[:, None]
    assert np.max(abs(answer[:, 0] - initial)) < 1e-11
    return answer


def solve_case(n, *, graded=True, rtol=1e-11, max_step=TF/256,
               mutant=False, temporal_reference=False):
    grid = Grid(n, graded=graded)
    initial = exact(grid.x, np.array([0.0]))[:, 0]
    lap = grid.jac_block(initial, D, 0.0, 1.0, 0.0, same=True)

    def r_used(t):
        return R0 if mutant else float(radius(t))

    def rhs(t, c):
        return grid.rate(c, D, r_used(t), 0.0, 0.0)

    result = solve_ivp(rhs, (0, TF), initial, method="BDF", rtol=rtol,
                       atol=rtol*0.01, max_step=max_step,
                       jac=lambda t, c: lap/r_used(t)**2, dense_output=True)
    if not result.success:
        raise RuntimeError(result.message)
    values = result.sol(TIMES)
    reference = exact(grid.x)
    mass0 = float(2*grid.w @ initial)
    mass_drift = max(float(np.max(abs(2*grid.w @ values-mass0))),
                     float(np.max(abs(2*grid.w @ result.y-mass0))))
    continuum_error = float(np.max(abs(values-reference)))
    temporal_error = None
    if temporal_reference:
        temporal_error = float(np.max(abs(values-semidiscrete_reference(grid, initial))))
    row = dict(n=n, graded=graded, rtol=rtol, atol=rtol*0.01,
               max_step_s=max_step, accepted_steps=len(result.t)-1,
               actual_max_step_s=float(np.diff(result.t).max()),
               continuum_max_error=continuum_error, temporal_max_error=temporal_error,
               mass_drift_max=mass_drift, initial_mean_quadrature_error=mass0-C0,
               numerical_center_final=float(values[0,-1]),
               exact_center_final=float(reference[0,-1]),
               mutant=mutant, nfev=result.nfev)
    return row


def regression_checks():
    arrays, snapshots = read_inputs()
    rows = []
    for q in (1, 3, 4):
        run = integrate(q, 80, arrays, horizon=600, stop_at_threshold=False)
        times = np.unique(np.concatenate([np.arange(0, 601)]
                                          + [s.t for s in run.segments]))
        values = run.values(times)
        m = run.model.grid.n+1
        residual = float(np.max(abs(2*run.model.grid.w@values[m:2*m]
                                    +values[-1]-2.55)))
        closed = run.model.grid.rate(np.full(m, 2.55), 1e-9, 0.015, 0, 0)
        rows.append(dict(question=q, mass_balance_max=residual,
                         closed_uniform_rate_max=float(np.max(abs(closed))),
                         radius_start_m=float(run.model.radius(0)),
                         radius_end_m=float(run.model.radius(600)),
                         passed=bool(residual < 1e-7 and np.max(abs(closed)) == 0)))
    return rows, snapshots


def run_suite():
    start = time.perf_counter()
    space = [solve_case(n) for n in (40, 80, 160, 320, 640)]
    temporal = [solve_case(160, graded=False, rtol=tol, max_step=step,
                           temporal_reference=True)
                for tol, step in ((1e-3, TF/8), (1e-5, TF/16), (1e-7, TF/32))]
    mutant = solve_case(640, mutant=True)
    regressions, snapshots = regression_checks()
    se = [r["continuum_max_error"] for r in space]
    te = [r["temporal_max_error"] for r in temporal]
    checks = dict(
        radius_changes=radius(TF) != radius(0),
        nonuniform_initial=np.ptp(exact(np.linspace(0,1,101),np.array([0]))) > 0.1,
        neumann_eigenvalue=abs(j1(LAM)) < 1e-12,
        correct_operator=se[-1] <= SPACE_LIMIT,
        space_errors_decrease=all(b < a for a,b in zip(se,se[1:])),
        time_errors_decrease=all(b < a for a,b in zip(te,te[1:])),
        time_grid_refines=all(b["accepted_steps"] > a["accepted_steps"]
                             for a,b in zip(temporal,temporal[1:])),
        finest_temporal_error=te[-1] <= TIME_LIMIT,
        mutant_rejected=mutant["continuum_max_error"] >= 100*SPACE_LIMIT,
        closed_mass_preserved=all(r["mass_drift_max"] <= MASS_LIMIT
                                  for r in space+temporal+[mutant]),
        original_mass_regression=all(r["passed"] for r in regressions),
    )
    checks = {k:bool(v) for k,v in checks.items()}
    return dict(status="PASS" if all(checks.values()) else "FAIL", checks=checks,
                parameters=dict(R0_m=R0, final_radius_m=float(radius(TF)), tf_s=TF,
                                D_m2_s=D, C0=C0, amplitude=AMP, beta_per_s=BETA,
                                eigenvalue=LAM, comparison_times_s=TIMES.tolist(),
                                spatial_error_limit=SPACE_LIMIT, time_error_limit=TIME_LIMIT,
                                closed_mass_limit=MASS_LIMIT),
                spatial=space, temporal=temporal, mutant=mutant,
                mass_regressions=regressions, inputs=snapshots, seed=SEED,
                production_solver_sha256=hashlib.sha256((PROJECT/"solve.py").read_bytes()).hexdigest(),
                benchmark_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                runtime_s=time.perf_counter()-start)


def main():
    warnings.filterwarnings("error", category=RuntimeWarning)
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path)
    parser.add_argument("--mutant-only", action="store_true")
    args = parser.parse_args()
    if args.mutant_only:
        row = solve_case(640, mutant=True)
        print(json.dumps(row, ensure_ascii=False, indent=2))
        return 1 if row["continuum_max_error"] > SPACE_LIMIT else 0
    report = run_suite()
    if args.output_dir:
        out = args.output_dir.resolve()
        if PROJECT not in out.parents:
            raise ValueError("Outputs must stay below A题")
        out.mkdir(parents=True, exist_ok=True)
        (out/"moving_domain_benchmark.json").write_text(
            json.dumps(report, ensure_ascii=False, indent=2), encoding="utf-8")
        with (out/"moving_domain_convergence.csv").open("w", newline="", encoding="utf-8") as stream:
            rows = [dict(kind=k, **r) for k, group in
                    (("space",report["spatial"]),("time",report["temporal"]),
                     ("mutant",[report["mutant"]])) for r in group]
            writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
            writer.writeheader(); writer.writerows(rows)
    print(json.dumps(report, ensure_ascii=False, indent=2))
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
