"""star27_chart.py: a print chart of star 27 (TIC 737546) for the town's planet league, from this bundle's own first-look curves
(Sector 108 and Sector 109's first half). Made 2026-11-02 at a town resident's ask, for the planet league's write-up.

What it does, in order:
  1. reads both sectors' curves (stars/TIC737546_s108_lc.csv, ..._s109_lc.csv) and their meta files (dates, method);
  2. times each primary eclipse that is covered on both sides: a parabola through the frames below 0.80 within 0.12 d of the
     predicted middle (from the box search's 3.0696 d), its vertex the time of minimum; an eclipse with fewer than 8 such frames, or
     whose vertex falls outside them, is left out and said so;
  3. fits T(n) = T0 + n P by least squares over those minima and prints P, T0, their standard errors and each eclipse's O-C;
  4. folds both sectors at P (phase 0 the primary), takes 0.005-phase bins (medians), and reads the primary and secondary depths
     from the binned curve;
  5. the control: folds again at P + 0.005 d and prints the binned primary depth there, which must come out shallower (smeared);
  6. draws two panels (Sector 108 against time; both sectors folded) at 300 dpi, with the credit and the method's limits on it.
Run with the astro venv: ~/claude/.venvs/astro/bin/python -I charts/star27_chart.py   (writes charts/star27_TIC737546.png)
"""
import csv, json, sys
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
STARS = HERE.parent / 'stars'
TIC = 737546
P_GUESS, T_GUESS = 3.0696, 4291.1117   # stars/TIC737546_s108_meta.json, box_search


def read(sector):
    t, f = [], []
    with open(STARS / f'TIC{TIC}_s{sector}_lc.csv', newline='') as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            t.append(float(row[0])); f.append(float(row[1]))
    meta = json.loads((STARS / f'TIC{TIC}_s{sector}_meta.json').read_text())
    return np.array(t), np.array(f), meta


def minima(t, f):
    out, skipped = [], []
    n0, n1 = int(np.floor((t.min() - T_GUESS) / P_GUESS)), int(np.ceil((t.max() - T_GUESS) / P_GUESS))
    for n in range(n0, n1 + 1):
        mid = T_GUESS + n * P_GUESS
        sel = (np.abs(t - mid) < 0.12) & (f < 0.80)
        if sel.sum() < 8:
            if (np.abs(t - mid) < 0.3).any(): skipped.append((n, int(sel.sum()), 'fewer than 8 frames below 0.80'))
            continue
        x, y = t[sel] - mid, f[sel]
        a, b, c = np.polyfit(x, y, 2)
        v = -b / (2 * a)
        if a <= 0 or not (x.min() < v < x.max()):
            skipped.append((n, int(sel.sum()), 'the vertex falls outside the frames')); continue
        out.append((n, mid + v, int(sel.sum())))
    return out, skipped


def binned(phase, f, width=0.005):
    edges = np.arange(-0.25, 0.75 + width / 2, width)
    centres, med = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (phase >= lo) & (phase < hi)
        if s.sum() >= 5: centres.append((lo + hi) / 2); med.append(np.median(f[s]))
    return np.array(centres), np.array(med)


def fold(t, P, T0):
    ph = ((t - T0) / P) % 1.0
    return np.where(ph >= 0.75, ph - 1.0, ph)


def main():
    t8, f8, m8 = read(108)
    t9, f9, m9 = read(109)
    mins, skipped = [], []
    for sec, t, f in ((108, t8, f8), (109, t9, f9)):
        a, b = minima(t, f); mins += [m + (sec,) for m in a]; skipped += b
    n = np.array([m[0] for m in mins], float); T = np.array([m[1] for m in mins])
    A = np.vstack([np.ones_like(n), n]).T
    coef, res, *_ = np.linalg.lstsq(A, T, rcond=None)
    T0, P = coef
    resid = T - A @ coef
    dof = len(T) - 2
    s2 = (resid @ resid) / dof
    cov = s2 * np.linalg.inv(A.T @ A)
    eT0, eP = np.sqrt(np.diag(cov))
    print(f'eclipses timed: {len(T)} (sector 108: {sum(1 for m in mins if m[3] == 108)}, sector 109: {sum(1 for m in mins if m[3] == 109)})')
    for (k, tm, nf, sec), r in zip(mins, resid):
        print(f'  n {k:3d}  sector {sec}  minimum BTJD {tm:.5f}  frames {nf:3d}  O-C {r * 1440:+6.2f} min')
    for k, nf, why in skipped: print(f'  n {k:3d}  left out ({nf} frames): {why}')
    print(f'period {P:.5f} +/- {eP:.5f} d ({eP * 86400:.0f} s); T0 BTJD {T0:.5f} +/- {eT0:.5f}; rms O-C {np.sqrt(s2) * 1440:.2f} min')
    t = np.concatenate([t8, t9]); f = np.concatenate([f8, f9])
    ph = fold(t, P, T0)
    bc, bm = binned(ph, f)
    base = np.median(f[(np.abs(ph - 0.25) < 0.05) | (np.abs(ph + 0.25) < 0.02) | (np.abs(ph - 0.70) < 0.04)])
    prim = base - bm[np.abs(bc) < 0.03].min()
    sec = base - bm[np.abs(bc - 0.5) < 0.04].min()
    print(f'binned depths, against the level near quarter phase ({base:.4f}): primary {prim * 100:.1f}%, secondary {sec * 100:.1f}%')
    bc2, bm2 = binned(fold(t, P + 0.005, T0), f)
    prim2 = base - bm2[np.abs(bc2) < 0.03].min()
    print(f'control: folded at P + 0.005 d the binned primary is {prim2 * 100:.1f}% deep (smeared), against {prim * 100:.1f}% at P')
    if not prim2 < prim: print('CONTROL FAILED: the wrong period folds as sharply as the fitted one'); return 3
    def scatter(Pt):
        ph_ = fold(t, Pt, T0); c_, m_ = binned(ph_, f); ins = np.abs(ph_) < 0.08
        return np.sqrt(np.mean((f[ins] - np.interp(ph_[ins], c_, m_)) ** 2))
    sc = {f'{Pt:.5f}': scatter(Pt) for Pt in (P, 3.06888, 3.0696)}
    print('in-eclipse scatter about the binned fold (|phase| < 0.08): ' + ', '.join(f'P {k}: {v * 100:.2f}%' for k, v in sc.items()))
    if min(sc, key=sc.get) != f'{P:.5f}': print('CHECK FAILED: an earlier period folds the eclipse more tightly than the fit'); return 3

    ink, accent, faint = '#1d2433', '#c2410c', '#8a94a6'
    plt.rcParams.update({'font.family': 'DejaVu Sans', 'font.size': 9, 'axes.edgecolor': ink, 'axes.labelcolor': ink,
                         'xtick.color': ink, 'ytick.color': ink, 'axes.spines.top': False, 'axes.spines.right': False})
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.0, 6.4), gridspec_kw={'height_ratios': [1, 1.25], 'hspace': 0.42})
    ax1.plot(t8 - t8.min(), f8, '.', ms=1.2, color=ink, alpha=0.55, rasterized=True)
    ax1.set_xlabel(f'days since {m8["date_obs"][:10]} (TESS Sector 108)')
    ax1.set_ylabel('brightness (relative)')
    ax1.set_title('Sector 108: four weeks of star 27', loc='left', fontsize=10, color=ink)
    ax2.plot(ph, f, '.', ms=1.0, color=faint, alpha=0.45, rasterized=True)
    ax2.plot(bc, bm, '-', lw=1.6, color=accent)
    ax2.set_xlim(-0.25, 0.75)
    ax2.set_xlabel(f'phase (one orbit = {P:.4f} days; 0 = the deep eclipse)')
    ax2.set_ylabel('brightness (relative)')
    ax2.set_title(f'Both sectors folded on the {P:.4f}-day period ({len(T)} deep eclipses timed)', loc='left', fontsize=10, color=ink)
    ax2.annotate(f'deep eclipse, {prim * 100:.0f}% fainter', xy=(0.0, base - prim), xytext=(0.08, base - prim + 0.02), color=ink, fontsize=8.5,
                 arrowprops={'arrowstyle': '-', 'color': ink, 'lw': 0.7})
    ax2.annotate(f'shallow eclipse at half phase, {sec * 100:.0f}% fainter', xy=(0.5, base - sec), xytext=(0.18, base - sec - 0.085), color=ink,
                 fontsize=8.5, arrowprops={'arrowstyle': '-', 'color': ink, 'lw': 0.7})
    fig.suptitle(f'Star 27 · TIC {TIC} · two stars eclipsing each other every {P:.3f} days', x=0.02, ha='left', fontsize=12,
                 fontweight='bold', color=ink)
    credit = (f'Data: NASA TESS, Sector 108 ({m8["date_obs"][:10]} to {m8["date_end"][:10]}) and the first half of Sector 109, full-frame-image '
              f'cutouts from MAST (TESSCut).\nPhotometry and chart: Kit (kit999.delve.town), github.com/kit-999/tess-s108-cutouts. '
              f'Depths against the brightness at quarter phase.\nA first look, not a vetted light curve: a simple aperture, no detrending; '
              f'the cutout curves carry small instrument steps shared across stars.')
    fig.text(0.02, 0.012, credit, fontsize=6.6, color=faint, ha='left', va='bottom', linespacing=1.35)
    fig.subplots_adjust(left=0.1, right=0.98, top=0.9, bottom=0.175)
    out = HERE / f'star27_TIC{TIC}.png'
    fig.savefig(out, dpi=300, facecolor='white')
    print(f'wrote {out.name}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
