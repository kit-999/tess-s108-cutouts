"""Lane cross-check for the S108 vet log (Day 273, assigned by kimik3 at 00:40 UTC Oct 31): antoineclaude's published 12-hour
binned values against my own light curves, on ranks drawn before looking with random.Random(20261031).sample(range(1, 31), 3)
= 8, 10, 25. The item under test is his reduction (a fixed 3x3 aperture, background from pixels at least 6 from the centre,
lightkurve's default quality mask, flux over its median), not mine.

Rule, fixed before any number was computed (this file is committed before its first run):
  my curve is re-normalized the way his is (flux_e_per_s over its own median), binned into his 55 bins (bin k covers BTJD
  4291.0 + 0.5(k-1) to 4291.0 + 0.5k, half-open, median of each), and compared with his three-decimal values on the bins
  both have.
  - consistent: Pearson r >= 0.9 over those bins AND RMS difference <= 2 ppt
  - shape agrees, level differs: r >= 0.9 and RMS > 2 ppt (two apertures can see the same light with different dilution)
  - inconsistent: r < 0.9
  The empty bins (his '-') must coincide with mine having no frames there; a mismatch is reported, not scored.

  python3 crosscheck.py PACK_TEXT CSV_DIR [ranks...]
"""
import csv, hashlib, math, re, statistics, sys
from pathlib import Path

MANIFEST = [589826, 589906, 590274, 592638, 593228, 593230, 595279, 607907, 608932, 610232,
            612908, 614129, 614556, 615133, 627296, 664269, 672060, 672624, 675445, 682256,
            682491, 683071, 709603, 710187, 711970, 737327, 737546, 740069, 767526, 767943]
T0, W, NBINS = 4291.0, 0.5, 55


def his_bins(pack_text):
    out = {}
    for line in pack_text.splitlines():
        m = re.match(r'^(\d{6,}): (.+)$', line.strip())
        if m:
            vals = [None if v == '-' else float(v) for v in m.group(2).split()]
            assert len(vals) == NBINS, (m.group(1), len(vals))
            out[int(m.group(1))] = vals
    return out


def my_bins(csv_path):
    t, f = [], []
    with open(csv_path, newline='') as fh:
        r = csv.reader(fh); next(r)
        for row in r:
            t.append(float(row[0])); f.append(float(row[2]))
    med = statistics.median(f)
    norm = [x / med for x in f]
    bins = []
    for k in range(1, NBINS + 1):
        lo, hi = T0 + W * (k - 1), T0 + W * k
        v = [n for tt, n in zip(t, norm) if lo <= tt < hi]
        bins.append(statistics.median(v) if v else None)
    return bins, len(t)


def compare(h, m):
    pairs = [(a, b) for a, b in zip(h, m) if a is not None and b is not None]
    gaps_his = [i + 1 for i, a in enumerate(h) if a is None]
    gaps_mine = [i + 1 for i, b in enumerate(m) if b is None]
    xs, ys = [p[1] for p in pairs], [p[0] for p in pairs]
    mx, my = statistics.fmean(xs), statistics.fmean(ys)
    cov = sum((x - mx) * (y - my) for x, y in zip(xs, ys))
    r = cov / math.sqrt(sum((x - mx) ** 2 for x in xs) * sum((y - my) ** 2 for y in ys))
    d = [x - y for x, y in zip(xs, ys)]
    rms = math.sqrt(sum(v * v for v in d) / len(d)) * 1000
    mad = max(abs(v) for v in d) * 1000
    verdict = ('consistent' if r >= 0.9 and rms <= 2 else 'shape agrees, level differs' if r >= 0.9 else 'inconsistent')
    return {'bins': len(pairs), 'r': r, 'rms_ppt': rms, 'max_ppt': mad, 'gaps_his': gaps_his, 'gaps_mine': gaps_mine,
            'verdict': verdict}


if __name__ == '__main__':
    pack, csvdir = Path(sys.argv[1]), Path(sys.argv[2])
    ranks = [int(x) for x in sys.argv[3:]] or [8, 10, 25]
    text = pack.read_text(encoding='utf-8')
    print(f'pack sha256 {hashlib.sha256(pack.read_bytes()).hexdigest()} ({len(text)} chars)')
    hb = his_bins(text)
    assert len(hb) == 30 and set(hb) == set(MANIFEST), sorted(set(MANIFEST) ^ set(hb))
    for rank in ranks:
        tic = MANIFEST[rank - 1]
        mb, n = my_bins(csvdir / f'TIC{tic}_s108_lc.csv')
        c = compare(hb[tic], mb)
        print(f"rank {rank} TIC {tic}: {c['bins']} bins compared ({n} frames of mine) | r = {c['r']:.3f} | RMS diff "
              f"{c['rms_ppt']:.2f} ppt | max diff {c['max_ppt']:.2f} ppt | his empty bins {c['gaps_his']} | mine "
              f"{c['gaps_mine']} | {c['verdict']}")
        print('   his :', ' '.join('-' if v is None else f'{v:.3f}' for v in hb[tic]))
        print('   mine:', ' '.join('-' if v is None else f'{v:.3f}' for v in mb))
