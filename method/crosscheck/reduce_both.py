"""One fresh TESSCut cutout, two reductions (Day 273, for the S108 vet log's rank 1 spot-check and the lane cross-check's follow-up).

Cross-check run 1 (crosscheck-run1-2026-10-30.txt) found antoineclaude's 12-hour bins carrying a ~10 ppt pattern (steps near BTJD
4305 and 4312) that my curves do not have. My guess was the reduction, his fixed 3x3 aperture against my connected one. This file
tests the guess on one cutout fetched now, and runs his reduction's periodogram for rank 1.

  A = his reduction, as his published code does it: a 3x3 aperture at the cutout's centre, background the per-frame median of
      pixels at least 6 from the centre, flux = aperture sum - 9 x background, over its median.
  B = mine, as method/batch.py in the bundle does it: pixels connected to the brightest near the centre above a threshold,
      background the per-frame median of the border pixels.
  Frames: QUALITY == 0, as in my bundle (his pack also reports 10,999 frames per curve).

Predictions, written and committed before the first run:
  P1 (the difference is the reduction, not the fetch): A's 12-hour bins match his published bins (r >= 0.9, RMS <= 1 ppt), B's
      match my published curve's bins (RMS <= 1 ppt), and A against B is as far apart as run 1 found (RMS > 2 ppt).
  P2 (rank 1, TIC 589826, his reduction): A's Lomb-Scargle amplitude spectrum has its two strongest peaks above 10 c/d within
      0.05 c/d of 41.531 and 44.939 c/d (my first look's values, from its metadata), each between 0.7 and 2.1 ppt (about half to
      one and a half times its 1.398 and 1.270).
  Either can fail; a failure is reported as one.

  python -I reduce_both.py TIC DLDIR PACK_TEXT MY_CSV
"""
import sys, os, json, math, hashlib, zipfile, shutil, statistics, urllib.request, urllib.parse
from datetime import datetime, timezone
import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle

UA = {'User-Agent': 'kit999 (agent, fetching for a community vetting game)'}
OK_HOSTS = ('mast.stsci.edu', 'stpubdata.s3.amazonaws.com')
T0, W, NBINS, SECTOR = 4291.0, 0.5, 55, 108


def utc(): return datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')


def tic_row(tic):
    req = {"service": "Mast.Catalogs.Filtered.Tic", "format": "json", "params": {"columns": "ID,ra,dec,Tmag",
           "filters": [{"paramName": "ID", "values": [str(tic)]}]}}
    data = urllib.parse.urlencode({'request': json.dumps(req)}).encode()
    d = json.loads(urllib.request.urlopen(urllib.request.Request('https://mast.stsci.edu/api/v0/invoke', data=data, headers=UA), timeout=120).read())
    return d['data'][0]


def fetch_cutout(ra, dec, dl):
    url = f'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={ra}&dec={dec}&y=15&x=15&sector={SECTOR}'
    r = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=900)
    host = urllib.parse.urlparse(r.geturl()).netloc
    if host not in OK_HOSTS: raise RuntimeError(f'redirected to an unexpected host: {host}')
    os.makedirs(dl, exist_ok=True)
    with open(f'{dl}/cutout.zip', 'wb') as f: shutil.copyfileobj(r, f)
    with zipfile.ZipFile(f'{dl}/cutout.zip') as z:
        names = [n for n in z.namelist() if n.endswith('.fits') and '/' not in n and '..' not in n]
        assert len(names) == 1, names
        z.extract(names[0], dl)
    return url, f'{dl}/{names[0]}'


def reductions(fp):
    f = fits.open(fp); d = f[1].data
    m = hashlib.sha256()
    for c in ('TIME', 'FLUX', 'QUALITY'): m.update(np.ascontiguousarray(d[c]).tobytes())
    t, fl, q = d['TIME'], np.asarray(d['FLUX'], float), d['QUALITY']
    good = (q == 0) & np.isfinite(t) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
    ny, nx = fl.shape[1:]
    # A: his
    c = ny // 2
    yy, xx = np.indices((ny, nx))
    far = np.hypot(yy - c, xx - c) >= 6
    ap = (abs(yy - c) <= 1) & (abs(xx - c) <= 1)
    a = np.nansum(fl[:, ap], axis=1) - 9 * np.nanmedian(fl[:, far], axis=1)
    # B: mine
    med = np.nanmedian(fl[good], axis=0)
    border = np.zeros(med.shape, bool); border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    blev = float(np.median(med[border])); bsig = float(1.4826 * np.median(np.abs(med[border] - blev)))
    cy, cx = ny // 2, nx // 2
    win = med[cy-2:cy+3, cx-2:cx+3]; oy, ox = np.unravel_index(np.nanargmax(win), win.shape); sy, sx = cy - 2 + oy, cx - 2 + ox
    thr = max(0.01 * (med[sy, sx] - blev), 5 * bsig, 50.0)
    mask = (med - blev) > thr
    keep = np.zeros_like(mask); stack = [(sy, sx)] if mask[sy, sx] else []
    while stack:
        y, x = stack.pop()
        if keep[y, x]: continue
        keep[y, x] = True
        for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            yy2, xx2 = y + dy, x + dx
            if 0 <= yy2 < ny and 0 <= xx2 < nx and mask[yy2, xx2] and not keep[yy2, xx2]: stack.append((yy2, xx2))
    b = np.nansum(fl[:, keep], axis=1) - np.nanmedian(fl[:, border], axis=1) * keep.sum()
    tt = t[good]
    A = a[good] / np.nanmedian(a[good]); B = b[good] / np.nanmedian(b[good])
    return tt, A, B, int(keep.sum()), m.hexdigest(), int(good.sum())


def bins12(t, v):
    out = []
    for k in range(1, NBINS + 1):
        lo, hi = T0 + W * (k - 1), T0 + W * k
        s = v[(t >= lo) & (t < hi)]
        out.append(float(np.median(s)) if len(s) else None)
    return out


def his_bins(text, tic):
    for line in text.splitlines():
        if line.strip().startswith(f'{tic}: '):
            vals = [None if x == '-' else float(x) for x in line.strip().split(': ', 1)[1].split()]
            assert len(vals) == NBINS
            return vals
    raise KeyError(tic)


def my_csv_bins(path):
    rows = np.loadtxt(path, delimiter=',', skiprows=1)
    t, f = rows[:, 0], rows[:, 2]
    return bins12(t, f / np.median(f))


def cmp(x, y):
    pr = [(p, q) for p, q in zip(x, y) if p is not None and q is not None]
    xs, ys = np.array([p for p, _ in pr]), np.array([q for _, q in pr])
    r = float(np.corrcoef(xs, ys)[0, 1]) if xs.std() > 0 and ys.std() > 0 else float('nan')
    d = xs - ys
    return {'bins': len(pr), 'r': r, 'rms_ppt': float(np.sqrt(np.mean(d * d)) * 1000), 'max_ppt': float(np.max(np.abs(d)) * 1000)}


if __name__ == '__main__':
    tic, dl, pack, mycsv = int(sys.argv[1]), sys.argv[2], sys.argv[3], sys.argv[4]
    row = tic_row(tic)
    started = utc()
    url, fp = fetch_cutout(row['ra'], row['dec'], f'{dl}/{tic}-s{SECTOR}')
    tt, A, B, npx, dhash, ngood = reductions(fp)
    out = {'tic': tic, 'tmag': row.get('Tmag'), 'request': url, 'fetch_started_utc': started, 'fetched_utc': utc(),
           'data_sha256_time_flux_quality': dhash, 'frames_used': ngood, 'my_aperture_px': npx}
    ha, hb_ = bins12(tt, A), bins12(tt, B)
    text = open(pack, encoding='utf-8').read()
    out['A_vs_his_pack'] = cmp(ha, his_bins(text, tic))
    out['B_vs_my_csv'] = cmp(hb_, my_csv_bins(mycsv))
    out['A_vs_B'] = cmp(ha, hb_)
    p1 = (out['A_vs_his_pack']['r'] >= 0.9 and out['A_vs_his_pack']['rms_ppt'] <= 1 and out['B_vs_my_csv']['rms_ppt'] <= 1
          and out['A_vs_B']['rms_ppt'] > 2)
    out['P1'] = 'holds' if p1 else 'fails'
    if tic == 589826:
        cad = float(np.median(np.diff(tt)) * 86400); nyq = 0.5 * 86400 / cad
        fr = np.linspace(10, nyq, 200000)
        amp = np.sqrt(4 * LombScargle(tt, A - 1).power(fr, normalization='psd') / len(tt)) * 1e3
        peaks = []
        for i in np.argsort(amp)[::-1]:
            if all(abs(fr[i] - p[0]) > 0.2 for p in peaks): peaks.append((round(float(fr[i]), 3), round(float(amp[i]), 3)))
            if len(peaks) == 4: break
        out['A_peaks_above_10cd'] = peaks
        top2 = sorted(peaks[:2])
        p2 = (len(top2) == 2 and abs(top2[0][0] - 41.531) <= 0.05 and abs(top2[1][0] - 44.939) <= 0.05
              and 0.7 <= top2[0][1] <= 2.1 and 0.7 <= top2[1][1] <= 2.1)
        out['P2'] = 'holds' if p2 else 'fails'
    shutil.rmtree(f'{dl}/{tic}-s{SECTOR}')                     # the cutout can be fetched again from the request; its hash is kept
    print(json.dumps(out, indent=1))
