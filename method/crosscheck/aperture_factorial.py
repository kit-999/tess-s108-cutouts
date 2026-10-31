"""Why does a fixed 3x3 aperture carry the steps near BTJD 4305 and 4312? (Day 273, the desk's open follow-up after cross-check #1;
volunteered at 01:0x UTC Oct 31.)

On one fresh S108 TESSCut cutout per star, three apertures crossed with two backgrounds, nothing else changed:
  apertures: A3 (fixed 3x3 at the centre, his), A5 (fixed 5x5 at the centre), AC (connected, mine, as in method/batch.py)
  backgrounds: Bh (per-frame median of pixels at least 6 from the centre, his), Bb (per-frame median of the border, mine)
  flux = aperture sum - npix x background, over its median; frames with QUALITY == 0.
Step sizes from the 12-hour bins of crosscheck.py (55 from BTJD 4291.0):
  S1 = median(bins 23..27) - median(bins 29..33)   the step near 4305, in ppt
  S2 = median(bins 44..48) - median(bins 38..42)   the step near 4312, in ppt

Predictions, written and committed before the first run (the guess on file is the aperture):
  APERTURE: for each star, under BOTH backgrounds, |S1|(A3) - |S1|(AC) >= 2 ppt, |S1|(A5) lies between |S1|(AC) and |S1|(A3),
            and swapping the background moves |S1| by less than 1 ppt for every aperture.
  BACKGROUND (the rival): swapping the background moves |S1| by 2 ppt or more for A3, while A3 and AC differ by less than 1 ppt
            under the same background.
  Anything else is reported as MIXED, with the numbers. Crowding is not tested here: the same steps on unrelated stars point at
  the instrument, not at a neighbour.

  python -I aperture_factorial.py DLDIR TIC [TIC ...]
"""
import sys, os, json, shutil, zipfile, urllib.request, urllib.parse
import numpy as np
from astropy.io import fits

UA = {'User-Agent': 'kit999 (agent, fetching for a community vetting game)'}
OK_HOSTS = ('mast.stsci.edu', 'stpubdata.s3.amazonaws.com')
T0, W, NBINS, SECTOR = 4291.0, 0.5, 55, 108


def tic_row(tic):
    req = {"service": "Mast.Catalogs.Filtered.Tic", "format": "json", "params": {"columns": "ID,ra,dec,Tmag",
           "filters": [{"paramName": "ID", "values": [str(tic)]}]}}
    data = urllib.parse.urlencode({'request': json.dumps(req)}).encode()
    return json.loads(urllib.request.urlopen(urllib.request.Request('https://mast.stsci.edu/api/v0/invoke', data=data, headers=UA),
                                             timeout=120).read())['data'][0]


def cutout(ra, dec, dl):
    url = f'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={ra}&dec={dec}&y=15&x=15&sector={SECTOR}'
    r = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=900)
    if urllib.parse.urlparse(r.geturl()).netloc not in OK_HOSTS: raise RuntimeError('unexpected host')
    os.makedirs(dl, exist_ok=True)
    with open(f'{dl}/c.zip', 'wb') as f: shutil.copyfileobj(r, f)
    with zipfile.ZipFile(f'{dl}/c.zip') as z:
        n = [x for x in z.namelist() if x.endswith('.fits') and '/' not in x and '..' not in x]
        assert len(n) == 1, n
        z.extract(n[0], dl)
    return f'{dl}/{n[0]}'


def connected(med, border):
    ny, nx = med.shape
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
            yy, xx = y + dy, x + dx
            if 0 <= yy < ny and 0 <= xx < nx and mask[yy, xx] and not keep[yy, xx]: stack.append((yy, xx))
    return keep


def bins12(t, v):
    return [float(np.median(v[(t >= T0 + W * (k - 1)) & (t < T0 + W * k)])) if np.any((t >= T0 + W * (k - 1)) & (t < T0 + W * k))
            else None for k in range(1, NBINS + 1)]


def steps(b):
    pick = lambda lo, hi: float(np.median([b[k - 1] for k in range(lo, hi + 1) if b[k - 1] is not None]))
    return round((pick(23, 27) - pick(29, 33)) * 1000, 2), round((pick(44, 48) - pick(38, 42)) * 1000, 2)


def run(tic, dl):
    row = tic_row(tic)
    fp = cutout(row['ra'], row['dec'], f'{dl}/{tic}')
    d = fits.open(fp)[1].data
    t, fl, q = d['TIME'], np.asarray(d['FLUX'], float), d['QUALITY']
    good = (q == 0) & np.isfinite(t) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
    ny, nx = fl.shape[1:]; c = ny // 2
    yy, xx = np.indices((ny, nx))
    med = np.nanmedian(fl[good], axis=0)
    border = np.zeros((ny, nx), bool); border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    aps = {'A3': (abs(yy - c) <= 1) & (abs(xx - c) <= 1), 'A5': (abs(yy - c) <= 2) & (abs(xx - c) <= 2), 'AC': connected(med, border)}
    bgs = {'Bh': np.hypot(yy - c, xx - c) >= 6, 'Bb': border}
    out = {'tic': tic, 'tmag': row.get('Tmag'), 'frames': int(good.sum()), 'npix': {k: int(v.sum()) for k, v in aps.items()}}
    tt = t[good]
    for an, ap in aps.items():
        for bn, bg in bgs.items():
            f = np.nansum(fl[:, ap], axis=1) - ap.sum() * np.nanmedian(fl[:, bg], axis=1)
            v = f[good] / np.nanmedian(f[good])
            out[f'{an}-{bn}'] = dict(zip(('S1_ppt', 'S2_ppt'), steps(bins12(tt, v))))
    shutil.rmtree(f'{dl}/{tic}')
    S = lambda k: abs(out[k]['S1_ppt'])
    ap_ok = all(S(f'A3-{b}') - S(f'AC-{b}') >= 2 and min(S(f'AC-{b}'), S(f'A3-{b}')) <= S(f'A5-{b}') <= max(S(f'AC-{b}'), S(f'A3-{b}'))
                for b in bgs) and all(abs(S(f'{a}-Bh') - S(f'{a}-Bb')) < 1 for a in aps)
    bg_ok = abs(S('A3-Bh') - S('A3-Bb')) >= 2 and all(abs(S(f'A3-{b}') - S(f'AC-{b}')) < 1 for b in bgs)
    out['verdict'] = 'APERTURE' if ap_ok and not bg_ok else 'BACKGROUND' if bg_ok and not ap_ok else 'MIXED'
    return out


if __name__ == '__main__':
    dl = sys.argv[1]
    for tic in [int(x) for x in sys.argv[2:]]:
        print(json.dumps(run(tic, dl)))
