"""TESS sector cutouts and first-look light curves for a list of TIC IDs (a town game's manifest, Oct 26, 2026).
For each TIC: position from the MAST TIC catalog, sector coverage from TESSCut, a 15x15 cutout from TESSCut (followed only to
MAST's own hosts), then simple aperture photometry, a periodogram and a box search for repeating dips. A first look, not a vet.
Usage: python -I batch.py SECTOR OUTDIR DLDIR TIC [TIC ...]"""
import sys, os, json, time, hashlib, zipfile, shutil, urllib.request, urllib.parse
from datetime import datetime, timezone
import numpy as np
from astropy.io import fits
from astropy.timeseries import LombScargle, BoxLeastSquares
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt

UA = {'User-Agent': 'kit999 (agent, fetching for a community vetting game)'}
OK_HOSTS = ('mast.stsci.edu', 'stpubdata.s3.amazonaws.com')
SECTOR, OUT, DL = int(sys.argv[1]), sys.argv[2], sys.argv[3]
TICS = [int(x) for x in sys.argv[4:]]
os.makedirs(OUT, exist_ok=True); os.makedirs(DL, exist_ok=True)
log = open(f'{OUT}/batch-s{SECTOR}.log', 'a')
def say(*a):
    line = ' '.join(str(x) for x in a); print(line, flush=True); log.write(line + '\n'); log.flush()
def utc(): return datetime.now(timezone.utc).strftime('%Y-%m-%dT%H:%M:%SZ')

def tic_rows(ids):
    req = {"service": "Mast.Catalogs.Filtered.Tic", "format": "json", "params": {"columns": "ID,ra,dec,Tmag,Teff,rad,contratio,HIP,GAIA",
           "filters": [{"paramName": "ID", "values": [str(i) for i in ids]}]}}
    data = urllib.parse.urlencode({'request': json.dumps(req)}).encode()
    d = json.loads(urllib.request.urlopen(urllib.request.Request('https://mast.stsci.edu/api/v0/invoke', data=data, headers=UA), timeout=120).read())
    return {int(r['ID']): r for r in d.get('data') or []}

def sectors(ra, dec):
    u = f'https://mast.stsci.edu/tesscut/api/v0.1/sector?ra={ra}&dec={dec}&radius=0m'
    return {int(x['sector']): (int(x['camera']), int(x['ccd'])) for x in json.loads(urllib.request.urlopen(urllib.request.Request(u, headers=UA), timeout=120).read()).get('results', [])}

def fetch(url, dest):
    r = urllib.request.urlopen(urllib.request.Request(url, headers=UA), timeout=900)
    host = urllib.parse.urlparse(r.geturl()).netloc
    if host not in OK_HOSTS: raise RuntimeError(f'redirected to an unexpected host: {host}')
    with open(dest, 'wb') as f: shutil.copyfileobj(r, f)
    return r.geturl()

def analyse(fpath, tic, meta):
    f = fits.open(fpath); d = f[1].data; h0 = f[0].header
    m = hashlib.sha256()                                   # survives a re-fetch: TESSCut stamps each copy's checksum cards
    for c in ('TIME', 'FLUX', 'QUALITY'): m.update(np.ascontiguousarray(d[c]).tobytes())
    meta['data_sha256_time_flux_quality'] = m.hexdigest()
    t, fl, q = d['TIME'], d['FLUX'], d['QUALITY']
    good = (q == 0) & np.isfinite(t) & np.all(np.isfinite(fl.reshape(len(t), -1)), axis=1)
    med = np.nanmedian(fl[good], axis=0); ny, nx = med.shape
    border = np.zeros(med.shape, bool); border[0, :] = border[-1, :] = border[:, 0] = border[:, -1] = True
    blev = float(np.median(med[border])); bsig = float(1.4826 * np.median(np.abs(med[border] - blev)))
    cy, cx = ny // 2, nx // 2
    win = med[cy-2:cy+3, cx-2:cx+3]; oy, ox = np.unravel_index(np.nanargmax(win), win.shape); sy, sx = cy - 2 + oy, cx - 2 + ox
    peak = med[sy, sx] - blev
    thr = max(0.01 * peak, 5 * bsig, 50.0)
    mask = (med - blev) > thr
    keep = np.zeros_like(mask); stack = [(sy, sx)] if mask[sy, sx] else []
    while stack:
        y, x = stack.pop()
        if keep[y, x]: continue
        keep[y, x] = True
        for dy, dx in ((1, 0), (-1, 0), (0, 1), (0, -1)):
            yy, xx = y + dy, x + dx
            if 0 <= yy < ny and 0 <= xx < nx and mask[yy, xx] and not keep[yy, xx]: stack.append((yy, xx))
    if keep.sum() == 0: raise RuntimeError('no aperture: the star is not above the background in the cutout')
    bkg = np.nanmedian(fl[:, border], axis=1)
    raw = np.nansum(fl[:, keep], axis=1) - bkg * keep.sum()
    tt, ff = t[good], raw[good]; norm = ff / np.nanmedian(ff)
    cad = float(np.median(np.diff(tt)) * 86400)
    rob = float(1.4826 * np.median(np.abs(norm - np.median(norm))) * 1e6)
    edges = np.arange(tt.min(), tt.max() + 1/24, 1/24); idx = np.digitize(tt, edges)
    bins = [(edges[i-1] + 1/48, float(np.mean(norm[idx == i]))) for i in np.unique(idx) if (idx == i).sum() >= 10]
    bt = np.array([b[0] for b in bins]); bv = np.array([b[1] for b in bins])
    brob = float(1.4826 * np.median(np.abs(bv - np.median(bv))) * 1e6); lo = int(np.argmin(bv))
    # periodogram: strongest sinusoids from 0.1 c/d to the Nyquist frequency
    nyq = 0.5 * 86400 / cad; fr = np.linspace(0.1, nyq, 150000)
    amp = np.sqrt(4 * LombScargle(tt, norm - 1).power(fr, normalization='psd') / len(tt)) * 1e3
    peaks = []
    for i in np.argsort(amp)[::-1]:
        if all(abs(fr[i] - p[0]) > 0.2 for p in peaks): peaks.append((round(float(fr[i]), 3), round(float(amp[i]), 3)))
        if len(peaks) == 3: break
    # box search for repeating dips, on 20-minute bins
    e2 = np.arange(tt.min(), tt.max() + 1/72, 1/72); i2 = np.digitize(tt, e2)
    b2 = [(e2[i-1] + 1/144, float(np.mean(norm[i2 == i]))) for i in np.unique(i2) if (i2 == i).sum() >= 3]
    bt2 = np.array([b[0] for b in b2]); bv2 = np.array([b[1] for b in b2])
    dy2 = np.full(len(bv2), 1.4826 * np.median(np.abs(bv2 - np.median(bv2))))   # one robust scatter for every 20-minute bin
    bls = BoxLeastSquares(bt2, bv2, dy=dy2)
    periods = np.exp(np.linspace(np.log(0.5), np.log((tt.max() - tt.min()) / 2), 6000))
    res = bls.power(periods, [0.04, 0.08, 0.125, 0.2], objective='snr')
    k = int(np.argmax(res.power)); sde = float((res.power[k] - np.mean(res.power)) / np.std(res.power))
    st = bls.compute_stats(res.period[k], res.duration[k], res.transit_time[k])
    ntr = int(np.sum(st['per_transit_count'] > 0))
    fig, ax = plt.subplots(3, 1, figsize=(10, 8.5))
    ax[0].plot(tt, norm, ',', color='0.45'); ax[0].plot(bt, bv, '.', color='C1', ms=3); ax[0].set_ylabel('normalized flux')
    ax[0].set_title(f'TIC {tic} · TESS sector {SECTOR} · TESSCut FFI cutout, simple aperture ({int(keep.sum())} px) · first look, not a vet')
    ax[0].set_xlabel('time (BTJD = BJD − 2457000)')
    ax[1].plot(fr, amp, lw=0.5, color='C0'); ax[1].set_xlabel('frequency (cycles/day)'); ax[1].set_ylabel('amplitude (ppt)')
    ph = ((bt2 - res.transit_time[k] + 0.5 * res.period[k]) % res.period[k]) / res.period[k] - 0.5
    ax[2].plot(ph * res.period[k] * 24, (bv2 - 1) * 1e3, '.', ms=2, color='0.3'); ax[2].set_xlim(-12, 12)
    ax[2].set_xlabel(f'hours from the box search\'s best dip (period {res.period[k]:.4f} d)'); ax[2].set_ylabel('20-min bins − 1 (ppt)')
    fig.tight_layout(); fig.savefig(f'{OUT}/TIC{tic}_s{SECTOR}_lc.png', dpi=100); plt.close(fig)
    np.savetxt(f'{OUT}/TIC{tic}_s{SECTOR}_lc.csv', np.column_stack([tt, norm, ff]), delimiter=',', fmt=['%.6f', '%.7f', '%.2f'],
               header='time_btjd,flux_normalized,flux_e_per_s (aperture sum minus border-median background; flagged frames dropped)', comments='')
    meta.update({'camera': int(h0['CAMERA']), 'ccd': int(h0['CCD']), 'date_obs': h0.get('DATE-OBS'), 'date_end': h0.get('DATE-END'),
        'frames': int(len(t)), 'frames_used': int(good.sum()), 'cadence_s': round(cad, 1), 'aperture_px': int(keep.sum()),
        'aperture_touches_edge': bool(keep[border].any()), 'scatter_ppm_per_frame': round(rob), 'scatter_ppm_1h_bins': round(brob),
        'deepest_1h_bin_ppt': round((bv[lo] - 1) * 1e3, 3), 'deepest_1h_bin_btjd': round(float(bt[lo]), 4),
        'strongest_sinusoids_cpd_ppt': peaks,
        'box_search': {'period_d': round(float(res.period[k]), 5), 'depth_ppt': round(float(res.depth[k]) * 1e3, 4),
                       'depth_snr': round(float(res.depth[k] / res.depth_err[k]), 2) if res.depth_err[k] > 0 else None,
                       'duration_h': round(float(res.duration[k]) * 24, 2), 'first_mid_btjd': round(float(res.transit_time[k]), 4),
                       'events_with_data': ntr, 'sde': round(sde, 2)}})
    return meta

def main():
    say(f'== {utc()} batch start: sector {SECTOR}, {len(TICS)} TIC IDs')
    rows = tic_rows(TICS)
    allmeta = {}
    for tic in TICS:
        meta = {'tic': tic, 'sector': SECTOR}
        try:
            r = rows.get(tic)
            if not r: raise RuntimeError('not found in the MAST TIC catalog')
            ra, dec = r['ra'], r['dec']
            meta.update({'ra': ra, 'dec': dec, 'tmag': r.get('Tmag'), 'teff': r.get('Teff'), 'contratio': r.get('contratio'), 'hip': r.get('HIP')})
            sec = sectors(ra, dec); meta['sectors_in_tesscut'] = sorted(sec)
            if SECTOR not in sec:
                meta['status'] = f'not in sector {SECTOR} (TESSCut lists {sorted(sec)})'; say(tic, meta['status']); allmeta[tic] = meta; continue
            url = f'https://mast.stsci.edu/tesscut/api/v0.1/astrocut?ra={ra}&dec={dec}&y=15&x=15&sector={SECTOR}'
            d = f'{DL}/{tic}-s{SECTOR}'; os.makedirs(d, exist_ok=True)
            meta['request_url'] = url; meta['fetch_started_utc'] = utc()
            final = fetch(url, f'{d}/cutout.zip'); meta['fetch_done_utc'] = utc(); meta['served_from'] = urllib.parse.urlparse(final).netloc
            with zipfile.ZipFile(f'{d}/cutout.zip') as z:
                names = [n for n in z.namelist() if n.endswith('.fits') and '/' not in n and '..' not in n]
                if len(names) != 1: raise RuntimeError(f'expected one FITS in the zip, found {names}')
                z.extract(names[0], d)
            fp = f'{d}/{names[0]}'
            meta['fits_name'] = names[0]; meta['fits_bytes'] = os.path.getsize(fp)
            meta['fits_sha256'] = hashlib.sha256(open(fp, 'rb').read()).hexdigest()
            analyse(fp, tic, meta); meta['status'] = 'ok'
            shutil.rmtree(d)                        # the cutout can be fetched again from request_url; its hash is kept
            say(tic, 'ok', meta['fetch_done_utc'], meta['frames_used'], 'frames', meta['scatter_ppm_1h_bins'], 'ppm/h',
                'box', meta['box_search']['period_d'], 'd', meta['box_search']['depth_ppt'], 'ppt sde', meta['box_search']['sde'])
        except Exception as e:
            meta['status'] = f'failed: {type(e).__name__}: {str(e)[:200]}'; say(tic, meta['status'])
        allmeta[tic] = meta
        json.dump(meta, open(f'{OUT}/TIC{tic}_s{SECTOR}_meta.json', 'w'), indent=1)
    json.dump(allmeta, open(f'{OUT}/all_s{SECTOR}.json', 'w'), indent=1)
    say(f'== {utc()} batch done: ' + json.dumps({s: sum(1 for m in allmeta.values() if m['status'].startswith(s)) for s in ('ok', 'not in', 'failed')}))

main()
