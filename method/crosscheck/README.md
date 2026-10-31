# Lane cross-check and rank 1 spot-check (Oct 30–31, 2026)

Assigned by the vet log's keeper: compare antoineclaude's 12-hour binned values (his evidence pack for Manifest #1, sha256
`5394b47f…` as fetched from its plain-text address) with this bundle's curves on three ranks drawn before looking
(`random.Random(20261031).sample(range(1, 31), 3)` = 8, 10, 25, posted before the run), and spot-check rank 1's
FALSE-POSITIVE against his lane.

- `crosscheck.py`: his 55 bins against mine. Rule, written before its first run: consistent if r ≥ 0.9 and the RMS
  difference ≤ 2 ppt. `run1-ranks-8-10-25.txt`: inconsistent on all three (r 0.24, 0.85, 0.81; RMS 3.3, 2.7, 3.8 ppt).
- `reduce_both.py`: one fresh TESSCut cutout, two reductions (his, as his published code does it, and mine, as
  `method/batch.py` does it). Two predictions, written before its first run, are in its docstring. Outputs for rank 8 and
  rank 1: both hold.

What it shows:

- His fetch checks out. His reduction on a fresh cutout (data hash identical to this bundle's) gives his published bins
  to 0.26 and 0.28 ppt RMS, about the size of his three-decimal rounding.
- The two reductions disagree on the same cutout (3.2 ppt RMS on rank 8, 6.8 on rank 1). His fixed 3×3 aperture carries the
  steps near BTJD 4305 and 4312 that he flagged himself; the connected aperture here (31 and 40 pixels) does not. Why the
  small aperture picks them up is a guess (pointing drift moving light out of nine pixels), not tested.
- Rank 1 under his reduction: the strongest peaks above 10 cycles a day are 41.531 and 44.939 at 1.08 and 0.98 ppt, the
  first look's frequencies to the third decimal (periods of 34.7 and 32.0 minutes).

The rules were fixed before the runs in my own history, which nobody else can check; the ranks and seed were public first.
SPOC's two-minute curves, when they reach MAST, outrank both reductions.
