# point2all_prop_dumper — 24c benchmarks + production sampling plan

Companion to `point2all_prop_dumper_remote_handoff_claude.md`. Records the pdebug single-config
benchmarks, the per-ensemble sampling strides, and the `PROP_TPT_SECONDS` choice used by the
production submits (`submit_prop_dumper_b<beta>_m<mass>_claude.sh`) + unified chain driver
(`run_prop_dumper_chain_claude.sh`), all under `/usr/workspace/lsd/matsumoto5/su4_32c`.

All runs: MI300A, `-N 8 -n 32 -g 1`, `--grid 24.24.24.48 --mpi 2.2.2.4 --split 1 1 1 4`
(nrhs=8, nbatch=1, all 8 corner sources solved concurrently), `W=3.0 N=40`. Output per config is one
0.68 GB SciDAC/LIME with 8 single-precision q00 records.

## Benchmark results (config lat.1000, pdebug)

Every run PASSED the handoff sanity checks: 8 `ildg-binary-data` records, `srcIndex` 0-7,
`<store>q00</store>` x8, `<prec>single</prec>`, 0.68 GB, q00 norms finite ~4-7e-4, no errors/asserts.

| mass | beta   | job id         | date       | per-config wall |
|------|--------|----------------|------------|-----------------|
| 0.4  | 11.045 | f3WP32KZ8Tfd   | 2026-09-04 | 61 s            |
| 0.3  | 11.035 | f3WXcf6Wd6Mu   | 2026-09-05 | 61 s            |
| 0.2  | 10.990 | f3WXceyTSWaj   | 2026-09-05 | 68 s            |
| 0.1  | 10.865 | f3WXcerxsf19   | 2026-09-05 | 78 s            |

Trend is monotonic: lighter quark -> more CG iterations -> slower; m0.1 (lightest of the heavy set)
is the slowest at 78 s. (m0.05, m0.01 deferred: they need deflation.)

## Sampling plan (production)

On-disk HMC checkpoint interval is **20** for every ensemble (`lat.1000, 1020, ...`), so any
measurement `CONF_STRIDE` must be a multiple of 20. `CONF_MIN=1000`; all listed indices exist (0
missing). Numbers can be increased later (finer stride) if more statistics are needed.

| mass | beta   | CONF_STRIDE | configs | index range   | on-disk (stride 20) |
|------|--------|-------------|---------|---------------|---------------------|
| 0.1  | 10.865 | 40  (1/2)   | 92      | 1000 .. 4640  | 183 (max 4640)      |
| 0.2  | 10.990 | 40  (1/2)   | 93      | 1000 .. 4680  | 185 (max 4680)      |
| 0.3  | 11.035 | 80  (1/4)   | 88      | 1000 .. 7960  | 350 (max 7980)      |
| 0.4  | 11.045 | 80  (1/4)   | 113     | 1000 .. 9960  | 450 (max 9980)      |
|      |        |             | **386** |               |                     |

## PROP_TPT_SECONDS = 150 (uniform), and why

`PROP_TPT_SECONDS` is the *estimated per-config wall* that feeds BOTH graceful-blocker layers in the
submit + binary:
- **shell-side**: before launching config i, if `flux job timeleft < PROP_TPT_SECONDS + BLOCKER_OVERHEAD`
  (`BLOCKER_OVERHEAD=300`), the loop stops cleanly (resubmit/chain continues).
- **in-binary** (`PROP_DEADLINE_EPOCH` exported from `flux job timeleft`): the dumper refuses to START
  the solve+write if `now + 1.2*PROP_TPT_SECONDS > deadline`, exiting before opening any file.

So the effective reservation before the wall is `max(PROP_TPT_SECONDS+300 shell, 1.2*PROP_TPT_SECONDS
in-binary)`. With 150: shell reserves 450 s, in-binary 180 s.

**Why 150 and not tighter (e.g. 80-100):**
- Measured walls span 61-78 s; 150 ~= 2x the slowest (78 s), leaving headroom for config-to-config
  variation (CG iteration count varies with the gauge field; occasional slow-converging configs,
  more likely at lighter mass), first-config lustre load/first-touch, and node-to-node jitter.
- The cost of over-reserving is tiny: at end-of-job the blocker idles at most ~450 s (~7.5 min) out of
  a 480 m job (< 1.6 %). Under-reserving risks launching a config that gets wall-killed mid-solve --
  wasted GPU-hours (the atomic write still prevents a corrupt .prop.lime, but the compute is lost).
- **Uniform vs per-mass**: the spread (61 vs 78 s) is small enough that one value for all four masses
  removes a per-script knob with no meaningful throughput cost. Override per job with
  `PROP_TPT_SECONDS=<n>` if a specific ensemble proves slower in production.

## Wall-time budget (measured rate x config count) — each fits ONE -t 480m job

| mass | configs x wall        | ~total  |
|------|-----------------------|---------|
| 0.1  | 92  x 78 s            | ~120 m  |
| 0.2  | 93  x 68 s            | ~105 m  |
| 0.3  | 88  x 61 s            | ~90 m   |
| 0.4  | 113 x 61 s            | ~115 m  |

So no chaining is strictly required (CHAINLEN=1); the chain driver still gives free resume if a job is
cut short. Skip a finished mass by commenting out its line in `run_prop_dumper_chain_claude.sh`
(`ensembles` array), or run one with `ONLY=<substr>`.
