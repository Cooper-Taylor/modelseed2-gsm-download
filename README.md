# ModelSEED v2 genome-scale model downloader

Downloads the 5,420-model genome-scale metabolic model sets from the
ModelSEED v2 paper's public KBase deposit.

## Provenance

Faria JP, Liu F, Edirisinghe JN, Gupta N, Seaver SMD, Freiburger AP, Zhang Q,
Weisenhorn P, Conrad N, Zarecki R, Song H-S, DeJongh M, Best AA, Cottingham RW,
Arkin AP, Henry CS. *ModelSEED v2: High-throughput genome-scale metabolic model
reconstruction with enhanced energy biosynthesis pathway prediction.*
bioRxiv 2023.10.04.556561. <https://doi.org/10.1101/2023.10.04.556561>

The paper's Data Availability statement points to the public **ModelSEED 2
Manuscript Data** KBase organization:
<https://narrative.kbase.us/#/orgs/ms2-manuscript>. The two genome-scale
model workspaces used here were identified by browsing that org (their IDs
are not named in the paper itself):

| workspace id | narrative | objects |
|---|---|---|
| [155807](https://narrative.kbase.us/narrative/155807) | ModelSEED 2 Genome-scale Models - Glucose Minimal Media | 5,420 `KBaseFBA.FBAModel` |
| [155808](https://narrative.kbase.us/narrative/155808) | ModelSEED 2 Genome-scale Models - Auxotrophy Media | 5,420 `KBaseFBA.FBAModel` |

Both workspaces are world-readable (`globalread = 'r'`), so **no KBase auth
token is required** to run this script.

## What it does

Talks directly to the public KBase Workspace JSON-RPC API
(`https://kbase.us/services/ws`) — no KBase SDK, no other dependency beyond
the Python standard library:

1. `Workspace.list_objects` enumerates every `KBaseFBA.FBAModel` object in a
   workspace (paginated), plus the `KBaseBiochem.Media` object the models
   were gap-filled against.
2. `Workspace.get_objects2` fetches objects in small batches across a pool
   of worker threads.
3. Each object's JSON payload is gzip-compressed and written to
   `<slug>/models/<object name>.json.gz`, via an atomic `.part` → rename so
   an interrupted run can resume safely (already-downloaded files are
   skipped unless `--force` is given).
4. A `manifest.tsv` is rebuilt by re-reading every file back off disk, so it
   reflects what actually landed rather than what the run thinks it wrote.

## Usage

```bash
python3 download_ms2_models.py                 # both workspaces
python3 download_ms2_models.py --ws 155807      # just the glucose-minimal set
python3 download_ms2_models.py --limit 25       # smoke test
python3 download_ms2_models.py --workers 12     # more parallelism
python3 download_ms2_models.py --manifest-only  # rebuild manifest.tsv from what's already on disk
```

Output defaults to `/scratch/ctaylor/modelseed2_gs_models` (override with
`--out` or the `MS2_GS_MODELS_DIR` env var). Expect roughly 2.1 GB total and
~15 minutes per workspace at the default worker count.

## Known gap

The downloaded JSON is written via `json.dump(..., sort_keys=True)` — a
canonical re-serialization, not the original bytes KBase returned. That
means the `ws_md5` checksum recorded in the manifest (from the object's
`info` tuple) can never be used to verify the download against KBase's own
checksum; it doesn't match the locally stored bytes by construction. Gzip's
own CRC and the JSON parse during manifest rebuild will still catch
truncation/corruption, just not a byte-for-byte integrity check against the
source.
