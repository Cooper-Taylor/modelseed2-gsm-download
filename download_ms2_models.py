#!/usr/bin/env python3
"""Download the ModelSEED v2 genome-scale metabolic models from KBase.

Source: the public "ModelSEED 2 Manuscript Data" KBase organization
        https://narrative.kbase.us/#/orgs/ms2-manuscript
Paper : Faria et al. (2023), "ModelSEED v2: High-throughput genome-scale
        metabolic model reconstruction with enhanced energy biosynthesis
        pathway prediction", bioRxiv 10.1101/2023.10.04.556561

Two public, world-readable workspaces (globalread='r', so NO auth token is
needed) each hold 5,420 ``KBaseFBA.FBAModel`` objects, one per KEGG
Bacteria/Archaea genome:

    155807  ModelSEED 2 Genome-scale Models - Glucose Minimal Media  -> gmm/
    155808  ModelSEED 2 Genome-scale Models - Auxotrophy Media       -> auxotrophy/

Objects are written one-per-file as gzipped JSON, byte-faithful to what the
workspace returns (``data`` payload only; the object ``info`` tuple is recorded
in the manifest). The run is resumable: files already present and non-empty are
skipped unless --force is given.

Usage
-----
    python3 download_ms2_models.py                    # both workspaces
    python3 download_ms2_models.py --ws 155807        # just GMM
    python3 download_ms2_models.py --limit 25         # smoke test
    python3 download_ms2_models.py --workers 12       # more parallelism
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import os
import queue
import random
import sys
import threading
import time
import urllib.error
import urllib.request
from pathlib import Path

WS_URL = "https://kbase.us/services/ws"

WORKSPACES = {
    155807: {
        "slug": "gmm",
        "narrative": "ModelSEED 2 Genome-scale Models - Glucose Minimal Media",
        "media": "Carbon-D-Glucose",
    },
    155808: {
        "slug": "auxotrophy",
        "narrative": "ModelSEED 2 Genome-scale Models - Auxotrophy Media",
        "media": "auxotrophy media",
    },
}

DEFAULT_OUT = Path(
    os.environ.get("MS2_GS_MODELS_DIR", "/scratch/ctaylor/modelseed2_gs_models")
)

_print_lock = threading.Lock()


def log(msg: str) -> None:
    with _print_lock:
        print(f"[{time.strftime('%H:%M:%S')}] {msg}", flush=True)


# --------------------------------------------------------------------------
# Workspace JSON-RPC (unauthenticated -- these workspaces are public)
# --------------------------------------------------------------------------
def rpc(method: str, params, timeout: int = 300, attempts: int = 6):
    """Call Workspace.<method>. Retries with exponential backoff + jitter."""
    body = json.dumps(
        {"version": "1.1", "method": f"Workspace.{method}", "params": params, "id": "1"}
    ).encode()
    last = None
    for i in range(attempts):
        try:
            req = urllib.request.Request(
                WS_URL, data=body, headers={"Content-Type": "application/json"}
            )
            with urllib.request.urlopen(req, timeout=timeout) as fh:
                payload = json.load(fh)
            if "error" in payload:
                raise RuntimeError(payload["error"].get("message", "unknown ws error"))
            return payload["result"][0]
        except Exception as exc:  # noqa: BLE001 -- retry anything transient
            last = exc
            if i == attempts - 1:
                break
            sleep = min(60.0, (2.0 ** i) + random.uniform(0, 1.5))
            time.sleep(sleep)
    raise RuntimeError(f"{method} failed after {attempts} attempts: {last}") from last


def list_models(wsid: int) -> list[list]:
    """Every KBaseFBA.FBAModel object-info tuple in the workspace."""
    out, min_id = [], 0
    while True:
        page = rpc(
            "list_objects",
            [{"ids": [wsid], "type": "KBaseFBA.FBAModel", "limit": 10000,
              "includeMetadata": 0, "minObjectID": min_id}],
        )
        if not page:
            break
        out.extend(page)
        last = max(o[0] for o in page)
        if last <= min_id:
            break
        min_id = last + 1
    out.sort(key=lambda o: o[0])
    return out


# --------------------------------------------------------------------------
# Worker
# --------------------------------------------------------------------------
def target_path(root: Path, slug: str, name: str) -> Path:
    # Object names look like GCF_000007605.1.RAST.GMM.mdl -> keep verbatim.
    return root / slug / "models" / f"{name}.json.gz"


def fetch_batch(wsid: int, infos: list[list], root: Path, slug: str) -> list[dict]:
    """Fetch and persist a batch of objects; return manifest rows."""
    refs = [{"ref": f"{wsid}/{o[0]}/{o[4]}"} for o in infos]
    res = rpc("get_objects2", [{"objects": refs}])
    rows = []
    for obj in res["data"]:
        info = obj["info"]
        objid, name, otype, ver, chsum, size = info[0], info[1], info[2], info[4], info[8], info[9]
        data = obj["data"]
        path = target_path(root, slug, name)
        path.parent.mkdir(parents=True, exist_ok=True)
        tmp = path.with_suffix(path.suffix + ".part")
        with gzip.open(tmp, "wt", encoding="utf-8", compresslevel=6) as fh:
            json.dump(data, fh, separators=(",", ":"), sort_keys=True)
        tmp.replace(path)
        rows.append(
            {
                "wsid": wsid,
                "workspace_slug": slug,
                "objid": objid,
                "version": ver,
                "name": name,
                "type": otype,
                "ws_md5": chsum,
                "ws_bytes": size,
                "gz_bytes": path.stat().st_size,
                "genome_id": name.split(".RAST.")[0] if ".RAST." in name else "",
                "genome_ref": data.get("genome_ref", ""),
                "template_ref": data.get("template_ref", ""),
                "n_reactions": len(data.get("modelreactions", [])),
                "n_compounds": len(data.get("modelcompounds", [])),
                "n_biomasses": len(data.get("biomasses", [])),
                "n_gapfillings": len(data.get("gapfillings", []) or []),
                "path": str(path.relative_to(root)),
            }
        )
    return rows


def download_workspace(
    wsid: int, root: Path, workers: int, batch: int, limit: int | None, force: bool
) -> None:
    meta = WORKSPACES[wsid]
    slug = meta["slug"]
    ws_root = root / slug
    ws_root.mkdir(parents=True, exist_ok=True)

    log(f"ws {wsid} ({slug}): listing FBAModel objects ...")
    infos = list_models(wsid)
    log(f"ws {wsid} ({slug}): {len(infos)} models in workspace")
    if limit:
        infos = infos[:limit]
        log(f"ws {wsid} ({slug}): --limit {limit} -> {len(infos)} to consider")

    if not force:
        todo = [o for o in infos
                if not (p := target_path(root, slug, o[1])).exists() or p.stat().st_size == 0]
        log(f"ws {wsid} ({slug}): {len(infos) - len(todo)} already on disk, {len(todo)} to fetch")
    else:
        todo = infos

    # Save the workspace's media object alongside the models -- the models are
    # gap-filled against it, so any FBA rerun needs it.
    try:
        media = rpc("list_objects", [{"ids": [wsid], "type": "KBaseBiochem.Media",
                                      "limit": 10, "includeMetadata": 0}])
        for m in media:
            mref = f"{wsid}/{m[0]}/{m[4]}"
            mdata = rpc("get_objects2", [{"objects": [{"ref": mref}]}])["data"][0]["data"]
            mpath = ws_root / "media" / f"{m[1]}.json"
            mpath.parent.mkdir(parents=True, exist_ok=True)
            mpath.write_text(json.dumps(mdata, indent=2, sort_keys=True))
            log(f"ws {wsid} ({slug}): saved media {m[1]}")
    except Exception as exc:  # noqa: BLE001
        log(f"ws {wsid} ({slug}): WARNING could not save media: {exc}")

    if not todo:
        log(f"ws {wsid} ({slug}): nothing to do")
    else:
        batches = [todo[i:i + batch] for i in range(0, len(todo), batch)]
        work: queue.Queue = queue.Queue()
        for b in batches:
            work.put(b)

        rows_lock = threading.Lock()
        all_rows: list[dict] = []
        failures: list[tuple[str, str]] = []
        done = [0]
        t0 = time.time()

        def worker(wid: int) -> None:
            while True:
                try:
                    b = work.get_nowait()
                except queue.Empty:
                    return
                try:
                    rows = fetch_batch(wsid, b, root, slug)
                    with rows_lock:
                        all_rows.extend(rows)
                except Exception as exc:  # noqa: BLE001
                    with rows_lock:
                        for o in b:
                            failures.append((o[1], str(exc)[:200]))
                    log(f"  ! batch starting {b[0][1]} FAILED: {str(exc)[:160]}")
                finally:
                    with rows_lock:
                        done[0] += len(b)
                        n = done[0]
                    if n % 200 < batch:
                        el = time.time() - t0
                        rate = n / el if el else 0
                        eta = (len(todo) - n) / rate if rate else 0
                        log(f"  ws {wsid} ({slug}): {n}/{len(todo)} "
                            f"({rate:.1f} mdl/s, ETA {eta/60:.1f} min)")
                    work.task_done()

        threads = [threading.Thread(target=worker, args=(i,), daemon=True)
                   for i in range(workers)]
        for t in threads:
            t.start()
        for t in threads:
            t.join()

        log(f"ws {wsid} ({slug}): fetched {len(all_rows)}, failed {len(failures)} "
            f"in {(time.time() - t0)/60:.1f} min")
        if failures:
            fp = ws_root / "failures.tsv"
            with fp.open("w", newline="") as fh:
                w = csv.writer(fh, delimiter="\t")
                w.writerow(["name", "error"])
                w.writerows(failures)
            log(f"ws {wsid} ({slug}): failures written to {fp}")

    write_manifest(wsid, root, slug, infos)


def write_manifest(wsid: int, root: Path, slug: str, infos: list[list]) -> None:
    """Manifest built by re-reading what is actually on disk (authoritative)."""
    ws_root = root / slug
    cols = ["wsid", "workspace_slug", "objid", "version", "name", "type", "ws_md5",
            "ws_bytes", "gz_bytes", "genome_id", "genome_ref", "template_ref",
            "n_reactions", "n_compounds", "n_biomasses", "n_gapfillings", "path"]
    by_name = {o[1]: o for o in infos}
    rows, missing = [], []
    for name, o in sorted(by_name.items()):
        p = target_path(root, slug, name)
        if not p.exists() or p.stat().st_size == 0:
            missing.append(name)
            continue
        try:
            with gzip.open(p, "rt", encoding="utf-8") as fh:
                d = json.load(fh)
        except Exception as exc:  # noqa: BLE001
            missing.append(f"{name} (unreadable: {exc})")
            continue
        rows.append({
            "wsid": wsid, "workspace_slug": slug, "objid": o[0], "version": o[4],
            "name": name, "type": o[2], "ws_md5": o[8], "ws_bytes": o[9],
            "gz_bytes": p.stat().st_size,
            "genome_id": name.split(".RAST.")[0] if ".RAST." in name else "",
            "genome_ref": d.get("genome_ref", ""), "template_ref": d.get("template_ref", ""),
            "n_reactions": len(d.get("modelreactions", [])),
            "n_compounds": len(d.get("modelcompounds", [])),
            "n_biomasses": len(d.get("biomasses", [])),
            "n_gapfillings": len(d.get("gapfillings", []) or []),
            "path": str(p.relative_to(root)),
        })
    mp = ws_root / "manifest.tsv"
    with mp.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    log(f"ws {wsid} ({slug}): manifest -> {mp} ({len(rows)} rows, {len(missing)} missing)")
    if missing:
        (ws_root / "missing.txt").write_text("\n".join(missing) + "\n")
        log(f"ws {wsid} ({slug}): {len(missing)} missing listed in {ws_root/'missing.txt'}")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--ws", type=int, action="append", choices=sorted(WORKSPACES),
                    help="workspace id (repeatable); default both")
    ap.add_argument("--workers", type=int, default=8)
    ap.add_argument("--batch", type=int, default=4, help="objects per get_objects2 call")
    ap.add_argument("--limit", type=int, default=None, help="cap models per workspace")
    ap.add_argument("--force", action="store_true", help="re-download existing files")
    ap.add_argument("--manifest-only", action="store_true",
                    help="skip fetching; just rebuild manifests from disk")
    args = ap.parse_args()

    root: Path = args.out
    root.mkdir(parents=True, exist_ok=True)
    targets = args.ws or sorted(WORKSPACES)

    # Provenance and layout are documented in README.md at the root of this
    # directory; it is hand-maintained and deliberately not overwritten here.

    for wsid in targets:
        if args.manifest_only:
            write_manifest(wsid, root, WORKSPACES[wsid]["slug"], list_models(wsid))
        else:
            download_workspace(wsid, root, args.workers, args.batch, args.limit, args.force)

    log("done")
    return 0


if __name__ == "__main__":
    sys.exit(main())
