"""
Membrane-frame benchmark: how accurately is the bilayer placed?

    python bench/membrane_frame_benchmark.py --results bench/results

For every structure of the evaluation cohort that is also covered by the
Orientation of Proteins in Membranes (OPM) database, the OPM/PPM-oriented coordinate
file is downloaded from OPM.  That file carries the reference bilayer as a pair of
`DUM` atom planes, from which the reference membrane frame is computed exactly:

  * the normal is the smallest singular vector of the two DUM point clouds;
  * the hydrophobic thickness is the distance between the two planes along that
    normal;
  * the mid-plane is the mean of the two plane centres.

The frame fitted by ProteinLab is then compared with that reference:

  * **normal error** - the angle between the fitted normal and the reference normal,
    as an axis (0-90 degrees; the sign of a membrane normal is not observable);
  * **mid-plane offset** - the distance between the two mid-planes measured along the
    fitted normal, in angstrom;
  * **thickness error** - the difference between the fitted hydrophobic thickness and
    the reference bilayer thickness, in angstrom.

The 3D slab geometry provider is run on the same files as a control: the question is
not only whether the frame is accurate, but whether the more elaborate energy model
places it better.

Reference: Lomize MA, Blagovestova PV, Pogozheva ID. Membrane orientation and
curvature in OPM/PPM.  OPM distributes the oriented files at
https://opm-assets.storage.googleapis.com/pdb/{pdb_id}.pdb
"""
from __future__ import annotations

import argparse
import json
import math
import sys
import urllib.request
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "backend"))

import numpy as np  # noqa: E402

from app.schemas.topology import TMParams  # noqa: E402
from app.services.topology import membrane_energy as ME  # noqa: E402
from app.services.topology.providers.tm.energy import EnergyTMProvider  # noqa: E402
from app.services.topology.providers.tm.uniprot import fetch_uniprot_entry  # noqa: E402
from app.services.topology.residues import load_residue_frame  # noqa: E402
from app.services.topology.topology_predictor import GeometryTMProvider  # noqa: E402

OPM_URL = "https://opm-assets.storage.googleapis.com/pdb/{pdb}.pdb"
UA = {"User-Agent": "ProteinLab-benchmark/1.0 (academic validation)"}
MIN_DUM_ATOMS = 3


def download(url: str, dest: Path) -> bool:
    if dest.exists() and dest.stat().st_size > 1000:
        return True
    try:
        req = urllib.request.Request(url, headers=UA)
        with urllib.request.urlopen(req, timeout=120) as resp:
            data = resp.read()
    except Exception:
        return False
    if len(data) < 1000:
        return False
    dest.parent.mkdir(parents=True, exist_ok=True)
    dest.write_bytes(data)
    return True


def read_dum_planes(path: Path):
    """The reference membrane frame from the OPM DUM bilayer atoms.

    Returns (normal, half_thickness, mid_plane) or None when the file has no usable
    bilayer.  The normal is oriented from the lower to the upper plane.
    """
    lower, upper = [], []
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.startswith(("ATOM  ", "HETATM")):
                continue
            if line[17:20].strip() != "DUM":       # residue name, columns 18-20
                continue
            atom = line[12:16].strip()              # atom name, columns 13-16
            try:
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            # OPM writes one dummy atom per bilayer atom: N = lower plane, O = upper
            (lower if atom.upper().startswith("N") else upper).append(xyz)
    if len(lower) < MIN_DUM_ATOMS or len(upper) < MIN_DUM_ATOMS:
        return None
    lo = np.asarray(lower, float)
    up = np.asarray(upper, float)
    c_lo, c_up = lo.mean(axis=0), up.mean(axis=0)
    cloud = np.vstack([lo - c_lo, up - c_up])
    _, _, vt = np.linalg.svd(cloud)
    normal = vt[-1]
    if normal @ (c_up - c_lo) < 0:
        normal = -normal
    half = float(abs((c_up - c_lo) @ normal) / 2.0)
    mid = (c_lo + c_up) / 2.0
    return normal, half, mid


def read_ca_pdb(path: Path):
    """C-alpha coordinates of a PDB file, keyed by (chain, resseq, icode)."""
    out = {}
    with open(path, errors="replace") as fh:
        for line in fh:
            if not line.startswith("ATOM  "):
                continue
            if line[12:16].strip() != "CA":
                continue
            try:
                out[(line[21].strip(), int(line[22:26]), line[26].strip())] = (
                    float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
    return out


def kabsch(mobile, target):
    """Rotation and translation that bring `mobile` onto `target` (both N x 3)."""
    cm = mobile.mean(axis=0)
    ct = target.mean(axis=0)
    h = (mobile - cm).T @ (target - ct)
    u, _, vt = np.linalg.svd(h)
    d = np.sign(np.linalg.det(vt.T @ u.T))
    r = vt.T @ np.diag([1.0, 1.0, d]) @ u.T
    t = ct - r @ cm
    return r, t


def map_reference_to_input(opm_path, input_path, n_opm, mid_opm):
    """Bring the OPM reference frame into the coordinate frame of the input file.

    OPM distributes a *reoriented* copy of the deposited structure.  The reference
    frame therefore has to be rotated back into the frame of the file the user
    actually provides; the rotation is taken from a least-squares superposition of
    the C-alpha atoms of the two files.
    """
    opm_ca = read_ca_pdb(opm_path)
    frame = load_residue_frame(input_path, None)
    in_ca = {}
    for i, res in enumerate(frame.residues):
        in_ca[(res.chain, res.resseq, res.icode)] = tuple(float(v) for v in frame.coords[i])
    keys = [k for k in in_ca if k in opm_ca]
    if len(keys) < 3:
        return None
    mobile = np.array([opm_ca[k] for k in keys], float)
    target = np.array([in_ca[k] for k in keys], float)
    r, t = kabsch(mobile, target)
    rmsd = float(np.sqrt(((mobile @ r.T + t - target) ** 2).sum(axis=1).mean()))
    return r @ n_opm, r @ mid_opm + t, rmsd, len(keys)


def axis_angle(a: np.ndarray, b: np.ndarray) -> float:
    """Angle between two axes in [0, 90] degrees (the sign of a normal is free)."""
    a = a / np.linalg.norm(a)
    b = b / np.linalg.norm(b)
    ang = math.degrees(math.acos(max(-1.0, min(1.0, float(a @ b)))))
    return min(ang, 180.0 - ang)


def analyse(job):
    acc, pdb_id, chain, opm_path, input_path = job
    ref = read_dum_planes(Path(opm_path))
    if ref is None:
        return None
    n_opm, half_ref, mid_opm = ref
    mapped = map_reference_to_input(Path(opm_path), Path(input_path), n_opm, mid_opm)
    if mapped is None:
        return None
    n_ref, mid_ref, fit_rmsd, n_atoms = mapped
    out = {"accession": acc, "pdb_id": pdb_id, "chain": chain,
           "reference_thickness_A": round(2 * half_ref, 2),
           "superposition_rmsd_A": round(fit_rmsd, 3),
           "n_superposed_atoms": n_atoms}

    # --- the energy model -----------------------------------------------------
    try:
        frame = load_residue_frame(Path(input_path), chain or None)
        fit = ME.fit_structure(Path(input_path), ca_coords=np.asarray(frame.coords, float),
                               ca_names=frame.names)
        pl = fit.placement
        n_fit = np.asarray(pl.normal, float)
        d = float((mid_ref - pl.center) @ n_fit)
        out["3d_energy"] = {
            "normal_error_deg": round(axis_angle(n_fit, n_ref), 2),
            "midplane_offset_A": round(abs(d), 2),
            "thickness_A": round(2 * pl.half_thickness, 2),
            "thickness_error_A": round(2 * pl.half_thickness - 2 * half_ref, 2),
            "dG_transfer": round(float(pl.energy), 2),
            "curved": bool(pl.curved),
            "curvature_radius_A": (round(float(pl.radius), 1)
                                   if pl.curved and pl.radius else None),
        }
    except Exception as exc:
        out["3d_energy"] = {"error": f"{type(exc).__name__}: {exc}"}

    # --- the slab geometry control -------------------------------------------
    try:
        pred = GeometryTMProvider().predict_tm(Path(input_path), params=TMParams(), frame=frame)
        n_slab = np.asarray(pred.membrane_normal, float)
        if n_slab.size == 3:
            out["3d_slab_geom"] = {
                "normal_error_deg": round(axis_angle(n_slab, n_ref), 2),
            }
    except Exception as exc:
        out["3d_slab_geom"] = {"error": f"{type(exc).__name__}: {exc}"}
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench/results")
    ap.add_argument("--cohort", default=None,
                    help="a frame_cohort.json built by bench/frame_cohort.py; the "
                         "default is the topology cohort in --results")
    ap.add_argument("--out", default=None)
    ap.add_argument("--limit", type=int, default=0, help="0 = all")
    args = ap.parse_args()
    results = Path(args.results)
    out_dir = Path(args.out) if args.out else results
    structs = results / "structures"
    structs.mkdir(parents=True, exist_ok=True)

    if args.cohort:
        # an independent cohort, built for the frame benchmark alone
        cohort = json.loads(Path(args.cohort).read_text())
        jobs = []
        for r in cohort:
            opm = Path(r["opm_file"])
            original = Path(r["structure_file"])
            if not opm.exists() or not original.exists():
                continue
            jobs.append((r["accession"], r["pdb_id"], r.get("chain") or "", str(opm),
                         str(original)))
        print(f"[frame-bench] {len(jobs)} structures from the independent frame cohort")
    else:
        records = json.loads((results / "benchmark.json").read_text())
        if args.limit:
            records = records[: args.limit]
        opm_dir = results / "opm"
        opm_dir.mkdir(parents=True, exist_ok=True)
        jobs = []
        n_missing = 0
        for r in records:
            dest = opm_dir / f"{r['pdb_id']}.pdb"
            if not download(OPM_URL.format(pdb=r["pdb_id"].lower()), dest):
                n_missing += 1
                continue
            original = structs / f"{r['pdb_id']}.cif"
            if not download(f"https://files.rcsb.org/download/{r['pdb_id'].upper()}.cif",
                            original):
                n_missing += 1
                continue
            jobs.append((r["accession"], r["pdb_id"], r["chain"], str(dest), str(original)))
        print(f"[frame-bench] {len(jobs)} structures with an OPM reference frame "
              f"({n_missing} of {len(records)} without one)")

    rows = []
    for job in jobs:
        row = analyse(job)
        if row is None:
            continue
        rows.append(row)
        if len(rows) % 20 == 0:
            print(f"[frame-bench] {len(rows)} analysed")
    (out_dir / "membrane_frame_benchmark.json").write_text(json.dumps(rows, indent=2))

    def stats(key, field):
        vals = [r[key][field] for r in rows if key in r and field in r[key]]
        if not vals:
            return None
        return (round(float(np.median(vals)), 2), round(float(np.mean(vals)), 2),
                round(float(np.percentile(vals, 90)), 2), round(float(np.max(vals)), 2))

    print(f"\n{'metric':38s} {'n':>4s} {'median':>8s} {'mean':>8s} {'p90':>8s} {'max':>8s}")
    for key, label in (("3d_energy", "3D Energy"), ("3d_slab_geom", "3D slab geometry")):
        for field, name in (("normal_error_deg", "normal error (deg)"),
                            ("midplane_offset_A", "mid-plane offset (A)"),
                            ("thickness_A", "fitted hydrophobic thickness (A)"),
                            ("thickness_error_A", "thickness error (A)")):
            s = stats(key, field)
            if s is None:
                continue
            n = len([r for r in rows if key in r and field in r[key]])
            print(f"{label + ' - ' + name:38s} {n:4d} {s[0]:8.2f} {s[1]:8.2f} {s[2]:8.2f} {s[3]:8.2f}")
    ref = [r["reference_thickness_A"] for r in rows]
    if ref:
        print(f"{'OPM reference bilayer thickness (A)':38s} {len(ref):4d} "
              f"{np.median(ref):8.2f} {np.mean(ref):8.2f} "
              f"{np.percentile(ref, 90):8.2f} {np.max(ref):8.2f}")
    curved = [r for r in rows if r.get("3d_energy", {}).get("curved")]
    print(f"\ncurved placements: {len(curved)}")
    for r in curved:
        print(f"   {r['pdb_id']} R={r['3d_energy']['curvature_radius_A']} A  "
              f"normal error {r['3d_energy']['normal_error_deg']} deg")
    rmsd = [r["superposition_rmsd_A"] for r in rows if "superposition_rmsd_A" in r]
    if rmsd:
        print(f"superposition RMSD (A): median {np.median(rmsd):.2f}, max {max(rmsd):.2f} "
              f"(OPM file vs deposited file)")


if __name__ == "__main__":
    main()
