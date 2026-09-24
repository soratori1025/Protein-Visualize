"""
Membrane geometry layer: WHERE the bilayer is, independent of secondary structure.

Why a separate layer: a DSSP helix boundary is not a membrane boundary. A TM helix
often continues well beyond the lipid slab, so snapping a TM segment to the helix
ends puts residues that sit in water "inside the membrane". Boundaries are decided
here, from geometry; DSSP/STRIDE only describe what the membrane part looks like.

Sources, most to least authoritative:
  1. "file planes"   - the structure already carries the membrane as dummy atoms
                       (residue DUM, atoms N and O on the two planes): OPM downloads,
                       PPM 2/3 server output, memembed output. The plane normal and
                       thickness are taken as given. Run PPM on a structure and upload
                       its output to get physics-based placement.
  2. "segment axes"  - normal = dominant direction of the CA(i+4)-CA(i) vectors of
                       the TM segments proposed by the TM block (TM helices/strands
                       point across the bilayer); centre = mean depth of their
                       residues; thickness = TMParams.membrane_thickness.
  3. "slab fit"      - the hydrophobic slab of GeometryTMProvider (used by that
                       provider only; on FucP 3O7Q it is 44 deg off the OPM normal
                       while the axis estimate is 6 deg off, hence rank 3).

Every residue gets a signed depth d (A, 0 = bilayer mid-plane) and a zone:
  CORE |d| <= half - EDGE_WIDTH,  EDGE up to half + EDGE_WIDTH,  OUT beyond.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Sequence

import numpy as np

from app.services.topology.residues import ResidueFrame

EDGE_WIDTH = 4.0            # A either side of the hydrophobic boundary = interface band
ENVELOPE_TOLERANCE = 2.0    # A beyond the boundary still counted "in the membrane"
MIN_DUM_ATOMS = 3


@dataclass
class MembraneFrame:
    normal: np.ndarray          # unit vector
    mid_point: np.ndarray       # a point on the mid-plane (absolute coordinates)
    half_thickness: float       # A
    source: str
    depth: np.ndarray           # signed distance of every CA from the mid-plane (A)

    def zone(self, i: int) -> Optional[str]:
        d = self.depth[i]
        if not np.isfinite(d):
            return None
        a = abs(d)
        if a <= self.half_thickness - EDGE_WIDTH:
            return "CORE"
        if a <= self.half_thickness + EDGE_WIDTH:
            return "EDGE"
        return "OUT"

    def zones(self) -> list[Optional[str]]:
        return [self.zone(i) for i in range(len(self.depth))]

    def in_envelope(self, i: int, tolerance: float = ENVELOPE_TOLERANCE) -> bool:
        d = self.depth[i]
        return bool(np.isfinite(d) and abs(d) <= self.half_thickness + tolerance)

    def span(self, a: int, b: int) -> float:
        """Extent of residues a..b along the normal (A)."""
        seg = self.depth[a:b + 1]
        seg = seg[np.isfinite(seg)]
        return float(seg.max() - seg.min()) if seg.size else 0.0

    def recentered(self, element_spans: Sequence[tuple[int, int]]) -> "MembraneFrame":
        """Re-place the mid-plane at the median mid-depth of membrane-spanning SS
        elements (every crossing is centred on the bilayer). Only for estimated
        placements - planes read from the file are kept as they are."""
        if self.source.startswith("file") or not element_spans:
            return self
        mids = [float(np.nanmean(self.depth[a:b + 1])) for a, b in element_spans]
        shift = float(np.median(mids))
        if not np.isfinite(shift):
            return self
        return MembraneFrame(self.normal, self.mid_point + shift * self.normal,
                             self.half_thickness, self.source, self.depth - shift)

    def as_dict(self) -> dict:
        return {"source": self.source,
                "normal": [round(float(x), 4) for x in self.normal],
                "half_thickness": round(float(self.half_thickness), 2)}


def _with_depth(frame: ResidueFrame, normal, mid_point, half, source) -> MembraneFrame:
    normal = np.asarray(normal, dtype=float)
    normal = normal / np.linalg.norm(normal)
    depth = (frame.coords - np.asarray(mid_point, dtype=float)) @ normal
    return MembraneFrame(normal, np.asarray(mid_point, dtype=float), float(half), source, depth)


# ----------------------------------------------------------------------------- 1
def read_dummy_planes(file_path: Path):
    """(normal, mid_point, half_thickness) from DUM atoms, or None.
    PDB (HETATM ... DUM) and mmCIF (_atom_site with comp_id DUM) are both read by
    column-free parsing of the coordinates Bio.PDB would see."""
    planes: dict[str, list] = {"N": [], "O": []}
    path = Path(file_path)
    if path.suffix.lower() in (".cif", ".mmcif"):
        try:
            from Bio.PDB.MMCIF2Dict import MMCIF2Dict
            d = MMCIF2Dict(str(path))
            comp = d.get("_atom_site.label_comp_id", [])
            atom = d.get("_atom_site.label_atom_id", [])
            xs, ys, zs = (d.get(f"_atom_site.Cartn_{c}", []) for c in "xyz")
            for c, a, x, y, z in zip(comp, atom, xs, ys, zs):
                if c == "DUM" and a in planes:
                    planes[a].append([float(x), float(y), float(z)])
        except (OSError, ValueError, KeyError):
            return None
    else:
        try:
            with open(path, errors="replace") as fh:
                for line in fh:
                    if line.startswith(("HETATM", "ATOM")) and line[17:20] == "DUM":
                        name = line[12:16].strip()
                        if name in planes:
                            planes[name].append([float(line[30:38]), float(line[38:46]),
                                                 float(line[46:54])])
        except (OSError, ValueError):
            return None
    if min(len(planes["N"]), len(planes["O"])) < MIN_DUM_ATOMS:
        return None
    n_pts, o_pts = np.array(planes["N"]), np.array(planes["O"])
    c_n, c_o = n_pts.mean(axis=0), o_pts.mean(axis=0)
    # plane normal = smallest singular direction of both point clouds together
    _, _, vt = np.linalg.svd(np.vstack([n_pts - c_n, o_pts - c_o]))
    normal = vt[-1]
    if normal @ (c_o - c_n) < 0:
        normal = -normal
    thickness = float((c_o - c_n) @ normal)
    if not np.isfinite(thickness) or thickness < 10.0:
        return None
    mid = (c_n + c_o) / 2.0
    return normal, mid, thickness / 2.0


def membrane_from_file(frame: ResidueFrame) -> Optional[MembraneFrame]:
    if frame.source_path is None or len(frame) == 0:
        return None
    planes = read_dummy_planes(frame.source_path)
    if planes is None:
        return None
    normal, mid, half = planes
    return _with_depth(frame, normal, mid, half, "file planes (DUM atoms: OPM/PPM/memembed)")


# ----------------------------------------------------------------------------- 2
def axis_normal(coords: np.ndarray, spans: Sequence[tuple[int, int]], breaks, k: int = 4):
    """Top eigenvector of the scatter of unit CA(i+k)-CA(i) vectors inside `spans`
    (sign-free, so up- and down-going elements reinforce each other)."""
    dirs = []
    for s, e in spans:
        for i in range(s, e - k + 1):
            if breaks is not None and np.any(breaks[i + 1:i + k + 1]):
                continue
            v = coords[i + k] - coords[i]
            n = np.linalg.norm(v)
            if n > 0 and np.isfinite(n):
                dirs.append(v / n)
    if len(dirs) < 3:
        return None
    d = np.array(dirs)
    _, vec = np.linalg.eigh(d.T @ d)
    return vec[:, -1]


def membrane_from_segments(frame: ResidueFrame, segments: Sequence[tuple[int, int]],
                           half_thickness: float) -> Optional[MembraneFrame]:
    if len(frame) == 0 or not np.isfinite(frame.coords).all() or not segments:
        return None
    breaks = frame.chain_breaks()
    normal = axis_normal(frame.coords, segments, breaks)
    if normal is None:
        normal = axis_normal(frame.coords, [(0, len(frame) - 1)], breaks)
    if normal is None:
        return None
    normal = normal / np.linalg.norm(normal)
    proj = frame.coords @ normal
    idx = [k for s, e in segments for k in range(s, e + 1)]
    # each TM element is centred on the bilayer -> average of segment mid-depths
    mids = [float(np.mean(proj[s:e + 1])) for s, e in segments]
    center = float(np.median(mids)) if mids else float(np.mean(proj[idx]))
    mid_point = normal * center
    return _with_depth(frame, normal, mid_point, half_thickness, "TM segment axes")


def build_membrane(frame: ResidueFrame, segments: Sequence[tuple[int, int]],
                   half_thickness: float) -> Optional[MembraneFrame]:
    """Most authoritative membrane placement available for this structure."""
    return membrane_from_file(frame) or membrane_from_segments(frame, segments, half_thickness)
