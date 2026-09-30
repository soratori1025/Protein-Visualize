"""
Implicit-membrane energy model for placing a lipid bilayer on a protein structure
(the "3D Energy" TM provider). Self-contained: numpy only.

Model (in the spirit of PPM [Lomize et al.] and TmDet [Tusnady et al.]):

  * Every protein residue i contributes with its solvent-accessible surface, as a
    relative accessibility RSA_i in [0, 1] (Shrake-Rupley SASA of the whole model /
    theoretical maximum). Buried residues - the hydrophobic core of a soluble domain,
    subunit interfaces - never see lipid, so they carry no weight.
  * Depth dependence: a bilayer of hydrophobic half-thickness D is described by two
    smooth occupancy functions of the depth z along the membrane normal,
        core(z)      = s((D - |z|) / L)                  hydrocarbon core
        interface(z) = s((D + W - |z|) / L) - core(z)    head-group / interface region
    (s = logistic function, L = edge softness, W = interface width).
  * Transfer free energies (kcal/mol) per residue type. Core: the GES scale (Engelman,
    Steitz & Goldman), derived for residues of helices in the membrane interior, i.e.
    with the backbone already hydrogen-bonded. Interface: the Wimley-White water -> POPC
    interface scale, which makes Trp/Tyr prefer the membrane boundary (aromatic belt)
    and Lys/Arg prefer the interface over the core. (The Wimley-White octanol scale for
    the core was tested too: it charges every exposed backbone unit ~1.15 kcal/mol,
    which shrinks the fitted thickness to ~24 A; the GES core gives 28-33 A for GPCRs.)

        E(n, c, D) = sum_i w_i * [ dG_core(i) * core(z_i) + dG_if(i) * interface(z_i) ]
        z_i = (x_i - centroid) . n - c,   x_i = side-chain centroid of residue i

  * Weight w_i = RSA_i * pore(i). pore() removes residues that line a water-filled pore
    or a deep cavity: they are solvent-accessible but never touch lipid. It is measured
    as the fraction of straight rays from the residue that escape the protein
    (`escape_fraction`); pore-lining residues only see along the pore axis. Without it
    the charged residues lining the lumen of a porin make the true placement unfavourable.
  * E is minimised over the normal n (hemisphere), the centre c and the half-thickness
    D in [MIN_HALF, MAX_HALF] (hydrophobic thickness 20-40 A by default). The search is
    a global grid search where, for a fixed n, E(c) for all D is a 1-D convolution of the
    residue "energy histogram" along n with a D-dependent kernel (done with FFTs), then
    a local refinement around the best normals.
  * dG_transfer = E_min (negative = favourable). A structure whose best placement is not
    below MAX_MEMBRANE_DG is reported as non-membrane.

`fit_structure` is the single entry point used by the provider (parsing, SASA, pore
correction, search; C-alpha-only fallback with half-sphere exposure).

Checked on 54 PDB structures (2026-09-30; small set, re-validate on OPM before publishing):
  * 7 membrane proteins with a reference membrane frame (PDBTM/OPM), randomly rotated:
    normal error median 8.7 deg (max 13), centre error 2.3 A, hydrophobic thickness
    27.7 A on average (OPM 34 A for 4ib4, fitted 31.7 A); all 33 TM segments of the
    analysed chains found (3D slab: 21/33, it misses the beta-barrel entirely).
  * dG_transfer: 22 membrane proteins -175..-28 kcal/mol, 32 soluble proteins -3.5..+6.5.
  * Known limits: beta-barrels settle at the 20 A lower bound of the thickness; the
    ATP-synthase c-ring (2x2v) at the 40 A upper bound; residues of a plug domain inside
    a beta-barrel can be labelled TM (they sit at membrane depth).

References (verify the tabulated values against the original tables before publishing):
  Engelman, Steitz & Goldman (1986) Annu Rev Biophys Biophys Chem 15:321 (GES scale);
  Wimley & White (1996) Nat Struct Biol 3:842 (interface scale, whole residue);
  Wimley, Creamer & White (1996) Biochemistry 35:5109 (octanol scale, not used by default);
  Tien et al. (2013) PLoS One 8:e80635 (maximum ASA, theoretical);
  Shrake & Rupley (1973) J Mol Biol 79:351; Hamelryck (2005) Proteins 59:38 (HSE).
"""
from __future__ import annotations

import gzip
import math
import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np

# -----------------------------------------------------------------------------
# Residue tables
# -----------------------------------------------------------------------------
STANDARD = ("ALA", "ARG", "ASN", "ASP", "CYS", "GLN", "GLU", "GLY", "HIS", "ILE",
            "LEU", "LYS", "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "VAL")
ALIASES = {"MSE": "MET", "SEC": "CYS", "CYX": "CYS", "CYM": "CYS", "HSD": "HIS",
           "HSE": "HIS", "HSP": "HIS", "HID": "HIS", "HIE": "HIS", "HIP": "HIS",
           "ASH": "ASP", "GLH": "GLU", "LYN": "LYS", "SEP": "SER", "TPO": "THR",
           "PTR": "TYR", "MLY": "LYS", "KCX": "LYS", "HYP": "PRO", "CSO": "CYS",
           "CME": "CYS", "PYL": "LYS"}

# Wimley-White whole-residue transfer free energies, kcal/mol (Asp/Glu/Lys/Arg charged,
# His neutral). Positive = unfavourable.
WW_OCTANOL = {"ALA": 0.50, "ARG": 1.81, "ASN": 0.85, "ASP": 3.64, "CYS": -0.02,
              "GLN": 0.77, "GLU": 3.63, "GLY": 1.15, "HIS": 0.11, "ILE": -1.12,
              "LEU": -1.25, "LYS": 2.80, "MET": -0.67, "PHE": -1.71, "PRO": 0.14,
              "SER": 0.46, "THR": 0.25, "TRP": -2.09, "TYR": -0.71, "VAL": -0.46}
WW_INTERFACE = {"ALA": 0.17, "ARG": 0.81, "ASN": 0.42, "ASP": 1.23, "CYS": -0.24,
                "GLN": 0.58, "GLU": 2.02, "GLY": 0.01, "HIS": 0.17, "ILE": -0.31,
                "LEU": -0.56, "LYS": 0.99, "MET": -0.23, "PHE": -1.13, "PRO": 0.45,
                "SER": 0.13, "THR": 0.14, "TRP": -1.85, "TYR": -0.94, "VAL": 0.07}
# Maximum accessible surface area, A^2 (Tien et al. 2013, theoretical)
MAX_ASA = {"ALA": 129.0, "ARG": 274.0, "ASN": 195.0, "ASP": 193.0, "CYS": 167.0,
           "GLN": 225.0, "GLU": 223.0, "GLY": 104.0, "HIS": 224.0, "ILE": 197.0,
           "LEU": 201.0, "LYS": 236.0, "MET": 224.0, "PHE": 240.0, "PRO": 159.0,
           "SER": 155.0, "THR": 172.0, "TRP": 285.0, "TYR": 263.0, "VAL": 174.0}
# GES hydrophobicity (kcal/mol, positive = favourable in the membrane interior); the
# energy uses the negative.
GES = {"PHE": 3.7, "MET": 3.4, "ILE": 3.1, "LEU": 2.8, "VAL": 2.6, "CYS": 2.0, "TRP": 1.9,
       "ALA": 1.6, "THR": 1.2, "GLY": 1.0, "SER": 0.6, "PRO": -0.2, "TYR": -0.7, "HIS": -3.0,
       "GLN": -4.1, "ASN": -4.8, "GLU": -8.2, "LYS": -8.8, "ASP": -9.2, "ARG": -12.3}
CORE_SCALE = {k: -v for k, v in GES.items()}
INTERFACE_SCALE = WW_INTERFACE
VDW_RADIUS = {"C": 1.70, "N": 1.55, "O": 1.52, "S": 1.80, "SE": 1.90}
BACKBONE = {"N", "CA", "C", "O", "OXT"}

# -----------------------------------------------------------------------------
# Model parameters
# -----------------------------------------------------------------------------
PROBE = 1.4               # A, water probe for SASA
SR_POINTS = 96            # Shrake-Rupley test points per atom
MIN_HALF = 10.0           # A, hydrophobic half-thickness search range (20-40 A total)
MAX_HALF = 20.0
INTERFACE_WIDTH = 7.0     # A, head-group/interface region beyond the core boundary
SOFTNESS = 1.5            # A, logistic edge width of the profile
N_NORMALS = 1200          # coarse normals on the hemisphere (~5 deg spacing)
REFINE_TOP = 5            # distinct coarse minima refined locally
REFINE_CAP_DEG = 7.0      # radius of the refinement cap
REFINE_NORMALS = 160      # normals per refinement cap
PORE_ESCAPE = (0.03, 0.12)  # ray-escape fraction mapped to pore weight 0 -> 1
MAX_MEMBRANE_DG = -10.0   # kcal/mol; a best placement above this = not a membrane protein
                          # (54 test structures: membrane <= -27, soluble >= -3.6)


# -----------------------------------------------------------------------------
# Structure parsing (first model, protein residues only)
# -----------------------------------------------------------------------------
@dataclass
class Residue:
    chain: str
    resseq: int
    icode: str
    resname: str               # parent 3-letter code (aliases resolved)
    atom_names: list[str] = field(default_factory=list)
    elements: list[str] = field(default_factory=list)
    xyz: list[tuple[float, float, float]] = field(default_factory=list)

    def coord(self, name):
        try:
            return np.asarray(self.xyz[self.atom_names.index(name)], float)
        except ValueError:
            return None


def _open_text(path: Path) -> str:
    raw = Path(path).read_bytes()
    if raw[:2] == b"\x1f\x8b":
        raw = gzip.decompress(raw)
    return raw.decode("utf-8", errors="replace")


def _element(name: str, element: str, resname: str) -> str:
    e = (element or "").strip().upper()
    if e:
        return e
    n = re.sub(r"[^A-Za-z]", "", name).upper()
    if resname == "MSE" and n.startswith("SE"):
        return "SE"
    return n[:1] if n else "C"


def _add_atom(residues, index, chain, resseq, icode, resname, name, element, x, y, z):
    parent = ALIASES.get(resname, resname)
    if parent not in STANDARD:
        return                                  # water, ligands, lipids, DUM atoms ...
    el = _element(name, element, resname)
    if el in ("H", "D"):
        return
    key = (chain, resseq, icode)
    res = index.get(key)
    if res is None:
        res = Residue(chain, resseq, icode, parent)
        index[key] = res
        residues.append(res)
    if name in res.atom_names:                  # alternate location: keep the first
        return
    res.atom_names.append(name)
    res.elements.append(el)
    res.xyz.append((x, y, z))


def _parse_pdb(text: str) -> list[Residue]:
    residues, index = [], {}
    for line in text.splitlines():
        rec = line[:6]
        if rec.startswith("ENDMDL"):
            break
        if rec not in ("ATOM  ", "HETATM"):
            continue
        alt = line[16:17]
        if alt not in (" ", "A", "1"):
            continue
        try:
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            resseq = int(line[22:26])
        except ValueError:
            continue
        _add_atom(residues, index, line[21:22].strip(), resseq, line[26:27].strip(),
                  line[17:20].strip(), line[12:16].strip(), line[76:78], x, y, z)
    return residues


_CIF_TOKEN = re.compile(r"'(?:[^']|'(?!\s|$))*'|\"(?:[^\"]|\"(?!\s|$))*\"|\S+")


def _parse_mmcif(text: str) -> list[Residue]:
    lines = text.splitlines()
    residues, index = [], {}
    i = 0
    while i < len(lines):
        if lines[i].strip() == "loop_" and i + 1 < len(lines) and lines[i + 1].startswith("_atom_site."):
            cols = []
            i += 1
            while i < len(lines) and lines[i].startswith("_atom_site."):
                cols.append(lines[i].split(".", 1)[1].strip())
                i += 1
            col = {c: k for k, c in enumerate(cols)}

            def pick(*names):
                for nm in names:
                    if nm in col:
                        return col[nm]
                return None
            c_grp, c_el = pick("group_PDB"), pick("type_symbol")
            c_atom = pick("auth_atom_id", "label_atom_id")
            c_alt = pick("label_alt_id")
            c_res = pick("auth_comp_id", "label_comp_id")
            c_chain = pick("auth_asym_id", "label_asym_id")
            c_seq = pick("auth_seq_id", "label_seq_id")
            c_ins = pick("pdbx_PDB_ins_code")
            c_x, c_y, c_z = pick("Cartn_x"), pick("Cartn_y"), pick("Cartn_z")
            c_model = pick("pdbx_PDB_model_num")
            if c_x is None or c_y is None or c_z is None or c_res is None or c_chain is None or c_seq is None or c_atom is None:
                i += 1
                continue
            first_model = None
            while i < len(lines):
                line = lines[i]
                if not line.strip() or line.startswith(("_", "loop_", "#", "data_")):
                    break
                tok = [t[1:-1] if t[:1] in "'\"" else t for t in _CIF_TOKEN.findall(line)]
                i += 1
                if len(tok) < len(cols):
                    continue
                if c_model is not None:
                    if first_model is None:
                        first_model = tok[c_model]
                    elif tok[c_model] != first_model:
                        continue
                if c_grp is not None and tok[c_grp] not in ("ATOM", "HETATM"):
                    continue
                alt = tok[c_alt] if c_alt is not None else "."
                if alt not in (".", "?", "A", "1"):
                    continue
                try:
                    x, y, z = float(tok[c_x]), float(tok[c_y]), float(tok[c_z])
                    resseq = int(tok[c_seq])
                except (ValueError, TypeError):
                    continue
                icode = tok[c_ins] if c_ins is not None else ""
                icode = "" if icode in ("?", ".") else icode
                _add_atom(residues, index, tok[c_chain], resseq, icode, tok[c_res],
                          tok[c_atom], tok[c_el] if c_el is not None else "", x, y, z)
            break
        i += 1
    return residues


def load_protein_residues(path) -> list[Residue]:
    """Protein residues (with a CA) of the first model, all chains, in file order."""
    text = _open_text(Path(path))
    name = str(path).lower().replace(".gz", "")
    is_cif = name.endswith((".cif", ".mmcif")) or "_atom_site." in text[:200000]
    res = _parse_mmcif(text) if is_cif else _parse_pdb(text)
    return [r for r in res if "CA" in r.atom_names]


# -----------------------------------------------------------------------------
# Solvent accessibility
# -----------------------------------------------------------------------------
def _sphere_points(n: int) -> np.ndarray:
    k = np.arange(n) + 0.5
    phi = np.arccos(1.0 - 2.0 * k / n)
    theta = np.pi * (1.0 + 5 ** 0.5) * k
    return np.stack([np.cos(theta) * np.sin(phi), np.sin(theta) * np.sin(phi), np.cos(phi)], 1)


def shrake_rupley(xyz: np.ndarray, radii: np.ndarray, probe: float = PROBE,
                  n_points: int = SR_POINTS) -> np.ndarray:
    """Per-atom solvent-accessible surface area (A^2), Shrake & Rupley (1973).
    Neighbour search with a uniform cell grid; one vectorised block per cell."""
    xyz = np.asarray(xyz, float)
    n = len(xyz)
    if n == 0:
        return np.zeros(0)
    R = np.asarray(radii, float) + probe
    sphere = _sphere_points(n_points)
    cell = 2.0 * R.max()
    keys = np.floor((xyz - xyz.min(0)) / cell).astype(np.int64)
    dims = keys.max(0) + 1
    flat = (keys[:, 0] * dims[1] + keys[:, 1]) * dims[2] + keys[:, 2]
    order = np.argsort(flat, kind="stable")
    uniq, starts, counts = np.unique(flat[order], return_index=True, return_counts=True)
    members = {int(u): order[s:s + c] for u, s, c in zip(uniq, starts, counts)}
    offsets = [(a, b, c) for a in (-1, 0, 1) for b in (-1, 0, 1) for c in (-1, 0, 1)]
    sasa = np.zeros(n)
    for u, atoms in members.items():
        kx, rem = divmod(u, dims[1] * dims[2])
        ky, kz = divmod(rem, dims[2])
        cand = []
        for a, b, c in offsets:
            x, y, z = kx + a, ky + b, kz + c
            if 0 <= x < dims[0] and 0 <= y < dims[1] and 0 <= z < dims[2]:
                m = members.get(int((x * dims[1] + y) * dims[2] + z))
                if m is not None:
                    cand.append(m)
        cand = np.concatenate(cand)
        # true neighbours of each atom of this cell (spheres overlap), padded to kmax
        dd = np.linalg.norm(xyz[atoms, None, :] - xyz[None, cand, :], axis=-1)       # (m,k)
        nb = (dd < R[atoms, None] + R[None, cand]) & (cand[None, :] != atoms[:, None])
        kmax = int(nb.sum(1).max()) if nb.size else 0
        if kmax == 0:
            sasa[atoms] = 4.0 * np.pi * R[atoms] ** 2
            continue
        order_nb = np.argsort(~nb, axis=1, kind="stable")[:, :kmax]                  # neighbours first
        valid = np.take_along_axis(nb, order_nb, axis=1)                             # (m,kmax)
        nbr = cand[order_nb]                                                          # (m,kmax)
        pts = xyz[atoms, None, :] + R[atoms, None, None] * sphere[None, :, :]        # (m,P,3)
        d2 = ((pts[:, :, None, :] - xyz[nbr][:, None, :, :]) ** 2).sum(-1)           # (m,P,kmax)
        inside = (d2 < (R[nbr] ** 2)[:, None, :]) & valid[:, None, :]
        free = 1.0 - inside.any(-1).mean(-1)
        sasa[atoms] = 4.0 * np.pi * R[atoms] ** 2 * free
    return sasa


def residue_exposure(residues: list[Residue]):
    """(rsa, side-chain centroids, CA coords) per residue, SASA computed on all residues
    together (so subunit interfaces count as buried)."""
    xyz, radii, owner = [], [], []
    for k, r in enumerate(residues):
        for el, p in zip(r.elements, r.xyz):
            xyz.append(p); radii.append(VDW_RADIUS.get(el, 1.80)); owner.append(k)
    xyz, owner = np.asarray(xyz, float), np.asarray(owner)
    atom_sasa = shrake_rupley(xyz, np.asarray(radii))
    res_sasa = np.bincount(owner, weights=atom_sasa, minlength=len(residues))
    rsa = np.array([min(1.0, s / MAX_ASA[r.resname]) for s, r in zip(res_sasa, residues)])
    centroids, ca = [], []
    for r in residues:
        ca_r = r.coord("CA")
        side = [np.asarray(p) for nm, p in zip(r.atom_names, r.xyz) if nm not in BACKBONE]
        centroids.append(np.mean(side, axis=0) if side else ca_r)
        ca.append(ca_r)
    return rsa, np.asarray(centroids), np.asarray(ca), np.asarray(atom_sasa)


def ca_exposure(ca: np.ndarray, radius: float = 13.0) -> np.ndarray:
    """Fallback for C-alpha-only input: half-sphere exposure (HSE-up, Hamelryck 2005)
    mapped linearly to an RSA estimate. The pseudo side-chain direction at residue i is
    the sum of the unit vectors CA(i-1)->CA(i) and CA(i+1)->CA(i). The linear map was
    fitted against Shrake-Rupley RSA on 54 PDB structures (r = -0.81)."""
    ca = np.asarray(ca, float)
    n = len(ca)
    u = np.zeros((n, 3))
    for i in range(n):
        v = np.zeros(3)
        for j in (i - 1, i + 1):
            if 0 <= j < n:
                d = ca[i] - ca[j]
                dn = np.linalg.norm(d)
                if 0 < dn < 4.2:
                    v += d / dn
        nv = np.linalg.norm(v)
        u[i] = v / nv if nv > 0 else 0.0
    hse_up = np.zeros(n)
    for s in range(0, n, 512):
        diff = ca[None, :, :] - ca[s:s + 512, None, :]
        dist = np.linalg.norm(diff, axis=-1)
        up = ((diff * u[s:s + 512, None, :]).sum(-1) > 0) & (dist < radius) & (dist > 0)
        hse_up[s:s + 512] = up.sum(1)
    return np.clip(HSE_RSA_INTERCEPT + HSE_RSA_SLOPE * hse_up, 0.0, 1.0)


HSE_RSA_INTERCEPT = 0.514
HSE_RSA_SLOPE = -0.0179


# -----------------------------------------------------------------------------
# Energy and search
# -----------------------------------------------------------------------------
def _logistic(x):
    return 1.0 / (1.0 + np.exp(-np.clip(x, -50, 50)))


def profiles(z, half, width=INTERFACE_WIDTH, soft=SOFTNESS):
    az = np.abs(z)
    core = _logistic((half - az) / soft)
    outer = _logistic((half + width - az) / soft)
    return core, outer - core


def residue_energies(z, half, oct_w, if_w):
    core, inter = profiles(z, half)
    return oct_w * core + if_w * inter


def _hemisphere(n: int) -> np.ndarray:
    pts = _sphere_points(2 * n)
    return pts[pts[:, 2] >= 0][:n]


def _cap(axis: np.ndarray, radius_deg: float, n: int) -> np.ndarray:
    """n directions within radius_deg of axis (plus the axis itself)."""
    k = np.arange(n) + 0.5
    cos_min = math.cos(math.radians(radius_deg))
    cz = 1.0 - (1.0 - cos_min) * k / n
    theta = np.pi * (1.0 + 5 ** 0.5) * k
    sz = np.sqrt(1.0 - cz ** 2)
    local = np.stack([np.cos(theta) * sz, np.sin(theta) * sz, cz], 1)
    a = axis / np.linalg.norm(axis)
    t = np.array([1.0, 0, 0]) if abs(a[0]) < 0.9 else np.array([0, 1.0, 0])
    e1 = np.cross(a, t); e1 /= np.linalg.norm(e1)
    e2 = np.cross(a, e1)
    dirs = local[:, 0:1] * e1 + local[:, 1:2] * e2 + local[:, 2:3] * a
    return np.vstack([a, dirs])


def _grid_search(P, oct_w, if_w, normals, halves, bin_w, chunk=64):
    """For every normal: min over centre and half-thickness of the binned energy.
    Returns (best energy per normal, centre per normal, half per normal)."""
    rmax = float(np.linalg.norm(P, axis=1).max()) + 1.0
    nb = int(math.ceil(2 * rmax / bin_w)) + 1
    K = int(math.ceil((halves.max() + INTERFACE_WIDTH + 8 * SOFTNESS) / bin_w))
    L = 1
    while L < nb + 2 * K + 1:
        L *= 2
    u = (np.arange(L) - K) * bin_w                         # kernel sampled at offsets -K..
    kern_c, kern_i = [], []
    for h in halves:
        c, i = profiles(u, h)
        c[2 * K + 1:] = 0.0; i[2 * K + 1:] = 0.0
        kern_c.append(c); kern_i.append(i)
    FC = np.fft.rfft(np.asarray(kern_c), axis=1)           # (nH, Lf)
    FI = np.fft.rfft(np.asarray(kern_i), axis=1)
    best_e = np.empty(len(normals)); best_c = np.empty(len(normals)); best_h = np.empty(len(normals))
    for s in range(0, len(normals), chunk):
        nrm = normals[s:s + chunk]
        proj = P @ nrm.T                                   # (n, m)
        idx = np.clip(np.rint((proj + rmax) / bin_w).astype(np.int64), 0, nb - 1)
        m = nrm.shape[0]
        flat = (idx + (np.arange(m) * L)[None, :]).ravel()
        A = np.bincount(flat, weights=np.repeat(oct_w, m), minlength=m * L).reshape(m, L)
        B = np.bincount(flat, weights=np.repeat(if_w, m), minlength=m * L).reshape(m, L)
        FA, FB = np.fft.rfft(A, axis=1), np.fft.rfft(B, axis=1)
        E = np.fft.irfft(FA[:, None, :] * FC[None] + FB[:, None, :] * FI[None], n=L, axis=2)
        E = E[:, :, K:K + nb]                              # E[m, h, centre-bin]
        flatE = E.reshape(m, -1)
        k = flatE.argmin(1)
        hi, ci = np.divmod(k, nb)
        best_e[s:s + m] = flatE[np.arange(m), k]
        best_c[s:s + m] = ci * bin_w - rmax
        best_h[s:s + m] = halves[hi]
    return best_e, best_c, best_h


@dataclass
class MembranePlacement:
    normal: np.ndarray          # unit vector
    center: np.ndarray          # a point on the mid-plane (absolute coordinates)
    half_thickness: float       # hydrophobic half-thickness D (A)
    energy: float               # dG_transfer at the optimum (kcal/mol, negative = favourable)
    n_residues: int

    def depth(self, xyz) -> np.ndarray:
        return (np.asarray(xyz, float) - self.center) @ self.normal


def place_membrane(positions: np.ndarray, rsa: np.ndarray, resnames: list[str],
                   min_half: float = MIN_HALF, max_half: float = MAX_HALF,
                   core_scale: Optional[dict[str, float]] = None,
                   interface_scale: Optional[dict[str, float]] = None) -> MembranePlacement:
    """Global minimum of the implicit-membrane energy for residues at `positions`
    (side-chain centroids) with relative accessibilities `rsa`."""
    core_scale = CORE_SCALE if core_scale is None else core_scale
    interface_scale = INTERFACE_SCALE if interface_scale is None else interface_scale
    pos = np.asarray(positions, float)
    rsa = np.asarray(rsa, float)
    oct_w = rsa * np.array([core_scale.get(r, 0.0) for r in resnames])
    if_w = rsa * np.array([interface_scale.get(r, 0.0) for r in resnames])
    centroid = pos.mean(0)
    P = pos - centroid

    # 1. coarse: hemisphere normals, 1 A half-thickness steps, 0.5 A bins
    normals = _hemisphere(N_NORMALS)
    halves = np.arange(min_half, max_half + 1e-6, 1.0)
    e, c, h = _grid_search(P, oct_w, if_w, normals, halves, bin_w=0.5)

    # 2. refine the REFINE_TOP best distinct minima (>= 15 deg apart)
    order = np.argsort(e)
    seeds = []
    for k in order:
        if all(abs(float(normals[k] @ normals[j])) < math.cos(math.radians(15)) for j in seeds):
            seeds.append(k)
        if len(seeds) >= REFINE_TOP:
            break
    fine_halves = np.arange(min_half, max_half + 1e-6, 0.25)
    best = (float(np.inf), normals[0], 0.0, 0.0)
    for k in seeds:
        cap = _cap(normals[k], REFINE_CAP_DEG, REFINE_NORMALS)
        ce, cc, ch = _grid_search(P, oct_w, if_w, cap, fine_halves, bin_w=0.25)
        j = int(np.argmin(ce))
        if ce[j] < best[0]:
            best = (float(ce[j]), cap[j], float(cc[j]), float(ch[j]))

    # 3. exact (unbinned) polish of centre and half-thickness at the best normal
    _, n, c0, h0 = best
    z = P @ n
    cs = c0 + np.arange(-1.0, 1.0001, 0.1)
    hs = np.clip(h0 + np.arange(-0.5, 0.5001, 0.1), min_half, max_half)
    grid = [(float(residue_energies(z - cc, hh, oct_w, if_w).sum()), cc, hh) for cc in cs for hh in hs]
    energy, c_best, h_best = min(grid)
    return MembranePlacement(normal=n / np.linalg.norm(n), center=centroid + c_best * n,
                             half_thickness=float(h_best), energy=float(energy),
                             n_residues=len(pos))


# -----------------------------------------------------------------------------
# Pore / cavity correction
# -----------------------------------------------------------------------------
def pore_weight(escape: np.ndarray, lo: float = PORE_ESCAPE[0], hi: float = PORE_ESCAPE[1]) -> np.ndarray:
    return np.clip((np.asarray(escape, float) - lo) / (hi - lo), 0.0, 1.0)


def escape_fraction(points: np.ndarray, atom_xyz: np.ndarray, n_dirs: int = 64,
                    start: float = 4.0, length: float = 30.0, step: float = 1.0,
                    grid: float = 1.5) -> np.ndarray:
    """Fraction of straight rays from each point that leave the protein: rays start
    `start` A from the point and are blocked by any voxel (edge `grid`) holding an atom
    or next to one. Lipid-facing surface residues see roughly a half-space (~0.4-0.5);
    residues lining a water-filled pore or a deep cavity only see along the pore axis."""
    atom_xyz = np.asarray(atom_xyz, float)
    lo = atom_xyz.min(0) - 3 * grid
    dims = np.ceil((atom_xyz.max(0) + 3 * grid - lo) / grid).astype(int) + 1
    occ = np.zeros(dims, dtype=bool)
    k = np.floor((atom_xyz - lo) / grid).astype(int)
    for a in (-1, 0, 1):
        for b in (-1, 0, 1):
            for c in (-1, 0, 1):
                kk = np.clip(k + (a, b, c), 0, dims - 1)
                occ[kk[:, 0], kk[:, 1], kk[:, 2]] = True
    dirs = _sphere_points(n_dirs)
    steps = np.arange(start, length + 1e-6, step)
    out = np.empty(len(points))
    for s in range(0, len(points), 256):
        p = np.asarray(points[s:s + 256], float)
        ray = p[:, None, None, :] + dirs[None, :, None, :] * steps[None, None, :, None]   # (m,D,S,3)
        idx = np.floor((ray - lo) / grid).astype(int)
        inside = np.all((idx >= 0) & (idx < dims), axis=-1)
        idx = np.clip(idx, 0, dims - 1)
        hit = occ[idx[..., 0], idx[..., 1], idx[..., 2]] & inside
        out[s:s + 256] = 1.0 - hit.any(-1).mean(-1)
    return out


# -----------------------------------------------------------------------------
# One call for the provider
# -----------------------------------------------------------------------------
@dataclass
class EnergyFit:
    placement: MembranePlacement
    mode: str                  # "all-atom" or "C-alpha"
    n_chains: int
    warnings: list[str]


def fit_structure(file_path: Optional[Path | str] = None,
                  ca_coords: Optional[np.ndarray] = None,
                  ca_names: Optional[list[str]] = None) -> EnergyFit:
    """Place the membrane on the first model of `file_path` (all protein chains, all-atom
    SASA + pore correction). Falls back to the given C-alpha trace (HSE exposure, no pore
    correction) when the file cannot be read or holds no side chains."""
    warnings = []
    residues = []
    if file_path is not None:
        try:
            residues = load_protein_residues(file_path)
        except Exception as error:                       # unreadable / unsupported format
            warnings.append(f"could not read atoms from the file ({error}); using C-alpha only")
    has_side_chains = bool(residues) and bool(np.mean([len(r.xyz) for r in residues]) > 3.0)
    if has_side_chains:
        rsa, centroids, _, _ = residue_exposure(residues)
        atoms = np.array([p for r in residues for p in r.xyz])
        pore = pore_weight(escape_fraction(centroids, atoms))
        placement = place_membrane(centroids, rsa * pore, [r.resname for r in residues])
        return EnergyFit(placement, "all-atom", len({r.chain for r in residues}), warnings)
    if ca_coords is None:
        raise ValueError("no coordinates to place a membrane on")
    ca_names_list = ca_names if ca_names is not None else []
    names = [ALIASES.get(n, n) for n in ca_names_list]
    if residues and not has_side_chains:
        warnings.append("model has no side chains")
    warnings.append("exposure estimated from the C-alpha trace (half-sphere exposure); "
                    "no pore correction, so beta-barrel porins are unreliable")
    placement = place_membrane(np.asarray(ca_coords, float), ca_exposure(ca_coords), names)
    return EnergyFit(placement, "C-alpha", 1, warnings)