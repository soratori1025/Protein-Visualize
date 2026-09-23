"""
Core constants for the Protein-Visualize backend.
"""

# --- Kyte-Doolittle hydropathy (higher = more hydrophobic) --------------------
KYTE_DOOLITTLE = {
    "ALA": 1.8, "ARG": -4.5, "ASN": -3.5, "ASP": -3.5, "CYS": 2.5,
    "GLN": -3.5, "GLU": -3.5, "GLY": -0.4, "HIS": -3.2, "ILE": 4.5,
    "LEU": 3.8, "LYS": -3.9, "MET": 1.9, "PHE": 2.8, "PRO": -1.6,
    "SER": -0.8, "THR": -0.7, "TRP": -0.9, "TYR": -1.3, "VAL": 4.2,
}
POSITIVE_RESIDUES = {"LYS", "ARG"}  # positive-inside rule

# --- Secondary-structure code sets (shared by DSSP and STRIDE) ----------------
HELIX_CODES = {"H", "G", "I"}          # alpha, 3-10, pi
STRAND_CODES = {"E", "B", "b"}         # extended strand / beta bridge

# --- Membrane geometry default parameters (Structure Predictor) ---------------
MEMBRANE_THICKNESS = 30.0   # Angstrom, typical hydrophobic-core thickness
N_AXIS_SAMPLES = 500        # candidate membrane-normal directions
N_CENTER_SAMPLES = 160      # slab positions scanned along each candidate normal
SMOOTH_WINDOW = 19          # residues, hydropathy smoothing window (~1 TM helix)
MAX_JITTER_LEN = 3          # max length of a slab-edge jitter dip that may be merged
JITTER_MARGIN = 3.0         # Angstrom a jitter residue may sit beyond the slab face
MIN_TM_CORE = 5             # drop TM runs whose in-slab core is shorter than this
MIN_FACE_RESIDUES = 5       # min residues required OUTSIDE the slab on each face
MIN_FACE_FRACTION = 0.08    # ...and each face must hold >= this fraction of residues
MAX_SNAP = 4                # max residues a boundary may snap past the slab edge
MIN_MEMBRANE_SCORE = 0.5    # below this the structure is treated as non-membrane
MIN_EXTRA_SS_LEN = 3        # extramembrane helix/strand shorter than this -> Coil

# --- SS-element-first pipeline defaults ---------------------------------------
MIN_TM_ELEMENT_IN_SLAB = 4  # a DSSP element needs this many residues in the slab
FULL_CROSS_FRAC = 0.66      # element whose z-span >= this*thickness crosses fully
BROKEN_GAP_MAX = 9          # max intramembrane break between two halves of one helix
MIN_CROSS_SPAN_FRAC = 0.45  # a crossing must reach across at least this*thickness

# --- Sequence Predictor Constants ---------------------------------------------
TM_HYDRO_THRESHOLD = 1.6
TM_MIN_LENGTH = 15
TM_MERGE_GAP = 3
RESCUE_HYDRO_THRESHOLD = 0.6
RESCUE_MIN_ALPHA_LENGTH = 12
RESCUE_MIN_BETA_LENGTH = 5
