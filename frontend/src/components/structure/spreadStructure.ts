/**
 * "Spread out" a membrane protein for display - without falsifying its structure.
 *
 * Three modes
 * -----------
 *  aligned   ONE rigid rotation + translation of the whole file so the membrane normal
 *            is +Y (extracellular up) and the bilayer mid-plane is y = 0. Nothing is
 *            deformed: every distance, angle and contact is the real one. Use this for
 *            any structural analysis.
 *  exploded  TM bodies pushed apart in the membrane plane (native arrangement kept);
 *            a loop that cannot follow keeps its two bodies together. No deformation.
 *  unrolled  (default) 2D-like topology drawing, like a textbook topology figure:
 *              - TM crossings stand vertical (least-squares axis of all CA atoms ->
 *                normal) at their real depth, in ONE row along X at z = 0; each is spun
 *                about the normal so its loop ends lie in the drawing plane;
 *              - the gap between two TMs follows the loop between them (long loop ->
 *                more room), never less than their backbone footprints + margin;
 *              - every loop follows a planar guide curve (z = 0): out of the membrane,
 *                horizontally inside a narrow band (loopHeight above the membrane
 *                surface; long loops meander sideways instead of sticking out), back into
 *                the next TM. N/C tails run sideways the same way;
 *              - loops are shaped and re-closed by changing ONLY backbone torsions (phi/psi)
 *                of coil residues: multi-target CCD towards the curve, then exact closure
 *                (Canutescu & Dunbrack 2003). Bond lengths and bond angles are untouched,
 *                loop helices / strands / turns stay rigid, Pro phi is never rotated, side
 *                chains keep their rotamers and move with their own backbone frame;
 *              - a loop with too few flexible residues to close is drawn as a flagged
 *                "schematic junction" (report + REMARK 999), never silently stretched;
 *              - two crossings joined by < 3 flexible residues (hairpin) stay together
 *                with their native geometry, stood up by their mean axis;
 *              - ligands follow the nearest residue; waters, OPM/PPM dummy atoms and the
 *                other chains are left out.
 *            The report gives, per loop, closure, distance from its curve and clashes.
 *
 * Residues are identified by (chain, resSeq, insertion code) of the FIRST model only.
 * The drop-in signature `spreadStructure(pdbData, consensusMap)` still works; pass
 * `options` (chainId, membrane, residues from CalculatedTopologyData) for the exact
 * membrane placement computed by the backend.
 */
import type { ConsensusResidue, ResidueAnnotation } from '../../types/secondaryStructure';

/* ------------------------------------------------------------------ *
 * Public API
 * ------------------------------------------------------------------ */

export type SpreadMode = 'unrolled' | 'exploded' | 'aligned';

export interface SpreadOptions {
  /** 'unrolled' (default, 2D-like drawing), 'exploded' (TMs pushed apart in the membrane
   *  plane, native arrangement kept) or 'aligned' (exact, rigid). */
  mode?: SpreadMode;
  /** Analysed chain (CalculatedTopologyData.chain_id). Default: best match to the map. */
  chainId?: string | null;
  /** CalculatedTopologyData.membrane - exact normal used by the backend. */
  membrane?: { normal: number[]; half_thickness?: number } | null;
  /** CalculatedTopologyData.residues - per-residue depth (exact mid-plane), ss, side. */
  residues?: ResidueAnnotation[] | null;
  /** rotate each crossing so its axis is parallel to the normal
   *  (default: true for 'unrolled', false for 'exploded'). */
  straighten?: boolean;
  /** exploded: in-plane expansion factor (default: smallest that separates the bodies, <= 3). */
  explode?: number;
  /** unrolled: minimum gap between the footprints of neighbouring bodies, Å (default 6).
   *  The actual gap grows with the loop between them. */
  margin?: number;
  /** unrolled: how far loops / tails may go beyond the membrane surface, Å (default 12).
   *  Smaller = flatter drawing (narrower y range); long loops then meander sideways. */
  loopHeight?: number;
}

export interface LoopReport {
  from: string;            // last residue of the body before the loop
  to: string;              // first residue of the body after it
  residues: number;
  flexible: number;        // coil residues whose phi/psi may change
  closed: boolean;
  closureRmsd: number | null;  // Å, CCD end-point error (null = chain break, nothing to close)
  /** unrolled only: the loop could not be closed with real torsions; its last peptide
   *  bond is stretched to `stretchedBond` Å and must be drawn as schematic. */
  schematic?: boolean;
  stretchedBond?: number;
  clashes?: number;
  /** unrolled: RMS distance (Å) of the loop CAs from their planar guide curve. */
  pathDeviation?: number;
  note?: string;
}

export interface TailReport {
  tail: 'N' | 'C';
  residues: number;
  guided: number;          // residues re-shaped to lie flat (the rest move rigidly)
  pathDeviation?: number;
  clashes?: number;
  note?: string;
}

export interface SpreadReport {
  mode: SpreadMode;
  chainId: string;
  membraneSource: 'backend' | 'estimated from TM axes';
  crossings: { crossing: number; start: string; end: string; tiltRemovedDeg: number }[];
  /** groups of crossings kept together as one rigid body (native geometry). */
  keptTogether: number[][];
  loops: LoopReport[];
  tails: TailReport[];
  maxPeptideBond: number | null;     // Å, longest C-N bond between bonded residues
  brokenPeptideBonds: number;        // C-N > 1.5 Å (original < 2 Å)
  clashes: number;                   // heavy-atom pairs < 2.6 Å (residues >= 3 apart) not in the input
  brokenDisulfides: number;
  keptLigands: string[];
  dropped: string[];
  warnings: string[];
}

/** Drop-in replacement for the old function: returns the new PDB text. */
export function spreadStructure(pdbData: string, consensusMap: ConsensusResidue[],
                                options: SpreadOptions = {}): string {
  return spreadStructureWithReport(pdbData, consensusMap, options).pdb;
}

/* ------------------------------------------------------------------ *
 * Small linear algebra (rows = basis vectors)
 * ------------------------------------------------------------------ */

type V = [number, number, number];
type M3 = [V, V, V];
interface Rigid { R: M3; t: V }            // p -> R p + t

const add = (a: V, b: V): V => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
const sub = (a: V, b: V): V => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
const scale = (a: V, s: number): V => [a[0] * s, a[1] * s, a[2] * s];
const dot = (a: V, b: V) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
const cross = (a: V, b: V): V => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
const norm = (a: V) => Math.sqrt(dot(a, a));
const unit = (a: V): V => { const n = norm(a); return n > 0 ? scale(a, 1 / n) : [0, 0, 0]; };
const dist = (a: V, b: V) => norm(sub(a, b));
const I3: M3 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
const matVec = (R: M3, v: V): V => [dot(R[0], v), dot(R[1], v), dot(R[2], v)];
const transpose = (A: M3): M3 => [[A[0][0], A[1][0], A[2][0]], [A[0][1], A[1][1], A[2][1]], [A[0][2], A[1][2], A[2][2]]];
const matMul = (A: M3, B: M3): M3 => { const Bt = transpose(B); return A.map((row) => [dot(row, Bt[0]), dot(row, Bt[1]), dot(row, Bt[2])]) as M3; };
const det = (A: M3) => dot(A[0], cross(A[1], A[2]));
const rotY = (phi: number): M3 => [[Math.cos(phi), 0, Math.sin(phi)], [0, 1, 0], [-Math.sin(phi), 0, Math.cos(phi)]];
const apply = (T: Rigid, p: V): V => add(matVec(T.R, p), T.t);
const compose = (A: Rigid, B: Rigid): Rigid => ({ R: matMul(A.R, B.R), t: add(matVec(A.R, B.t), A.t) });

/** Proper rotation (det +1) by angle th about unit axis u (Rodrigues). */
function axisAngle(u: V, th: number): M3 {
  const [x, y, z] = u;
  const c = Math.cos(th), s = Math.sin(th), C = 1 - c;
  return [
    [c + x * x * C, x * y * C - z * s, x * z * C + y * s],
    [y * x * C + z * s, c + y * y * C, y * z * C - x * s],
    [z * x * C - y * s, z * y * C + x * s, c + z * z * C],
  ];
}

/** Proper rotation taking direction a onto b. (The old code returned -I for the
 *  antiparallel case: that is an inversion, det -1, i.e. a MIRROR IMAGE - L-amino
 *  acids become D and right-handed helices left-handed.) */
function rotationBetween(a: V, b: V): M3 {
  const ua = unit(a), ub = unit(b);
  const v = cross(ua, ub), s = norm(v), c = dot(ua, ub);
  if (s < 1e-9) {
    if (c > 0) return I3;
    const perp = unit(cross(ua, Math.abs(ua[0]) < 0.9 ? [1, 0, 0] : [0, 1, 0]));
    return axisAngle(perp, Math.PI);
  }
  return axisAngle(scale(v, 1 / s), Math.atan2(s, c));
}

function topEigenvector(S: M3): V {
  let v: V = unit([1, 0.7, 0.3]);
  for (let k = 0; k < 200; k++) {
    const w = unit(matVec(S, v));
    if (norm(w) === 0) break;
    if (dist(w, v) < 1e-12) { v = w; break; }
    v = w;
  }
  return v;
}

const centroid = (pts: V[]): V => scale(pts.reduce((a, b) => add(a, b), [0, 0, 0] as V), 1 / Math.max(1, pts.length));

/** Least-squares line through the points (largest principal axis). */
function principalAxis(pts: V[]): V {
  const c = centroid(pts);
  const S: M3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of pts) {
    const d = sub(p, c);
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) S[i][j] += d[i] * d[j];
  }
  return topEigenvector(S);
}

/** Orthonormal frame (rows) of three points; origin o. */
function frameOf(o: V, a: V, b: V): M3 {
  const e1 = unit(sub(a, o));
  const w = sub(b, o);
  const e2 = unit(sub(w, scale(e1, dot(w, e1))));
  return [e1, e2, cross(e1, e2)];
}

/** Rigid transform taking frame (o0,a0,b0) onto (o1,a1,b1). */
function frameTransform(o0: V, a0: V, b0: V, o1: V, a1: V, b1: V): Rigid {
  const R = matMul(transpose(frameOf(o1, a1, b1)), frameOf(o0, a0, b0));
  return { R, t: sub(o1, matVec(R, o0)) };
}

const median = (xs: number[]) => {
  const s = [...xs].sort((a, b) => a - b);
  return s.length ? (s.length % 2 ? s[(s.length - 1) / 2] : (s[s.length / 2 - 1] + s[s.length / 2]) / 2) : 0;
};

/* ------------------------------------------------------------------ *
 * PDB parsing
 * ------------------------------------------------------------------ */

const ONE: Record<string, string> = {
  ALA: 'A', ARG: 'R', ASN: 'N', ASP: 'D', CYS: 'C', GLN: 'Q', GLU: 'E', GLY: 'G', HIS: 'H', ILE: 'I',
  LEU: 'L', LYS: 'K', MET: 'M', PHE: 'F', PRO: 'P', SER: 'S', THR: 'T', TRP: 'W', TYR: 'Y', VAL: 'V',
  // modified / force-field names -> parent residue
  MSE: 'M', SEP: 'S', TPO: 'T', PTR: 'Y', HYP: 'P', MLY: 'K', CSO: 'C', SEC: 'C', PYL: 'K', FME: 'M',
  HSD: 'H', HSE: 'H', HSP: 'H', HID: 'H', HIE: 'H', HIP: 'H', CYX: 'C', CYM: 'C', ASH: 'D', GLH: 'E', LYN: 'K',
};
const WATER = new Set(['HOH', 'WAT', 'DOD', 'H2O', 'TIP', 'TIP3', 'SOL']);
const DUMMY = new Set(['DUM']);                     // OPM / PPM / memembed membrane planes
const N_H_NAMES = new Set(['H', 'HN', 'H1', 'H2', 'H3']);
const O_NAMES = new Set(['O', 'OXT', 'OT1', 'OT2']);
const CLASH = 2.6;             // Å; heavy atoms of residues >= 3 apart closer than this overlap
                               //    (backbone H-bonds are >= 2.8 Å, so they never count)
const CCD_RESTARTS = 64;       // max CCD attempts per loop (first from the native conformation)
const CCD_REFINE = 200;        // Monte Carlo steps that remove the remaining clashes
// 2D ('unrolled') drawing
const STEP_COIL = 3.3;         // Å of drawing per coil residue (CA-CA 3.8, not fully extended)
const STEP_HELIX = 1.5;        // Å per helix residue (rise along the axis)
const LOOP_SPAN_PER_A = 0.4;   // horizontal room given to a loop per Å of its contour
const LOOP_SPAN_MAX = 60;      // Å, widest gap a loop can open between two TMs
const TAIL_SPAN_MAX = 40;      // Å, how far a tail runs sideways
const TAIL_GUIDE = 40;         // tail residues next to the membrane that are laid flat
const END_WEIGHT = 20;
const GUIDED_RESTARTS = 24;    // attempts per loop in the 2D layout
const GUIDED_REFINE = 80;      // clash-removal steps per loop in the 2D layout
const BACKBONE_NAMES = new Set(['N', 'CA', 'C', 'O']);
const BACKBONE_PAD = 2;        // Å added to the backbone footprint of a TM in the 2D layout         // closure vs path when guiding a loop
const GOOD_PATH_DEV = 2.5;     // Å, stop searching once the loop follows its curve this well

interface Atom { li: number; rec: string; name: string; alt: string; resName: string; chain: string;
                 resSeq: number; iCode: string; p: V; model: number }
interface Res { chain: string; resSeq: number; iCode: string; resName: string; atoms: Atom[];
                N?: Atom; CA?: Atom; C?: Atom }

const resLabel = (r: { chain: string; resSeq: number; iCode: string }) => `${r.chain}:${r.resSeq}${r.iCode}`;

function parsePdb(text: string) {
  const lines = text.split(/\r?\n/);
  const atoms: Atom[] = [];
  let model = 0, seenModel = false;
  lines.forEach((line, li) => {
    if (line.startsWith('MODEL')) { if (seenModel) model++; seenModel = true; return; }
    if (!(line.startsWith('ATOM  ') || line.startsWith('HETATM'))) return;
    const x = parseFloat(line.substring(30, 38)), y = parseFloat(line.substring(38, 46)), z = parseFloat(line.substring(46, 54));
    const resSeq = parseInt(line.substring(22, 26), 10);
    if (!isFinite(x) || !isFinite(y) || !isFinite(z) || isNaN(resSeq)) return;
    atoms.push({
      li, rec: line.substring(0, 6).trim(), name: line.substring(12, 16).trim(), alt: line.substring(16, 17),
      resName: line.substring(17, 20).trim(), chain: line.substring(21, 22), resSeq,
      iCode: line.substring(26, 27).trim(), p: [x, y, z], model,
    });
  });
  // residues of the first model, file order, keyed by chain + number + insertion code
  const residues: Res[] = [];
  const byKey = new Map<string, Res>();
  for (const a of atoms) {
    if (a.model !== 0) continue;
    const key = `${a.chain}|${a.resSeq}|${a.iCode}`;
    let r = byKey.get(key);
    if (!r) { r = { chain: a.chain, resSeq: a.resSeq, iCode: a.iCode, resName: a.resName, atoms: [] }; byKey.set(key, r); residues.push(r); }
    r.atoms.push(a);
    if (a.name === 'N') r.N = pickAlt(r.N, a);
    if (a.name === 'CA') r.CA = pickAlt(r.CA, a);
    if (a.name === 'C') r.C = pickAlt(r.C, a);
  }
  return { lines, atoms, residues };
}

/** Backbone atom of the main conformer: blank altLoc, else 'A' / '1', else the first. */
const preferredAlt = (a: Atom) => a.alt === ' ' || a.alt === 'A' || a.alt === '1';
const pickAlt = (cur: Atom | undefined, a: Atom) => (!cur || (!preferredAlt(cur) && preferredAlt(a)) ? a : cur);

const STANDARD = new Set(['ALA', 'ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'GLY', 'HIS', 'ILE', 'LEU', 'LYS',
  'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']);
/** Polymer residue: ATOM record, or a HETATM modified residue (MSE ...). A standard
 *  amino acid written as HETATM is a free ligand, not part of the chain. */
const isAminoAcid = (r: Res) => !!r.CA && r.resName in ONE &&
  (r.atoms[0].rec === 'ATOM' || !STANDARD.has(r.resName));

/* ------------------------------------------------------------------ *
 * Main
 * ------------------------------------------------------------------ */

export function spreadStructureWithReport(pdbData: string, consensusMap: ConsensusResidue[],
                                          options: SpreadOptions = {}): { pdb: string; report: SpreadReport } {
  const mode: SpreadMode = options.mode ?? 'unrolled';
  const straighten = options.straighten ?? mode === 'unrolled';
  const margin = options.margin ?? 6;
  const loopHeight = options.loopHeight ?? 12;
  const half = options.membrane?.half_thickness ?? 15;
  const warnings: string[] = [];
  const { lines, atoms, residues: allRes } = parsePdb(pdbData);

  // ---- 1. analysed chain: the one whose residue numbers AND amino acids match the map
  const chains = [...new Set(allRes.filter(isAminoAcid).map((r) => r.chain))];
  const mapKey = (n: number, ic?: string | null) => `${n}|${(ic ?? '').trim()}`;
  let chainId = options.chainId ?? '';
  if (!chains.includes(chainId)) {
    if (options.chainId) warnings.push(`chain '${options.chainId}' not in the file; best match used`);
    let best = -1;
    for (const ch of chains) {
      const have = new Map(allRes.filter((r) => r.chain === ch && isAminoAcid(r)).map((r) => [mapKey(r.resSeq, r.iCode), r]));
      const score = consensusMap.filter((m) => have.get(mapKey(m.residue_number, m.insertion_code))
        && ONE[have.get(mapKey(m.residue_number, m.insertion_code))!.resName] === m.aa).length;
      if (score > best) { best = score; chainId = ch; }
    }
  }
  const res = allRes.filter((r) => r.chain === chainId && isAminoAcid(r));
  const n = res.length;
  const posOf = new Map(res.map((r, i) => [mapKey(r.resSeq, r.iCode), i]));
  const empty = (why: string) => ({
    pdb: pdbData,
    report: { mode, chainId, membraneSource: 'estimated from TM axes' as const, crossings: [], keptTogether: [],
              loops: [], tails: [], maxPeptideBond: null, brokenPeptideBonds: 0, clashes: 0, brokenDisulfides: 0,
              keptLigands: [], dropped: [],
              warnings: [...warnings, why] },
  });
  if (n === 0) return empty('no amino-acid residues found for the analysed chain');

  // ---- 2. crossings from the residue map (by crossing; old maps: by tm_segment)
  const hasCrossing = consensusMap.some((m) => m.crossing != null);
  const crossRanges = new Map<number, [number, number]>();
  let unmatched = 0;
  for (const m of consensusMap) {
    if (m.label !== 'TM_in') continue;
    const id = hasCrossing ? m.crossing : m.tm_segment;
    if (id == null) continue;
    const i = posOf.get(mapKey(m.residue_number, m.insertion_code));
    if (i === undefined) { unmatched++; continue; }
    const cur = crossRanges.get(id);
    crossRanges.set(id, cur ? [Math.min(cur[0], i), Math.max(cur[1], i)] : [i, i]);
  }
  if (!hasCrossing) warnings.push('residue map has no crossing numbers (old backend): grouped by TM segment - a hairpin under one segment is treated as one body');
  if (unmatched) warnings.push(`${unmatched} TM_in residues of the map are not in chain ${chainId} of this file`);
  const bodies = [...crossRanges.entries()].map(([id, [s, e]]) => ({ id, s, e })).sort((a, b) => a.s - b.s);
  if (bodies.length === 0) return empty('no TM crossing in the residue map');

  // ---- 3. membrane frame: Y = normal (extracellular up), y = 0 at the mid-plane
  const caOf = (i: number) => res[i].CA!.p;
  const annot = new Map((options.residues ?? []).map((r) => [mapKey(r.residue_number, r.insertion_code), r]));
  let normal: V;
  let membraneSource: SpreadReport['membraneSource'];
  if (options.membrane?.normal?.length === 3) {
    normal = unit(options.membrane.normal as V);
    membraneSource = 'backend';
  } else {
    const S: M3 = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
    for (const b of bodies) {
      const ax = principalAxis(res.slice(b.s, b.e + 1).map((r) => r.CA!.p));
      for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) S[i][j] += ax[i] * ax[j];
    }
    normal = topEigenvector(S);
    membraneSource = 'estimated from TM axes';
  }
  const depthRows = res.map((r, i) => ({ i, d: annot.get(mapKey(r.resSeq, r.iCode))?.depth }))
    .filter((x): x is { i: number; d: number } => typeof x.d === 'number');
  let offset = membraneSource === 'backend' && depthRows.length
    ? median(depthRows.map(({ i, d }) => dot(normal, caOf(i)) - d))
    : median(bodies.map((b) => dot(normal, centroid(res.slice(b.s, b.e + 1).map((r) => r.CA!.p)))));
  // extracellular up: evidence from the backend labels, else TM_E / TM_C of the map
  const ext: number[] = [], cyt: number[] = [];
  res.forEach((r, i) => {
    const lab = annot.get(mapKey(r.resSeq, r.iCode))?.label ?? '';
    if (lab.startsWith('Extracellular')) ext.push(dot(normal, caOf(i)) - offset);
    if (lab.startsWith('Cytoplasmic')) cyt.push(dot(normal, caOf(i)) - offset);
  });
  if (!ext.length && !cyt.length) {
    for (const m of consensusMap) {
      const i = posOf.get(mapKey(m.residue_number, m.insertion_code));
      if (i === undefined) continue;
      if (m.label === 'TM_E' || m.label === 'Turn_E') ext.push(dot(normal, caOf(i)) - offset);
      if (m.label === 'TM_C' || m.label === 'Turn_C') cyt.push(dot(normal, caOf(i)) - offset);
    }
  }
  const mean = (xs: number[]) => xs.reduce((a, b) => a + b, 0) / xs.length;
  if (ext.length && cyt.length && mean(ext) < mean(cyt)) { normal = scale(normal, -1); offset = -offset; }
  else if (!ext.length || !cyt.length) warnings.push('no side information: extracellular may be drawn down');
  // in-plane X = longest spread of the crossing centroids (side view along the protein)
  const cents = bodies.map((b) => centroid(res.slice(b.s, b.e + 1).map((r) => r.CA!.p)));
  const inPlane = (v: V): V => sub(v, scale(normal, dot(v, normal)));
  let ex: V = cents.length > 1 ? unit(inPlane(principalAxis(cents))) : [0, 0, 0];
  if (norm(ex) < 0.5) ex = unit(inPlane(Math.abs(normal[0]) < 0.9 ? [1, 0, 0] : [0, 0, 1]));
  const firstToLast = sub(cents[cents.length - 1], cents[0]);
  if (dot(ex, firstToLast) < 0) ex = scale(ex, -1);           // TM1 on the left
  let ez = cross(ex, normal);
  const Rm: M3 = [ex, normal, ez];
  if (det(Rm) < 0) { ez = scale(ez, -1); Rm[2] = ez; }
  const c0 = centroid(res.map((r) => r.CA!.p));
  const toMembrane: Rigid = { R: Rm, t: [-dot(ex, c0), -offset, -dot(ez, c0)] };

  const reportBase = {
    mode, chainId, membraneSource,
    crossings: [] as SpreadReport['crossings'], keptTogether: [] as number[][], loops: [] as LoopReport[],
    tails: [] as TailReport[],
    maxPeptideBond: null as number | null, brokenPeptideBonds: 0, clashes: 0, brokenDisulfides: 0,
    keptLigands: [] as string[], dropped: [] as string[], warnings,
  };

  // ---- aligned mode: one rigid transform of everything (all chains, models, ligands)
  if (mode === 'aligned') {
    const out = lines.flatMap((line) => {
      if (line.startsWith('ANISOU')) return [];                 // tensors would need rotating
      if (!(line.startsWith('ATOM  ') || line.startsWith('HETATM'))) return [line];
      const p = apply(toMembrane, [parseFloat(line.substring(30, 38)), parseFloat(line.substring(38, 46)), parseFloat(line.substring(46, 54))]);
      return [writeCoords(line, p)];
    });
    insertRemarks(out, [
      'SPREAD MODE: ALIGNED - one rigid rotation/translation, no deformation.',
      'Y = membrane normal (extracellular up), y = 0 at the bilayer mid-plane.',
    ]);
    reportBase.crossings = bodies.map((b) => ({ crossing: b.id, start: resLabel(res[b.s]), end: resLabel(res[b.e]), tiltRemovedDeg: 0 }));
    return { pdb: out.join('\n'), report: reportBase };
  }

  // ---- 4. unrolled: everything in membrane coordinates
  const q = new Map<Atom, V>();
  for (const a of atoms) if (a.model === 0) q.set(a, apply(toMembrane, a.p));
  const Q = (a: Atom) => q.get(a)!;
  const bonded = (i: number) => {                  // residue i bonded to i-1 in the model?
    if (i <= 0) return false;
    const c = res[i - 1].C, nn = res[i].N;
    if (c && nn) return dist(c.p, nn.p) < 2.0;
    return dist(res[i - 1].CA!.p, res[i].CA!.p) < 4.2;
  };
  const rigidByMap = new Set<number>();
  for (const m of consensusMap) {
    const i = posOf.get(mapKey(m.residue_number, m.insertion_code));
    if (i !== undefined && !m.label.startsWith('Turn')) rigidByMap.add(i);
  }
  const ss = (i: number): 'H' | 'E' | 'C' => {
    const a = annot.get(mapKey(res[i].resSeq, res[i].iCode));
    if (a) return a.ss === 'H' || a.ss === 'E' ? a.ss : 'C';
    return rigidByMap.has(i) ? 'H' : 'C';
  };
  const hasBackbone = (i: number) => !!(res[i].N && res[i].CA && res[i].C);
  const flexible = (i: number) => ss(i) === 'C' && hasBackbone(i);

  // body geometry: straightened (least-squares axis -> normal) or native
  interface Group { members: number[]; first: number; last: number; R: M3; c: V }
  const range = (s: number, e: number) => Array.from({ length: Math.max(0, e - s + 1) }, (_, k) => s + k);
  const makeGroup = (members: number[], first: number, last: number, straight: boolean): Group => {
    const cas = res.slice(first, last + 1).map((r) => Q(r.CA!));
    const c = centroid(cas);
    let R = I3;
    if (straight) {
      let u = principalAxis(cas);
      if (dot(sub(cas[cas.length - 1], cas[0]), u) < 0) u = scale(u, -1);
      R = rotationBetween(u, [0, u[1] >= 0 ? 1 : -1, 0]);
    }
    return { members, first, last, R, c };
  };
  // tilt of each crossing against the normal (what straightening removes)
  const tiltOf = new Map<number, number>(bodies.map((b) => {
    const u = principalAxis(res.slice(b.s, b.e + 1).map((r) => Q(r.CA!)));
    return [b.id, straighten ? Math.acos(Math.min(1, Math.abs(u[1]))) * 180 / Math.PI : 0];
  }));
  const merge = (g: Group, h: Group) => {
    const out = makeGroup([...g.members, ...h.members], g.first, h.last, false);
    if (!straighten) return out;
    // stand the group up by the mean axis of its crossings (antiparallel ones flipped),
    // so a hairpin that cannot be pulled apart is still drawn vertical
    let sum: V = [0, 0, 0], ref: V | null = null;
    for (const b of bodies.filter((x) => out.members.includes(x.id))) {
      let u = principalAxis(res.slice(b.s, b.e + 1).map((r) => Q(r.CA!)));
      if (ref && dot(u, ref) < 0) u = scale(u, -1);
      ref = ref ?? u;
      sum = add(sum, u);
    }
    if (norm(sum) > 0) { const u = unit(sum); out.R = rotationBetween(u, [0, u[1] >= 0 ? 1 : -1, 0]); }
    return out;
  };

  // start: every crossing is its own body, except where the loop between two crossings
  // has < 3 flexible residues (and no chain break): it cannot let them apart
  let groups: Group[] = [];
  for (const b of bodies) {
    const g = groups[groups.length - 1];
    if (g) {
      const loop = range(g.last + 1, b.s - 1);
      const gap = [...loop, b.s].some((i) => !bonded(i));
      if (!gap && loop.filter(flexible).length < 3) {
        groups[groups.length - 1] = merge(g, makeGroup([b.id], b.s, b.e, false));  // straightened as a group
        continue;
      }
    }
    groups.push(makeGroup([b.id], b.s, b.e, straighten));
  }

  // exploded: a loop that cannot be closed without a clash keeps its two bodies together
  // (native geometry) and the layout is redone. unrolled: the row is kept; such a loop is
  // placed as well as CCD can and flagged as a schematic junction.
  let result = layout(groups);
  for (let round = 0; mode === 'exploded' && result.failed.length && round <= bodies.length; round++) {
    const next: Group[] = [];
    groups.forEach((g, k) => {
      if (k > 0 && result.failed.includes(k - 1)) next[next.length - 1] = merge(next[next.length - 1], g);
      else next.push(g);
    });
    groups = next;
    result = layout(groups);
  }
  const { placed, resT, loops, tails } = result;

  function layout(gs: Group[]) {
    const flat = mode === 'unrolled';
    // footprint radius around each body's vertical axis. exploded: the first / last body
    // carry the N / C tail rigidly, so the tails count; unrolled: tails are laid out flat
    const span = (k: number): [number, number] => flat ? [gs[k].first, gs[k].last]
      : [k === 0 ? 0 : gs[k].first, k === gs.length - 1 ? n - 1 : gs[k].last];
    // unrolled (a cartoon): the backbone footprint + BACKBONE_PAD - side chains of two
    // neighbouring TMs may interleave, as they do in the membrane
    const radius = gs.map((g, k) => {
      let r = 0;
      const [s0, e0] = span(k);
      for (let i = s0; i <= e0; i++) for (const a of res[i].atoms) {
        if (a.model !== 0 || (flat && !BACKBONE_NAMES.has(a.name))) continue;
        const p = matVec(g.R, sub(Q(a), g.c));
        r = Math.max(r, Math.hypot(p[0], p[2]));
      }
      return flat ? r + BACKBONE_PAD : r;
    });
    // spin phi about the normal (Y) is free: it does not change depth or tilt
    const bodyT = (g: Group, target: V, phi: number): Rigid => {
      const R = matMul(rotY(phi), g.R);
      return { R, t: sub(target, matVec(R, g.c)) };
    };
    const T: Rigid[] = [];
    if (mode === 'exploded') {
      // push every body outwards in the membrane plane; arrangement (who neighbours whom) kept
      let f = options.explode ?? 1;
      if (options.explode == null) {
        for (let i = 0; i < gs.length; i++) for (let j = i + 1; j < gs.length; j++) {
          const d = Math.hypot(gs[i].c[0] - gs[j].c[0], gs[i].c[2] - gs[j].c[2]);
          if (d > 1e-6) f = Math.max(f, (radius[i] + radius[j] + margin) / d);
        }
        f = Math.min(f, 3);
      }
      gs.forEach((g) => T.push(bodyT(g, [g.c[0] * f, g.c[1], g.c[2] * f], 0)));
    } else {
      // one row along X at z = 0. The gap to the next body follows the loop between them
      // (long loop -> more room to lie flat), never less than the two footprints + margin.
      const xs: number[] = [0];
      for (let k = 1; k < gs.length; k++) {
        const loop = range(gs[k - 1].last + 1, gs[k].first - 1);
        const want = Math.min(LOOP_SPAN_MAX, LOOP_SPAN_PER_A * contour(loop));
        xs.push(xs[k - 1] + Math.max(radius[k - 1] + radius[k] + margin, want));
      }
      const shift = (xs[0] + xs[xs.length - 1]) / 2;
      gs.forEach((g, k) => {
        const target: V = [xs[k] - shift, g.c[1], 0];
        const entry = Q(res[g.first].N ?? res[g.first].CA!);
        const exit = Q(res[g.last].C ?? res[g.last].CA!);
        const prevExit = k ? apply(T[k - 1], Q(res[gs[k - 1].last].C ?? res[gs[k - 1].last].CA!)) : null;
        const memberCentres = bodies.filter((b) => g.members.includes(b.id))
          .map((b) => centroid(res.slice(b.s, b.e + 1).map((r) => Q(r.CA!))));
        let best = { phi: 0, cost: Infinity };
        for (let step = 0; step < 72; step++) {
          const phi = (step / 72) * 2 * Math.PI;
          const Tk = bodyT(g, target, phi);
          const en = apply(Tk, entry), ex = apply(Tk, exit);
          // loop ends in the drawing plane (z = 0), entry on the left, exit on the right,
          // and the crossings of a kept-together group spread along X, not in depth
          // + the loop from the previous body must be able to reach this entry
          let cost = (prevExit ? dist(en, prevExit) : 0)
            + 0.3 * (Math.abs(en[2]) + Math.abs(ex[2])) - 0.1 * (ex[0] - en[0]);
          for (const mc of memberCentres) cost += 0.3 * Math.abs(apply(Tk, mc)[2]);
          if (cost < best.cost) best = { phi, cost };
        }
        T.push(bodyT(g, target, best.phi));
      });
    }

    const resT: (Rigid | null)[] = new Array(n).fill(null);
    const placed = new Map<Atom, V>();
    const placeResidue = (i: number, Tr: Rigid) => {
      resT[i] = Tr;
      for (const a of res[i].atoms) if (a.model === 0) placed.set(a, apply(Tr, Q(a)));
    };
    gs.forEach((_, k) => { const [s0, e0] = span(k); for (let i = s0; i <= e0; i++) placeResidue(i, T[k]); });
    const gridWithout = (lo: number, hi: number) => {
      const grid = new Grid(3.0);
      for (let i = 0; i < n; i++) {
        if (!resT[i] || (i >= lo && i <= hi)) continue;
        for (const a of res[i].atoms) if (a.model === 0 && !isHydrogen(a)) grid.add(placed.get(a)!, i);
      }
      return grid;
    };
    const sideOf = (y: number) => (y >= 0 ? 1 : -1);

    const loops: LoopReport[] = [];
    const failed: number[] = [];
    for (let k = 0; k + 1 < gs.length; k++) {
      const g = gs[k], h = gs[k + 1];
      const loopIdx = range(g.last + 1, h.first - 1);
      const entry: LoopReport = { from: resLabel(res[g.last]), to: resLabel(res[h.first]), residues: loopIdx.length,
                                  flexible: loopIdx.filter(flexible).length, closed: false, closureRmsd: null };
      loops.push(entry);
      const breakAt = [...loopIdx, h.first].find((i) => !bonded(i));
      if (breakAt !== undefined) {
        for (const i of loopIdx) placeResidue(i, i < breakAt ? T[k] : T[k + 1]);
        entry.closed = true;
        entry.note = `chain break before ${resLabel(res[breakAt])} (missing residues): nothing to close`;
        continue;
      }
      if (!loopIdx.every(hasBackbone) || !hasBackbone(h.first) || !res[g.last].C) {
        for (const i of loopIdx) placeResidue(i, T[k]);
        entry.note = 'no full N/CA/C backbone: cannot be rebuilt';
        entry.schematic = mode === 'unrolled';
        failed.push(k);
        continue;
      }
      const ch = buildChain(loopIdx, false, T[k], h.first);
      const F = [res[h.first].N!, res[h.first].CA!, res[h.first].C!].map((a) => apply(T[k + 1], Q(a)));
      let targets: V[] | null = null;
      if (flat) {
        // planar guide: out of the membrane, along a narrow band, back into the next TM
        const A = placed.get(res[g.last].C!)!, B = F[0];
        targets = residueTargets(loopIdx, guideCurve(A, B, sideOf(A[1] + B[1]), contour(loopIdx) + STEP_COIL,
                                                     1, half + 3, half + loopHeight, LOOP_SPAN_MAX), true);
      }
      const best = solveChain(ch, targets, F, gridWithout(loopIdx[0], loopIdx[loopIdx.length - 1]), placed, 1000 + k);
      entry.closureRmsd = Math.round(best.rmsd * 100) / 100;
      entry.closed = best.rmsd <= 0.3;
      entry.clashes = Number.isFinite(best.clashes) ? best.clashes : undefined;
      if (targets) entry.pathDeviation = Math.round(best.dev * 10) / 10;
      const ok = entry.closed && best.clashes === 0;
      if (!ok) failed.push(k);
      if (!ok && mode === 'exploded') {                 // will be kept together next round
        for (const i of loopIdx) placeResidue(i, T[k]);
        continue;
      }
      if (!entry.closed && flat && targets && loopIdx.length >= 2) {
        // cannot be closed with real torsions. Alternative: grow the two halves flat along
        // the curve, one from each TM, so the (flagged) gap sits at the top of the arch.
        // Kept only if that gap is smaller than the stretched bond of the closure attempt.
        const endGap = dist(backbone(ch, best.P, loopIdx.length - 1).C, F[0]);
        const tryPlaced = new Map(placed), tryResT = [...resT];
        const mid = Math.floor(loopIdx.length / 2);
        const left = loopIdx.slice(0, mid), right = loopIdx.slice(mid).reverse();
        const gridFor = (lo: number, hi: number) => {
          const grid = new Grid(3.0);
          for (let i = 0; i < n; i++) {
            if (!tryResT[i] || (i >= lo && i <= hi)) continue;
            for (const a of res[i].atoms) if (a.model === 0 && !isHydrogen(a)) grid.add(tryPlaced.get(a)!, i);
          }
          return grid;
        };
        if (left.length) {
          const chL = buildChain(left, false, T[k]);
          const bL = solveChain(chL, targets.slice(0, mid), null, gridFor(loopIdx[0], loopIdx[loopIdx.length - 1]), tryPlaced, 2000 + k);
          placeChain(chL, bL.P, tryPlaced, tryResT);
        }
        const chR = buildChain(right, true, T[k + 1]);
        const bR = solveChain(chR, targets.slice(mid).reverse(), null, gridFor(loopIdx[mid], loopIdx[loopIdx.length - 1]), tryPlaced, 3000 + k);
        placeChain(chR, bR.P, tryPlaced, tryResT);
        const gapRes = loopIdx[mid - 1] ?? g.last;
        const midGap = dist(tryPlaced.get(res[gapRes].C!)!, tryPlaced.get(res[loopIdx[mid]].N!)!);
        if (midGap < endGap) {
          for (const [atom, p2] of tryPlaced) placed.set(atom, p2);
          tryResT.forEach((t2, i) => { resT[i] = t2; });
          entry.schematic = true;
          entry.stretchedBond = Math.round(midGap * 10) / 10;
          entry.pathDeviation = Math.round(Math.sqrt(loopIdx.reduce((s2, i, j) =>
            s2 + dist(placed.get(res[i].CA!)!, targets![j]) ** 2, 0) / loopIdx.length) * 10) / 10;
          entry.note = `schematic junction: too few flexible residues to close; gap of ${entry.stretchedBond} Å between ${resLabel(res[gapRes])} and ${resLabel(res[loopIdx[mid]])}`;
          continue;
        }
      }
      placeChain(ch, best.P, placed, resT);
      if (!entry.closed) {
        const lastC = placed.get(res[loopIdx[loopIdx.length - 1]].C!)!;
        entry.schematic = true;
        entry.stretchedBond = Math.round(dist(lastC, F[0]) * 10) / 10;
        entry.note = `schematic junction: the peptide bond into ${entry.to} is stretched to ${entry.stretchedBond} Å`;
      } else if (best.clashes > 0) {
        entry.note = `closed, but ${best.clashes} heavy-atom contact(s) < ${CLASH} Å remain`;
      }
    }

    // unrolled: N / C tails laid flat too (outwards, same narrow band). Only the TAIL_GUIDE
    // residues next to the membrane are re-shaped; anything beyond (a folded domain) moves
    // rigidly with the last of them and keeps its own structure.
    const tails: TailReport[] = [];
    const first = gs[0], last = gs[gs.length - 1];
    const doTail = (which: 'N' | 'C') => {
      if (!flat) return;                                // exploded: tails ride with their TM
      const reverse = which === 'N';
      const all = reverse ? range(0, first.first - 1).reverse() : range(last.last + 1, n - 1);
      if (!all.length) return;
      const anchorT = reverse ? T[0] : T[T.length - 1];
      const report: TailReport = { tail: which, residues: all.length, guided: 0 };
      tails.push(report);
      // guided part: contiguous (bonded) residues with a full backbone next to the body
      const guided: number[] = [];
      for (const i of all) {
        const linked = reverse ? bonded(i + 1) : bonded(i);
        if (!linked || !hasBackbone(i) || guided.length >= TAIL_GUIDE) break;
        guided.push(i);
      }
      const anchorRes = reverse ? first.first : last.last;
      if (!guided.length || !hasBackbone(anchorRes)) {
        for (const i of all) placeResidue(i, anchorT);
        report.note = 'kept rigid with its TM (chain break or no backbone)';
        return;
      }
      const ch = buildChain(guided, reverse, anchorT);
      const A = placed.get(reverse ? res[anchorRes].N! : res[anchorRes].C!)!;
      const targets = residueTargets(guided, guideCurve(A, null, sideOf(A[1]), contour(guided),
                                                        reverse ? -1 : 1, half + 3, half + loopHeight, TAIL_SPAN_MAX), false);
      const lo = Math.min(...guided), hi = Math.max(...guided);
      const best = solveChain(ch, targets, null, gridWithout(Math.min(lo, ...all), Math.max(hi, ...all)), placed, reverse ? 7 : 11);
      placeChain(ch, best.P, placed, resT);
      const tip = guided[guided.length - 1];
      for (const i of all.slice(guided.length)) placeResidue(i, resT[tip]!);
      report.guided = guided.length;
      report.pathDeviation = Math.round(best.dev * 10) / 10;
      report.clashes = best.clashes;
      if (guided.length < all.length) report.note = `${all.length - guided.length} residue(s) beyond move rigidly with ${resLabel(res[tip])}`;
    };
    doTail('N');
    doTail('C');
    return { placed, resT, loops, failed, tails };
  }

  /* ---- kinematic chains: N/CA/C of consecutive residues, rotatable phi/psi ---- */

  interface Tor { a: number; b: number; from: number }
  /** Residues `idx` in order AWAY from the fixed body (reverse = an N-terminal tail,
   *  growing towards lower residue numbers). 3 points per residue ([N,CA,C], or [C,CA,N]
   *  when reversed) + the next body's N/CA/C when the chain must close onto it. */
  interface Chain { idx: number[]; reverse: boolean; start: V[]; tors: Tor[]; end: boolean }

  function buildChain(idx: number[], reverse: boolean, T0: Rigid, endRes?: number): Chain {
    const start: V[] = [];
    for (const i of idx) {
      const r = res[i];
      for (const a of reverse ? [r.C!, r.CA!, r.N!] : [r.N!, r.CA!, r.C!]) start.push(apply(T0, Q(a)));
    }
    if (endRes !== undefined) for (const a of [res[endRes].N!, res[endRes].CA!, res[endRes].C!]) start.push(apply(T0, Q(a)));
    const tors: Tor[] = [];
    idx.forEach((i, j) => {
      if (!flexible(i)) return;                        // helices / strands stay rigid
      const pro = res[i].resName === 'PRO';            // the ring fixes Pro phi
      if (!reverse) {
        if (!pro) tors.push({ a: 3 * j, b: 3 * j + 1, from: 3 * j + 2 });       // phi  N-CA
        tors.push({ a: 3 * j + 1, b: 3 * j + 2, from: 3 * j + 3 });             // psi  CA-C
      } else {
        tors.push({ a: 3 * j, b: 3 * j + 1, from: 3 * j + 2 });                 // psi  C-CA
        if (!pro) tors.push({ a: 3 * j + 1, b: 3 * j + 2, from: 3 * j + 3 });   // phi  CA-N
      }
    });
    return { idx, reverse, start, tors: tors.filter((t) => t.from < start.length), end: endRes !== undefined };
  }

  function backbone(ch: Chain, P: V[], j: number) {
    return ch.reverse ? { N: P[3 * j + 2], CA: P[3 * j + 1], C: P[3 * j] } : { N: P[3 * j], CA: P[3 * j + 1], C: P[3 * j + 2] };
  }

  /** Rigid transform of every chain residue from its own backbone frame: side chain and
   *  N/CA/C from (CA, C, N), carbonyl O from (C, CA, N+1), amide H from (N, CA, C-1). */
  function chainTransforms(ch: Chain, P: V[], placed: Map<Atom, V>) {
    const m = ch.idx.length;
    const placedAt = (a?: Atom) => (a && placed.has(a) ? placed.get(a)! : null);
    return ch.idx.map((i, j) => {
      const r = res[i];
      const { N, CA, C } = backbone(ch, P, j);
      const nextN = !ch.reverse ? (j + 1 < m ? backbone(ch, P, j + 1).N : ch.end ? P[3 * m] : null)
        : (j === 0 ? placedAt(res[i + 1]?.N) : backbone(ch, P, j - 1).N);
      const prevC = !ch.reverse ? (j === 0 ? placedAt(res[i - 1]?.C) : backbone(ch, P, j - 1).C)
        : (j + 1 < m ? backbone(ch, P, j + 1).C : null);
      const main = frameTransform(Q(r.CA!), Q(r.C!), Q(r.N!), CA, C, N);
      const oT = nextN && i + 1 < n && res[i + 1].N && bonded(i + 1)
        ? frameTransform(Q(r.C!), Q(r.CA!), Q(res[i + 1].N!), C, CA, nextN) : main;
      const hT = prevC && i > 0 && res[i - 1].C && bonded(i)
        ? frameTransform(Q(r.N!), Q(r.CA!), Q(res[i - 1].C!), N, CA, prevC) : main;
      return { main, oT, hT };
    });
  }

  function placeChain(ch: Chain, P: V[], placed: Map<Atom, V>, resT: (Rigid | null)[]) {
    const tf = chainTransforms(ch, P, placed);
    ch.idx.forEach((i, j) => {
      resT[i] = tf[j].main;
      for (const a of res[i].atoms) {
        if (a.model !== 0) continue;
        placed.set(a, apply(O_NAMES.has(a.name) ? tf[j].oT : N_H_NAMES.has(a.name) ? tf[j].hT : tf[j].main, Q(a)));
      }
    });
  }

  /** Re-shape a chain by phi/psi only (bond lengths/angles untouched).
   *   F (loops)       : close exactly onto the next body's N/CA/C (CCD, Canutescu & Dunbrack 2003)
   *   targets (2D)    : pull every CA towards its point on the planar guide curve
   *                     (multi-target CCD), then re-close; the closure always wins
   *  Restarts from random torsions (fixed seed: same input, same output), then Monte
   *  Carlo moves that remove heavy-atom clashes with everything already placed. */
  function solveChain(ch: Chain, targets: V[] | null, F: V[] | null, grid: Grid,
                      placed: Map<Atom, V>, seed: number) {
    const m = ch.idx.length;
    const tors = ch.tors;
    const rotate = (P: V[], tor: Tor, th: number) => {
      const O = P[tor.a], Rt = axisAngle(unit(sub(P[tor.b], O)), th);
      for (let j = tor.from; j < P.length; j++) P[j] = add(O, matVec(Rt, sub(P[j], O)));
    };
    // one CCD sweep towards weighted targets (point index -> target)
    const sweep = (P: V[], goals: { k: number; t: V; w: number }[], backwards: boolean) => {
      for (const tor of backwards ? [...tors].reverse() : tors) {
        const O = P[tor.a], u = unit(sub(P[tor.b], O));
        let B = 0, C = 0;
        for (const gl of goals) {
          if (gl.k < tor.from) continue;
          const mv = sub(P[gl.k], O);
          const r = sub(mv, scale(u, dot(mv, u)));
          const f = sub(gl.t, O);
          B += gl.w * dot(f, r); C += gl.w * dot(f, cross(u, r));
        }
        const th = Math.atan2(C, B);
        if (Math.abs(th) > 1e-7) rotate(P, tor, th);
      }
    };
    const endGoals = F ? F.map((t, k) => ({ k: 3 * m + k, t, w: 1 })) : [];
    const pathGoals = targets ? targets.map((t, j) => ({ k: 3 * j + 1, t, w: 1 })) : [];
    const guideGoals = [...pathGoals, ...endGoals.map((g2) => ({ ...g2, w: END_WEIGHT }))];
    const endErr = (P: V[]) => (F ? Math.sqrt(endGoals.reduce((s2, g2) => s2 + dist(P[g2.k], g2.t) ** 2, 0) / 3) : 0);
    const close = (P: V[], maxSweeps: number) => {
      let err = endErr(P);
      if (!F) return 0;
      for (let s = 0; s < maxSweeps && err > 0.05; s++) {
        sweep(P, endGoals, s % 2 === 1);
        const e2 = endErr(P);
        if (s > 40 && err - e2 < 1e-4) return e2;
        err = e2;
      }
      return err;
    };
    const guide = (P: V[], sweeps: number) => { for (let s = 0; s < sweeps; s++) sweep(P, guideGoals, s % 2 === 1); };
    const dev = (P: V[]) => (targets ? Math.sqrt(targets.reduce((s2, t, j) => s2 + dist(P[3 * j + 1], t) ** 2, 0) / m) : 0);
    const heavyOf = ch.idx.map((i) => res[i].atoms.filter((a) => a.model === 0 && !isHydrogen(a)));
    const clashes = (P: V[]) => {
      const tf = chainTransforms(ch, P, placed);
      const pts: { i: number; p: V }[] = [];
      ch.idx.forEach((i, j) => { for (const a of heavyOf[j]) pts.push({ i, p: apply(O_NAMES.has(a.name) ? tf[j].oT : tf[j].main, Q(a)) }); });
      let c = 0;
      for (const x of pts) c += grid.count(x.p, CLASH, (i) => Math.abs(i - x.i) >= 3);
      for (let u = 0; u < pts.length; u++) for (let v = u + 1; v < pts.length; v++)
        if (Math.abs(pts[v].i - pts[u].i) >= 3 && dist(pts[u].p, pts[v].p) < CLASH) c++;
      return c;
    };
    const run = (P: V[], sweeps: number, closeSweeps = 300) => {
      if (targets) guide(P, sweeps);
      let err = close(P, closeSweeps);
      if (targets) {
        for (let r = 0; r < 3; r++) { guide(P, Math.ceil(sweeps / 4)); err = close(P, Math.min(80, closeSweeps)); }
        if (err > 0.3) err = close(P, closeSweeps);          // closure has the last word
      }
      return err;
    };
    let rng = seed >>> 0;
    const random = () => {                                  // mulberry32: reproducible
      rng = (rng + 0x6d2b79f5) >>> 0;
      let t = rng;
      t = Math.imul(t ^ (t >>> 15), t | 1);
      t ^= t + Math.imul(t ^ (t >>> 7), t | 61);
      return ((t ^ (t >>> 14)) >>> 0) / 4294967296;
    };
    type Cand = { P: V[]; rmsd: number; clashes: number; dev: number };
    const score = (P: V[], err: number): Cand => ({ P, rmsd: err, clashes: err <= 0.3 ? clashes(P) : Infinity, dev: dev(P) });
    const better = (x: Cand, y: Cand) =>
      ((x.rmsd <= 0.3 ? 0 : 1) - (y.rmsd <= 0.3 ? 0 : 1) || x.clashes - y.clashes
        || (x.rmsd <= 0.3 ? x.dev - y.dev : x.rmsd - y.rmsd)) < 0;
    const done = (b: Cand) => b.rmsd <= 0.3 && b.clashes === 0 && (!targets || b.dev < GOOD_PATH_DEV);
    // fewer restarts for long chains (cost grows with length^2)
    const restarts = Math.max(6, Math.min(targets ? GUIDED_RESTARTS : CCD_RESTARTS,
                                          Math.round(CCD_RESTARTS * 30 / Math.max(30, m))));
    let best: Cand = { P: ch.start, rmsd: Infinity, clashes: Infinity, dev: Infinity };
    for (let attempt = 0; attempt < restarts && !done(best); attempt++) {
      const P = ch.start.map((p) => [...p] as V);
      if (attempt > 0) {
        const amp = Math.PI * Math.min(1, attempt / 8);
        for (const tor of tors) rotate(P, tor, (random() * 2 - 1) * amp);
      }
      const cand = score(P, run(P, 30));
      if (better(cand, best)) best = cand;
      if (attempt >= 15 && best.rmsd > 2) break;         // not closable while following the curve
    }
    // closure beats looks: if the guided search could not close the loop, search again
    // for closure alone (the loop then lies where the torsions allow)
    if (targets && F && best.rmsd > 0.3) {
      for (let attempt = 0; attempt < CCD_RESTARTS / 2 && best.rmsd > 0.3; attempt++) {
        const P = ch.start.map((p) => [...p] as V);
        if (attempt > 0) {
          const amp = Math.PI * Math.min(1, attempt / 8);
          for (const tor of tors) rotate(P, tor, (random() * 2 - 1) * amp);
        }
        const cand = score(P, close(P, 300));
        if (better(cand, best)) best = cand;
      }
    }
    const refineSteps = targets ? GUIDED_REFINE : CCD_REFINE;
    for (let step = 0; step < refineSteps && best.rmsd <= 0.3 && best.clashes > 0 && tors.length; step++) {
      const P = best.P.map((p) => [...p] as V);
      const moves = 1 + Math.floor(random() * 3);
      for (let k = 0; k < moves; k++) rotate(P, tors[Math.floor(random() * tors.length)], (random() * 2 - 1) * Math.PI / 3);
      const cand = score(P, run(P, 12, 100));
      if (better(cand, best)) best = cand;
    }
    return best;
  }

  /** Contour length a stretch of residues can cover along the drawing (Å). */
  function contour(idx: number[]) {
    return idx.reduce((s, i) => s + (ss(i) === 'H' ? STEP_HELIX : STEP_COIL), 0);
  }

  /** CA target of every residue: its share of the chain's contour along the curve. */
  function residueTargets(idx: number[], curve: Curve, scaleToCurve: boolean): V[] {
    const total = contour(idx) + (scaleToCurve ? STEP_COIL : 0);
    const k = scaleToCurve ? curve.length / Math.max(1e-6, total) : 1;
    let s = 0;
    return idx.map((i) => { s += ss(i) === 'H' ? STEP_HELIX : STEP_COIL; return curveAt(curve, s * k); });
  }

  // ---- 6. ligands follow their nearest residue; waters / dummies / other chains dropped
  const chainAtoms = res.flatMap((r, i) => r.atoms.filter((a) => a.model === 0).map((a) => ({ a, i })));
  const hetero = new Map<string, Atom[]>();
  const droppedCount = new Map<string, number>();
  const drop = (what: string) => droppedCount.set(what, (droppedCount.get(what) ?? 0) + 1);
  const chainSet = new Set(res);
  for (const r of allRes) {
    if (chainSet.has(r)) continue;
    if (WATER.has(r.resName)) { drop('water'); continue; }
    if (DUMMY.has(r.resName)) { drop('membrane dummy atoms (DUM)'); continue; }
    if (isAminoAcid(r)) { drop(`chain ${r.chain} residues`); continue; }
    hetero.set(resLabel(r) + r.resName, r.atoms.filter((a) => a.model === 0));
  }
  const keptLigands: string[] = [];
  for (const [label, latoms] of hetero) {
    let best = { d: 4.5, i: -1 };
    for (const la of latoms) for (const { a, i } of chainAtoms) {
      const d = dist(Q(la), Q(a));
      if (d < best.d) best = { d, i };
    }
    if (best.i < 0 || !resT[best.i]) { drop('ligands not in contact with the chain'); continue; }
    for (const la of latoms) placed.set(la, apply(resT[best.i]!, Q(la)));
    keptLigands.push(`${label} (with ${resLabel(res[best.i])})`);
  }

  // ---- 7. checks: peptide bonds and clashes in the layout
  let maxBond = 0, broken = 0;
  for (let i = 1; i < n; i++) {
    const c = res[i - 1].C, nn = res[i].N;
    if (!c || !nn || dist(c.p, nn.p) >= 2.0) continue;
    const d = dist(placed.get(c)!, placed.get(nn)!);
    maxBond = Math.max(maxBond, d);
    if (d > 1.5) broken++;
  }
  // heavy-atom overlaps between residues >= 3 apart, in the layout but not in the input
  const heavy = res.flatMap((r, i) => r.atoms.filter((a) => a.model === 0 && !isHydrogen(a)).map((a) => ({ a, i })));
  const gridIn = new Grid(3.0), gridOut = new Grid(3.0);
  heavy.forEach(({ a, i }) => { gridIn.add(Q(a), i); gridOut.add(placed.get(a)!, i); });
  let clashes = 0;
  for (const { a, i } of heavy) {
    const far = (j: number) => j - i >= 3;
    clashes += Math.max(0, gridOut.count(placed.get(a)!, CLASH, far) - gridIn.count(Q(a), CLASH, far));
  }
  // disulfides are covalent bonds too
  const sg = heavy.filter(({ a }) => a.name === 'SG' && ONE[a.resName] === 'C');
  let brokenSS = 0;
  for (let u = 0; u < sg.length; u++) for (let v = u + 1; v < sg.length; v++) {
    if (dist(Q(sg[u].a), Q(sg[v].a)) < 2.5 && dist(placed.get(sg[u].a)!, placed.get(sg[v].a)!) > 2.5) {
      brokenSS++;
      warnings.push(`disulfide ${resLabel(res[sg[u].i])}-${resLabel(res[sg[v].i])} is broken in the layout`);
    }
  }
  if (broken) warnings.push(`${broken} peptide bond(s) stretched beyond 1.5 Å: schematic junction(s) - draw them dashed`);
  if (clashes) warnings.push(`${clashes} heavy-atom overlap(s) < ${CLASH} Å created by the layout`);

  // ---- 8. write
  const header = lines.filter((l, li) => li < (atoms[0]?.li ?? lines.length)
    && !/^(END|MASTER|CONECT|ANISOU|MODEL|ENDMDL)/.test(l) && l.trim() !== '');
  const body: string[] = [];
  for (const r of res) for (const a of r.atoms) if (a.model === 0) body.push(writeCoords(lines[a.li], placed.get(a)!));
  body.push(`TER   ${' '.repeat(5)}      ${res[n - 1].resName.padStart(3)} ${chainId}${String(res[n - 1].resSeq).padStart(4)}${res[n - 1].iCode || ' '}`);
  for (const latoms of hetero.values()) for (const a of latoms) if (placed.has(a)) body.push(writeCoords(lines[a.li], placed.get(a)!));
  const crossings = bodies.map((b) => ({ crossing: b.id, start: resLabel(res[b.s]), end: resLabel(res[b.e]),
    tiltRemovedDeg: groups.find((g) => g.members.includes(b.id))!.members.length > 1 ? 0 : Math.round(tiltOf.get(b.id)! * 10) / 10 }));
  const remarks = [
    mode === 'exploded'
      ? 'SPREAD MODE: EXPLODED - bodies pushed apart in the membrane plane.'
      : 'SPREAD MODE: UNROLLED - SCHEMATIC ROW LAYOUT, NOT FOR INTERACTION ANALYSIS.',
    'Y = membrane normal (extracellular up), y = 0 at the bilayer mid-plane.',
    'TM crossings, tails and loop helices/strands are moved as rigid bodies;',
    'loops are re-closed by changing coil phi/psi only (CCD): bond lengths and',
    'angles are those of the input model.',
    ...(mode === 'unrolled' ? ['2D layout: loops and tails follow planar (z = 0) guide curves within',
                               `${loopHeight} A of the membrane surface.`] : []),
    ...loops.map((l) => `LOOP ${l.from}-${l.to}: ${l.residues} res, closure ${l.closureRmsd ?? 'n/a'} A` +
      (l.schematic ? ` SCHEMATIC (C-N ${l.stretchedBond ?? '?'} A)` : '')),
  ];
  const out = [...header, ...body, 'END'];
  insertRemarks(out, remarks);
  return {
    pdb: out.join('\n'),
    report: {
      ...reportBase, crossings,
      keptTogether: groups.filter((g) => g.members.length > 1).map((g) => g.members),
      loops, tails, maxPeptideBond: Math.round(maxBond * 100) / 100, brokenPeptideBonds: broken, clashes,
      brokenDisulfides: brokenSS,
      keptLigands, dropped: [...droppedCount.entries()].map(([k, v]) => `${v} ${k}`),
    },
  };
}

/* ------------------------------------------------------------------ *
 * Planar guide curves for the 2D layout
 * ------------------------------------------------------------------ */

interface Curve { pts: V[]; cum: number[]; length: number }

function makeCurve(pts: V[]): Curve {
  const cum = [0];
  for (let i = 1; i < pts.length; i++) cum.push(cum[i - 1] + dist(pts[i - 1], pts[i]));
  return { pts, cum, length: cum[cum.length - 1] };
}

function curveAt(c: Curve, s: number): V {
  const x = Math.max(0, Math.min(c.length, s));
  let lo = 0, hi = c.cum.length - 1;
  while (hi - lo > 1) { const mid = (lo + hi) >> 1; if (c.cum[mid] <= x) lo = mid; else hi = mid; }
  const t = (x - c.cum[lo]) / (c.cum[hi] - c.cum[lo] || 1);
  return add(c.pts[lo], scale(sub(c.pts[hi], c.pts[lo]), t));
}

/**
 * Curve in the drawing plane (z = 0) for a loop from A to B, or a tail from A going
 * sideways (dir -1 left / +1 right). It leaves the membrane vertically, runs horizontally
 * in the band yLo..yHi (|y|, on the loop's side of the membrane) and comes back down to
 * B. Its length matches the chain's contour: short loops arch low, longer ones arch
 * higher, and what does not fit under yHi meanders sideways (a horizontal sine wave)
 * instead of sticking out of the drawing.
 */
function guideCurve(A: V, B: V | null, side: number, length: number, dir: number,
                    yLo: number, yHi: number, maxWidth: number): Curve {
  const s = side >= 0 ? 1 : -1;
  const legs = (h: number) => Math.abs(s * h - A[1]) + (B ? Math.abs(s * h - B[1]) : 0);
  const x0 = A[0];
  const x1 = B ? B[0] : x0 + dir * Math.min(maxWidth, Math.max(6, length - legs(yLo)));
  const D = Math.abs(x1 - x0), sg = Math.sign(x1 - x0) || 1;
  const flatLen = (h: number) => legs(h) + D;
  let h = yLo, m = 0, a = 0;
  if (flatLen(yLo) < length) {
    if (flatLen(yHi) >= length) {                         // a flat top at the right height
      let lo = yLo, hi = yHi;
      for (let it = 0; it < 40; it++) { const mid = (lo + hi) / 2; if (flatLen(mid) < length) lo = mid; else hi = mid; }
      h = (lo + hi) / 2;
    } else {                                              // meander inside the band
      h = (yLo + yHi) / 2;
      const amax = (yHi - yLo) / 2, need = length - legs(h);
      const waveLen = (mm: number, aa: number) => {
        let L = 0, px = 0, py = 0;
        for (let t = 1; t <= 200; t++) {
          const x = (D * t) / 200, y = aa * Math.sin((Math.PI * mm * t) / 200);
          L += Math.hypot(x - px, y - py); px = x; py = y;
        }
        return L;
      };
      const mMax = Math.max(1, Math.floor(D / 4));      // half-waves at least 4 Å wide
      for (let mm = 1; mm <= mMax && !m; mm++) {
        if (waveLen(mm, amax) < need) continue;
        m = mm;
        let lo = 0, hi = amax;
        for (let it = 0; it < 40; it++) { const mid = (lo + hi) / 2; if (waveLen(mm, mid) < need) lo = mid; else hi = mid; }
        a = (lo + hi) / 2;
      }
      if (!m) { m = mMax; a = amax; }                     // longer than the band allows
    }
  }
  const pts: V[] = [A, [x0, s * h, 0]];
  const steps = Math.max(8, 24 * Math.max(1, m));
  for (let t = 1; t <= steps; t++) {
    const u = t / steps;
    pts.push([x0 + sg * D * u, s * (h + a * Math.sin(Math.PI * m * u)), 0]);
  }
  if (B) pts.push(B);
  return makeCurve(pts);
}

/* ------------------------------------------------------------------ *
 * Spatial hash for clash tests
 * ------------------------------------------------------------------ */

const isHydrogen = (a: Atom) => /^(\d?H|D)/.test(a.name) && a.resName !== 'HG';

class Grid {
  private cells = new Map<string, { p: V; i: number }[]>();
  constructor(private size: number) {}
  private key(x: number, y: number, z: number) { return `${x}|${y}|${z}`; }
  add(p: V, i: number) {
    const k = this.key(Math.floor(p[0] / this.size), Math.floor(p[1] / this.size), Math.floor(p[2] / this.size));
    let c = this.cells.get(k);
    if (!c) { c = []; this.cells.set(k, c); }
    c.push({ p, i });
  }
  /** number of stored points closer than r to p whose residue index passes `keep` */
  count(p: V, r: number, keep: (i: number) => boolean) {
    const [cx, cy, cz] = [Math.floor(p[0] / this.size), Math.floor(p[1] / this.size), Math.floor(p[2] / this.size)];
    let c = 0;
    for (let x = cx - 1; x <= cx + 1; x++) for (let y = cy - 1; y <= cy + 1; y++) for (let z = cz - 1; z <= cz + 1; z++) {
      for (const q of this.cells.get(this.key(x, y, z)) ?? []) if (keep(q.i) && dist(p, q.p) < r) c++;
    }
    return c;
  }
}

/* ------------------------------------------------------------------ *
 * Output helpers
 * ------------------------------------------------------------------ */

/** Coordinates in PDB columns 31-54 (8.3f); falls back to fewer decimals rather than
 *  shifting the columns when a value does not fit. */
function writeCoords(line: string, p: V): string {
  const f = (v: number) => {
    for (const d of [3, 2, 1, 0]) { const s = v.toFixed(d); if (s.length <= 8) return s.padStart(8); }
    return (v < 0 ? '-9999999' : '99999999');
  };
  return line.substring(0, 30).padEnd(30) + f(p[0]) + f(p[1]) + f(p[2]) + line.substring(54);
}

function insertRemarks(out: string[], texts: string[]) {
  const at = Math.max(0, out.findIndex((l) => /^(ATOM  |HETATM|MODEL)/.test(l)));
  out.splice(at, 0, ...texts.map((t) => `REMARK 999 ${t}`.slice(0, 80)));
}
