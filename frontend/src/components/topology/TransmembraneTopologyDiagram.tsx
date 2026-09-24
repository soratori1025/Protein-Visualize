import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import type { Chain } from '../../types/protein';
import type {
  CalculatedTopologyData,
  SecondaryStructureResult,
  UniProtTopologyData,
  UniProtTopologyRegion,
} from '../../types/secondaryStructure';
import { compareResidues, isCalculatedTopology } from '../../types/secondaryStructure';
import { exportSvgAsImage } from '../structure/exportDiagram';
import { API_URL, runSecondaryStructure } from '../../services/api';
import './TransmembraneTopologyDiagram.css';

type FigureTheme = 'publication' | 'lab';
type TopologySource = 'uniprot' | 'calculated';
type MembraneSide = 'out' | 'in'; // out = extracellular/lumenal, in = cytoplasmic

/* ------------------------------------------------------------------ *
 * Color utilities
 * ------------------------------------------------------------------ */

function hexToRgb(hex: string): [number, number, number] {
  const h = hex.replace('#', '').trim();
  const full = h.length === 3 ? h.split('').map((c) => c + c).join('') : h;
  return [
    parseInt(full.slice(0, 2), 16) || 0,
    parseInt(full.slice(2, 4), 16) || 0,
    parseInt(full.slice(4, 6), 16) || 0,
  ];
}

function rgbToHex(r: number, g: number, b: number): string {
  const c = (v: number) => Math.max(0, Math.min(255, Math.round(v))).toString(16).padStart(2, '0');
  return `#${c(r)}${c(g)}${c(b)}`;
}

/** Linear blend between two hex colors. amount = 0 → a, 1 → b. */
function mixHex(a: string, b: string, amount: number): string {
  const [r1, g1, b1] = hexToRgb(a);
  const [r2, g2, b2] = hexToRgb(b);
  const t = Math.max(0, Math.min(1, amount));
  return rgbToHex(r1 + (r2 - r1) * t, g1 + (g2 - g1) * t, b1 + (b2 - b1) * t);
}

const lighten = (hex: string, amount: number) => mixHex(hex, '#ffffff', amount);
const darken = (hex: string, amount: number) => mixHex(hex, '#000000', amount);

/** WCAG relative luminance — decides whether a label should be white or near-black. */
function relativeLuminance(hex: string): number {
  const [r, g, b] = hexToRgb(hex).map((v) => {
    const s = v / 255;
    return s <= 0.03928 ? s / 12.92 : Math.pow((s + 0.055) / 1.055, 2.4);
  });
  return 0.2126 * r + 0.7152 * g + 0.0722 * b;
}

/**
 * Text color is decided purely by the background it sits on.
 * (The old version forced dark text in publication mode, which made labels on
 * dark helices such as TM11/TM12 almost unreadable.)
 */
function getContrastTextColor(bgHex: string): string {
  return relativeLuminance(bgHex) > 0.45 ? '#0f172a' : '#ffffff';
}

/**
 * Resample a base palette to exactly `count` colors by interpolating between
 * neighbouring stops. Works for 3 TMs or 24 TMs — no modulo wrap-around, so two
 * helices never share the same color.
 */
function buildPalette(base: string[], count: number): string[] {
  if (count <= 0) return [];
  if (count === 1) return [base[0]];
  return Array.from({ length: count }, (_, i) => {
    const t = (i / (count - 1)) * (base.length - 1);
    const lo = Math.floor(t);
    const hi = Math.min(base.length - 1, lo + 1);
    return mixHex(base[lo], base[hi], t - lo);
  });
}

/* ------------------------------------------------------------------ *
 * Topology helpers
 * ------------------------------------------------------------------ */

// Same vocabulary as the backend (app/services/topology/labels.py). OUT is checked
// first so "Perinuclear space" / "Mitochondrial intermembrane" are not read as inside.
const OUT_WORDS = [
  'extracellular', 'lumenal', 'luminal', 'periplasm', 'exoplasmic', 'vesicular',
  'intermembrane', 'perinuclear', 'vacuolar', 'peroxisomal', 'thylakoid', 'virion surface',
];
const IN_WORDS = ['cytoplasm', 'cytosol', 'intracellular', 'mitochondrial matrix', 'matrix', 'stromal', 'intravirion'];

/** Map a UniProt topological-domain description to a membrane side, or null if unknown. */
function descriptionToSide(description?: string): MembraneSide | null {
  const d = (description ?? '').toLowerCase();
  if (!d) return null;
  if (OUT_WORDS.some((w) => d.includes(w))) return 'out';
  if (IN_WORDS.some((w) => d.includes(w)) || d.trim() === 'nuclear') return 'in';
  return null;
}

/** Side of a region: the backend's explicit `side` first, the description otherwise. */
function regionSide(region: Pick<UniProtTopologyRegion, 'side' | 'description'>): MembraneSide | null {
  if (region.side === 'Cytoplasmic') return 'in';
  if (region.side === 'Extracellular') return 'out';
  return descriptionToSide(region.description);
}

const opposite = (s: MembraneSide): MembraneSide => (s === 'out' ? 'in' : 'out');

/**
 * Side on which the FIRST crossing starts. Taken from the first topological domain
 * whose side is known, flipped once per crossing that lies before it. (The old
 * version used that domain's side directly, which draws the whole map upside down
 * whenever the first annotated domain comes after TM1 — e.g. a chain that starts
 * inside the membrane, or a UniProt entry without an N-terminal domain.)
 */
function getNTermSide(data: UniProtTopologyData | null, crossingStarts: number[]): MembraneSide {
  const domains = (data?.regions ?? [])
    .filter((r) => r.type === 'Topological domain' && regionSide(r) !== null)
    .sort((a, b) => a.start - b.start);
  const first = domains[0];
  if (!first) return 'in';
  const side = regionSide(first) as MembraneSide;
  const crossingsBefore = crossingStarts.filter((start) => start < first.start).length;
  return crossingsBefore % 2 === 0 ? side : opposite(side);
}

/** Keywords UniProt uses for helices that are broken in the middle of the bilayer. */
const DISCONTINUOUS_WORDS = [
  'discontinuous',
  'broken',
  'kink',
  'unwound',
  'interrupted',
  'part 1',
  'part 2',
  'first part',
  'second part',
];

const LONG_TM_RESIDUES = 30; // a single straight TM helix is ~18–25 residues

interface RawTM {
  start: number;
  end: number;
  name?: string;
  description?: string;
  /** 'Helix' | 'Strand' | 'Loop' from the calculated endpoint. */
  ss?: string | null;
}

function mentionsDiscontinuity(tm: RawTM): boolean {
  const text = `${tm.name ?? ''} ${tm.description ?? ''}`.toLowerCase();
  return DISCONTINUOUS_WORDS.some((w) => text.includes(w));
}

/**
 * Decide whether a transmembrane segment is a β-strand rather than an α-helix,
 * from its description. UniProt writes beta-barrel TM segments as "Beta stranded"
 * (occasionally "Beta-stranded"); the structure-based predictor writes
 * "Transmembrane Beta Strand". Handles all of those plus "beta barrel"/"sheet",
 * and never fires on "Helical".
 */
function isBetaStrandDescription(description?: string): boolean {
  const d = (description ?? '').toLowerCase();
  if (!d) return false;
  const hasBeta = d.includes('beta') || d.includes('β');
  return hasBeta && (d.includes('strand') || d.includes('barrel') || d.includes('sheet'));
}

/** β segment: explicit `ss` from the calculated endpoint, else the description. */
function isBetaSegment(tm: { ss?: string | null; description?: string }): boolean {
  return tm.ss === 'Strand' || isBetaStrandDescription(tm.description);
}

/**
 * Two consecutive annotated segments that are really the two halves of one
 * broken helix: nearly touching, each too short to span the bilayer on its own,
 * or explicitly named "…a"/"…b" or "part 1"/"part 2".
 */
function looksLikeHalfPair(prev: RawTM, curr: RawTM, trustSegments = false): boolean {
  const gap = curr.start - prev.end - 1;
  if (gap < 0 || gap > 7) return false;

  const prevText = `${prev.name ?? ''} ${prev.description ?? ''}`.toLowerCase();
  const currText = `${curr.name ?? ''} ${curr.description ?? ''}`.toLowerCase();

  const prevNum = prevText.match(/(\d+)\s*a\b/);
  const currNum = currText.match(/(\d+)\s*b\b/);
  if (prevNum && currNum && prevNum[1] === currNum[1]) return true;
  if (prevText.includes('part 1') && currText.includes('part 2')) return true;

  // Calculated topology: the backend already decided what is one crossing (it fuses
  // kinked / broken helices and splits antiparallel hairpins using the 3D chain
  // direction). Guessing again from lengths would fuse two real crossings.
  if (trustSegments) return false;

  // β-strands are legitimately short and pack close together in a barrel — never
  // fuse two of them into "the two halves of one broken helix" (that heuristic is
  // only meant for a single α-helix that UniProt split across the bilayer).
  if (isBetaSegment(prev) || isBetaSegment(curr)) return false;

  const prevLen = prev.end - prev.start + 1;
  const currLen = curr.end - curr.start + 1;
  return prevLen <= 16 && currLen <= 16;
}

/* ------------------------------------------------------------------ *
 * Types
 * ------------------------------------------------------------------ */

interface Props {
  chain: Chain | undefined;
  secondaryResult?: SecondaryStructureResult | null;
  selectedResidue?: number | null;
  onSelectResidue?: (residueNumber: number) => void;
  uniprotId?: string | null;
  filename?: string | null;
  tmAlgorithm?: string;
  onTmAlgorithmChange?: (algo: string) => void;
  topologySource?: TopologySource;
  onTopologySourceChange?: (source: TopologySource) => void;
  onTopologyDataChange?: (data: UniProtTopologyData | null) => void;
  triggerTmRecalc?: number;
}

export interface TMHelix {
  id: string;
  helixNumber: number;
  subLabel: string; // "1a", "1b", "2", "6a"…
  startRes: number;
  endRes: number;
  length: number;
  color: string;
  isSplit?: boolean;
  partIndex?: number; // 0 = first half (entry side), 1 = second half
  description?: string;
  /** β-strand crossing (drawn as an arrow) rather than an α-helix. */
  isBeta?: boolean;
  /** Side of the membrane where this segment's N-terminal end sits. */
  entrySide: MembraneSide;
  /** Side of the membrane where this segment's C-terminal end sits. */
  exitSide: MembraneSide;
  /** Column index — both halves of a broken helix share one column. */
  column: number;
}

export interface ExtraFeature {
  label: string;
  startRes: number;
  endRes: number;
  type: string;
  offsetFactor?: number;
}

export interface TMLoop {
  id: string;
  type: 'EL' | 'IL';
  label: string;
  startRes: number;
  endRes: number;
  length: number;
  prevHelixId: string;
  nextHelixId: string;
  hasExtraFeature?: boolean;
  extraFeatures?: ExtraFeature[];
  domainName?: string;
}

export interface CustomResidueColorRule {
  id: string;
  label: string;
  startRes: number;
  endRes: number;
  color: string;
}

/* ------------------------------------------------------------------ *
 * Palettes
 * ------------------------------------------------------------------ */

const PALETTES: Record<string, { label: string; colors: string[] }> = {
  PAPER_DEFAULT: {
    label: 'Figure paper default',
    colors: ['#ef4444', '#f97316', '#eab308', '#22c55e', '#06b6d4', '#2563eb', '#7c3aed', '#ec4899', '#64748b', '#1f2937'],
  },
  RAINBOW: {
    label: 'Rainbow spectrum',
    colors: ['#ef4444', '#f97316', '#eab308', '#10b981', '#06b6d4', '#3b82f6', '#8b5cf6', '#d946ef'],
  },
  HYDROPHOBIC: {
    label: 'Hydrophobicity ramp',
    colors: ['#1d4ed8', '#60a5fa', '#e2e8f0', '#fdba74', '#f97316', '#b45309'],
  },
  NEON: {
    label: 'Vivid neon',
    colors: ['#ff0055', '#ff6600', '#ffcc00', '#00ff66', '#00ffff', '#0099ff', '#9900ff', '#ff00ff'],
  },
  PASTEL: {
    label: 'Pastel bio',
    colors: ['#fca5a5', '#fdba74', '#fde047', '#86efac', '#67e8f9', '#93c5fd', '#c084fc', '#f472b6'],
  },
};

const DEFAULT_PRESETS = [
  { id: 'P31645', label: 'P31645 (hSERT)' },
  { id: 'P23975', label: 'P23975 (hNET)' },
  { id: 'P08183', label: 'P08183 (P-gp)' },
  { id: 'P00533', label: 'P00533 (EGFR)' },
];

/* ------------------------------------------------------------------ *
 * Re-entrant / extra features inside loops and terminals
 * ------------------------------------------------------------------ */

export function getExtraFeatures(
  lStart: number,
  lEnd: number,
  topologyData: UniProtTopologyData | null,
  secondaryResult: SecondaryStructureResult | null | undefined,
  chainId: string | undefined
): ExtraFeature[] | undefined {
  return undefined; // Tạm thời ẩn các cấu trúc phụ ngoài màng
  /*
  if (lEnd < lStart) return undefined;

  const found: ExtraFeature[] = [];

  // 1. UniProt features
  for (const region of topologyData?.regions ?? []) {
    if (
      ['Intramembrane', 'Signal', 'Domain', 'Region', 'Motif', 'Propeptide', 'Repeat'].includes(region.type) &&
      region.start >= lStart &&
      region.end <= lEnd
    ) {
      found.push({
        startRes: region.start,
        endRes: region.end,
        type: region.type,
        label: region.name || region.description || region.type,
      });
    }
  }

  // 1b. Extramembrane secondary structure straight from the predictor's region
  //     schema. The structure-based predictor now splits each loop into its own
  //     Helix / Strand / Coil regions and tags every region with `side` + `ss`,
  //     so a helix or strand sitting OUTSIDE the bilayer (e.g. EL2, EL3a/3b) is a
  //     "Topological domain" region carrying ss:"Helix"|"Strand". Draw those
  //     directly — this is what makes loop structure appear for Calculated
  //     topology without needing a separate DSSP/STRIDE overlay run.
  for (const region of topologyData?.regions ?? []) {
    const ss = (region as { ss?: string }).ss;
    if (
      region.type === 'Topological domain' &&
      (ss === 'Helix' || ss === 'Strand') &&
      region.start >= lStart &&
      region.end <= lEnd
    ) {
      found.push({
        startRes: region.start,
        endRes: region.end,
        type: ss, // 'Helix' | 'Strand'
        label: ss === 'Helix' ? 'α-helix' : 'β-strand',
      });
    }
  }

  // 2. Helical runs from assigned secondary structure (fallback when the region
  //    schema carries no `ss`, e.g. UniProt data or an older backend). The overlap
  //    guard below skips anything already captured from the region schema above.
  if (secondaryResult?.residues && chainId) {
    const loopResidues = secondaryResult.residues
      .filter((r) => r.chain_id === chainId && r.residue_number >= lStart && r.residue_number <= lEnd)
      .sort(compareResidues);

    let run: { startRes: number; endRes: number } | null = null;
    const flush = () => {
      // Only keep helices >= 6 residues, and avoid overlapping heavily with UniProt features
      if (run && run.endRes - run.startRes + 1 >= 6) {
        const r = run;
        const overlap = found.some(f => Math.max(0, Math.min(f.endRes, r.endRes) - Math.max(f.startRes, r.startRes)) > 4);
        if (!overlap) {
          found.push({ ...r, type: 'Helix', label: 'α-helix' });
        }
      }
      run = null;
    };
    for (const r of loopResidues) {
      if (['H', 'G', 'I'].includes(r.code.toUpperCase())) {
        if (!run) run = { startRes: r.residue_number, endRes: r.residue_number };
        else run.endRes = r.residue_number;
      } else {
        flush();
      }
    }
    flush();
  }

  if (found.length === 0) return undefined;

  // Keep the longest 3 to avoid cluttering
  const kept = found
    .sort((a, b) => b.endRes - b.startRes - (a.endRes - a.startRes))
    .slice(0, 3)
    .sort((a, b) => a.startRes - b.startRes);

  return kept.map((sh, idx) => ({
    ...sh,
    offsetFactor: kept.length > 1 ? idx - (kept.length - 1) / 2 : 0,
  }));
  */
}

/* ------------------------------------------------------------------ *
 * Layout model (pure: exported so it can be unit-tested without React)
 * ------------------------------------------------------------------ */

export interface TopologyModelInput {
  chain: Chain | undefined;
  /** Topology being drawn (UniProt or calculated); null = derive from SS. */
  activeTopologyData: UniProtTopologyData | null;
  /** SS of the uploaded structure (used only for the no-topology fallback). */
  secondaryResult?: SecondaryStructureResult | null;
  /** SS in the SAME numbering as the topology, or null (UniProt numbering). */
  structureSS?: SecondaryStructureResult | null;
  topologySource: TopologySource;
  selectedPaletteKey?: string;
  customHelixColors?: Record<string, string>;
}

export function buildTopologyModel({
  chain,
  activeTopologyData,
  secondaryResult = null,
  structureSS = null,
  topologySource,
  selectedPaletteKey = 'PAPER_DEFAULT',
  customHelixColors = {},
}: TopologyModelInput): { helices: TMHelix[]; loops: TMLoop[]; helixCount: number } {
  const empty = { helices: [] as TMHelix[], loops: [] as TMLoop[], helixCount: 0 };
  if (!chain) return empty;

  /* --- Step 1: collect raw transmembrane segments --- */
  let rawTMs: RawTM[] = [];
  let domainRegions: { start: number; end: number; description: string }[] = [];
  let usingAnnotation = false;

  // A loaded topology with no Transmembrane region means "no TM segments" (e.g. a
  // soluble protein: the backend returns regions: []). Only when there is no
  // topology at all are candidate helices derived from the SS assignment.
  if (activeTopologyData) {
    usingAnnotation = true;
    rawTMs = (activeTopologyData.regions ?? [])
      .filter((r) => r.type === 'Transmembrane')
      .map((r) => ({ start: r.start, end: r.end, name: r.name, description: r.description, ss: r.ss }))
      .sort((a, b) => a.start - b.start);
    domainRegions = activeTopologyData.regions
      .filter((r) => r.type === 'Topological domain')
      .map((r) => ({ start: r.start, end: r.end, description: r.description ?? '' }));
  } else {
    // Fallback: derive candidate TM segments from the assigned secondary structure.
    const ssResidues = (secondaryResult?.residues ?? [])
      .filter((r) => r.chain_id === chain.id)
      .sort(compareResidues);

    let run: RawTM | null = null;
    const flush = () => {
      if (run && run.end - run.start + 1 >= 15) rawTMs.push(run);
      run = null;
    };
    for (const r of ssResidues) {
      if (['H', 'G', 'I'].includes(r.code.toUpperCase())) {
        // a gap in the numbering (unresolved residues) ends the helix
        if (run && r.residue_number > run.end + 1) flush();
        if (!run) run = { start: r.residue_number, end: r.residue_number };
        else run.end = r.residue_number;
      } else {
        flush();
      }
    }
    flush();
  }

  if (rawTMs.length === 0) return empty;

  /* --- Step 2: group segments into helices (a broken helix = one group of 2) --- */
  const trustSegments = topologySource === 'calculated';
  const groups: RawTM[][] = [];
  for (const tm of rawTMs) {
    const last = groups[groups.length - 1];
    if (last && last.length === 1 && looksLikeHalfPair(last[0], tm, trustSegments)) last.push(tm);
    else groups.push([tm]);
  }

  /* --- Step 3: palette sized to the actual number of helices --- */
  const base = PALETTES[selectedPaletteKey]?.colors ?? PALETTES.PAPER_DEFAULT.colors;
  const palette = buildPalette(base, groups.length);
  const colorFor = (subLabel: string, hNum: number) =>
    customHelixColors[subLabel] ?? customHelixColors[`TM${hNum}`] ?? palette[hNum - 1] ?? base[0];

  /* --- Step 4: which side does each helix start on? --- */
  const nTermSide = usingAnnotation
    ? getNTermSide(activeTopologyData, groups.map((g) => g[0].start))
    : 'in';

  // Full transmembrane helices MUST strictly alternate sides. 
  // Relying on incomplete 'Topological domain' annotations can result in consecutive ELs.
  const entrySides: MembraneSide[] = groups.map((g, i) => {
    return i % 2 === 0 ? nTermSide : opposite(nTermSide);
  });

  /* --- Step 5: expand groups into drawable segments --- */
  const tmHelices: TMHelix[] = [];

  groups.forEach((group, gi) => {
    const hNum = gi + 1;
    const entrySide = entrySides[gi];
    const exitSide = opposite(entrySide);
    const first = group[0];
    const last = group[group.length - 1];
    const spanLen = last.end - first.start + 1;

    const preSplit = group.length === 2;
    // Auto-splitting a long segment into "a/b" halves is a UniProt heuristic. A
    // calculated crossing is already one element (snapped helices of 30+ residues
    // are normal there), so it is never cut at its midpoint.
    const shouldSplit =
      preSplit || mentionsDiscontinuity(first) || (!trustSegments && spanLen >= LONG_TM_RESIDUES);
    const isBeta = isBetaSegment(first);

    if (!shouldSplit) {
      const subLabel = `${hNum}`;
      tmHelices.push({
        id: `tm-${hNum}`,
        helixNumber: hNum,
        subLabel,
        startRes: first.start,
        endRes: last.end,
        length: spanLen,
        color: colorFor(subLabel, hNum),
        isSplit: false,
        description: first.description,
        isBeta,
        entrySide,
        exitSide,
        column: gi,
      });
      return;
    }

    const halves = preSplit
      ? [
          { start: first.start, end: first.end, description: first.description },
          { start: last.start, end: last.end, description: last.description },
        ]
      : (() => {
          const mid = Math.floor((first.start + last.end) / 2);
          return [
            { start: first.start, end: mid, description: first.description },
            { start: mid + 1, end: last.end, description: first.description },
          ];
        })();

    halves.forEach((half, idx) => {
      const subLabel = `${hNum}${idx === 0 ? 'a' : 'b'}`;
      tmHelices.push({
        id: `tm-${subLabel}`,
        helixNumber: hNum,
        subLabel,
        startRes: half.start,
        endRes: half.end,
        length: half.end - half.start + 1,
        color: colorFor(subLabel, hNum),
        isSplit: true,
        partIndex: idx,
        description: half.description,
        isBeta,
        // Part a runs from the entry side inwards; part b carries on to the exit side.
        entrySide: idx === 0 ? entrySide : entrySide,
        exitSide: idx === 0 ? entrySide : exitSide,
        column: gi,
      });
    });
  });

  /* --- Step 6: loops between consecutive helices --- */
  const tmLoops: TMLoop[] = [];
  let elCount = 1;
  let ilCount = 1;

  for (let i = 0; i < tmHelices.length - 1; i++) {
    const hCurr = tmHelices[i];
    const hNext = tmHelices[i + 1];
    // Two halves of the same helix are joined inside the bilayer, not by a loop.
    if (hCurr.column === hNext.column) continue;

    const lStart = hCurr.endRes + 1;
    const lEnd = hNext.startRes - 1;
    const lLen = Math.max(0, lEnd - lStart + 1);
    const side = hCurr.exitSide;
    const label = side === 'out' ? `EL${elCount++}` : `IL${ilCount++}`;

    const midRes = Math.floor((lStart + lEnd) / 2);
    const matchedDomain = domainRegions.find((d) => midRes >= d.start && midRes <= d.end);

    const extraFeatures = getExtraFeatures(
      lStart,
      lEnd,
      usingAnnotation ? activeTopologyData : null,
      structureSS,
      chain.id
    );

    tmLoops.push({
      id: `loop-${label}`,
      type: side === 'out' ? 'EL' : 'IL',
      label,
      startRes: lStart,
      endRes: lEnd,
      length: lLen,
      prevHelixId: hCurr.id,
      nextHelixId: hNext.id,
      hasExtraFeature: extraFeatures && extraFeatures.length > 0,
      extraFeatures,
      domainName: matchedDomain?.description,
    });
  }

  return { helices: tmHelices, loops: tmLoops, helixCount: groups.length };
}

/* ------------------------------------------------------------------ *
 * Component
 * ------------------------------------------------------------------ */

export function TransmembraneTopologyDiagram({
  chain,
  secondaryResult,
  filename,
  onSelectResidue,
  selectedResidue,
  uniprotId,
  tmAlgorithm: propsTmAlgorithm,
  onTmAlgorithmChange,
  topologySource: propsTopologySource,
  onTopologySourceChange,
  onTopologyDataChange,
  triggerTmRecalc,
}: Props) {
  const [uniprotIdInput, setUniprotIdInput] = useState<string>('P31645');
  const [uniprotData, setUniprotData] = useState<UniProtTopologyData | null>(null);
  const [loadingUniProt, setLoadingUniProt] = useState<boolean>(false);
  const [uniprotError, setUniprotError] = useState<string | null>(null);

  const [internalTopologySource, setInternalTopologySource] = useState<TopologySource>('uniprot');
  const topologySource = propsTopologySource ?? internalTopologySource;
  
  const handleTopologySourceChange = (source: TopologySource) => {
    if (onTopologySourceChange) onTopologySourceChange(source);
    else setInternalTopologySource(source);
  };

  const [internalTmAlgorithm, setInternalTmAlgorithm] = useState<string>('3d_slab_geom');
  const tmAlgorithm = propsTmAlgorithm ?? internalTmAlgorithm;
  
  const handleTmAlgorithmChange = (algo: string) => {
    if (onTmAlgorithmChange) onTmAlgorithmChange(algo);
    else setInternalTmAlgorithm(algo);
  };

  const [ssAlgorithm, setSsAlgorithm] = useState<string>('dssp');
  const [flowType, setFlowType] = useState<string>('ss_then_tm');
  const [customUniprotId, setCustomUniprotId] = useState<string>('');
  const [overlayType, setOverlayType] = useState<'none' | 'dssp' | 'stride'>('none');
  const [calculatedData, setCalculatedData] = useState<CalculatedTopologyData | null>(null);
  const [overlayResult, setOverlayResult] = useState<SecondaryStructureResult | null>(null);
  const [overlayError, setOverlayError] = useState<string | null>(null);
  const [loadingCalculated, setLoadingCalculated] = useState<boolean>(false);
  const [calculatedError, setCalculatedError] = useState<string | null>(null);
  const [showUniProtInfo, setShowUniProtInfo] = useState<boolean>(false);

  // --- Advanced TM parameters (user-tunable biological thresholds) ---
  // Empty = use the backend default (app.core.constants). The values the backend
  // actually used come back in `parameters_used` and are shown as placeholders, so
  // the form can never silently disagree with the server.
  const [showAdvancedParams, setShowAdvancedParams] = useState(false);
  const [tmThickness, setTmThickness] = useState<string>('');
  const [tmMinElement, setTmMinElement] = useState<string>('');
  const [tmMinCrossSpan, setTmMinCrossSpan] = useState<string>('');
  const [tmFullCrossFrac, setTmFullCrossFrac] = useState<string>('');
  const [tmBrokenGapMax, setTmBrokenGapMax] = useState<string>('');
  const [tmMinMembraneScore, setTmMinMembraneScore] = useState<string>('');

  const [colorDrawerOpen, setColorDrawerOpen] = useState<boolean>(false);
  const [activeTab, setActiveTab] = useState<'preset' | 'helices' | 'residues' | 'effects'>('preset');
  const [selectedPaletteKey, setSelectedPaletteKey] = useState<string>('PAPER_DEFAULT');
  const [customHelixColors, setCustomHelixColors] = useState<Record<string, string>>({});
  const [customResidueRules, setCustomResidueRules] = useState<CustomResidueColorRule[]>([]);
  const [hoverBrightness, setHoverBrightness] = useState<number>(1.1);
  const [hoverShadow, setHoverShadow] = useState<number>(0.3);

  const [resStartInput, setResStartInput] = useState<string>('');
  const [resEndInput, setResEndInput] = useState<string>('');
  const [resColorInput, setResColorInput] = useState<string>('#ff0055');

  const [hoveredElement, setHoveredElement] = useState<{
    title: string;
    range: string;
    length: number;
    details?: string;
  } | null>(null);

  const [figureTheme, setFigureTheme] = useState<FigureTheme>('publication');
  const [exporting, setExporting] = useState<boolean>(false);
  const svgRef = useRef<SVGSVGElement>(null);

  const fetchUniProtTopology = async (uniprotIdToFetch: string) => {
    if (!uniprotIdToFetch.trim()) return;
    setLoadingUniProt(true);
    setUniprotError(null);
    try {
      const response = await fetch(`${API_URL}/api/secondary-structure/uniprot/${uniprotIdToFetch.trim()}`);
      if (!response.ok) {
        const errData = await response.json().catch(() => ({}));
        throw new Error(errData.detail || `Could not load UniProt entry ${uniprotIdToFetch}`);
      }
      setUniprotData((await response.json()) as UniProtTopologyData);
    } catch (err: any) {
      setUniprotError(err.message || 'Could not load UniProt data');
      setUniprotData(null);
    } finally {
      setLoadingUniProt(false);
    }
  };

  const fetchCalculatedTopology = useCallback(async (
    filenameToFetch: string,
    tm: string,
    ss: string,
    flow: string,
    cUni: string,
    chainIdToFetch?: string
  ) => {
    if (!filenameToFetch.trim()) return;
    setLoadingCalculated(true);
    setCalculatedError(null);
    try {
      const qp = new URLSearchParams({ tm_algo: tm, ss_algo: ss, flow_type: flow });
      // Without chain_id the backend analyses the first protein chain, which is not
      // necessarily the chain selected here (SS overlay and 3D selection use this one).
      if (chainIdToFetch) qp.set('chain_id', chainIdToFetch);
      // Empty UniProt ID is allowed: the backend reads it from DBREF / _struct_ref.
      if (tm === 'uniprot_api' && cUni.trim()) qp.set('uniprot_id', cUni.trim().toUpperCase());
      const setNumber = (key: string, raw: string, integer = false) => {
        if (!raw.trim()) return;
        const value = integer ? parseInt(raw, 10) : parseFloat(raw);
        if (!isNaN(value)) qp.set(key, String(value));
      };
      // Read by the TM x SS flows. (min_tm_element / min_cross_span / full_cross_frac
      // are only read by the SS-element-first predictor, so they are not sent.)
      setNumber('thickness', tmThickness);
      setNumber('broken_gap_max', tmBrokenGapMax, true);
      setNumber('min_membrane_score', tmMinMembraneScore);

      const response = await fetch(
        `${API_URL}/api/secondary-structure/predict-topology/${encodeURIComponent(filenameToFetch.trim())}?${qp.toString()}`
      );
      if (!response.ok) {
        const errData = await response.json().catch(() => ({}));
        throw new Error(errData.detail || `Could not compute topology for ${filenameToFetch}`);
      }
      setCalculatedData((await response.json()) as CalculatedTopologyData);
    } catch (err: any) {
      setCalculatedError(err.message || 'Could not compute topology');
      setCalculatedData(null);
    } finally {
      setLoadingCalculated(false);
    }
  }, [tmThickness, tmBrokenGapMax, tmMinMembraneScore]);

  useEffect(() => {
    if (uniprotId) {
      setUniprotIdInput(uniprotId);
      fetchUniProtTopology(uniprotId);
    } else {
      fetchUniProtTopology('P31645');
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [uniprotId]);

  useEffect(() => {
    if (triggerTmRecalc && triggerTmRecalc > 0) {
      if (filename) {
        setCalculatedData(null);
        setCalculatedError(null);
        fetchCalculatedTopology(filename, tmAlgorithm, ssAlgorithm, flowType, customUniprotId, chain?.id);
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [triggerTmRecalc]);

  // The calculated topology belongs to one chain: recompute when another chain is selected.
  useEffect(() => {
    if (topologySource !== 'calculated' || !filename || !chain?.id || loadingCalculated) return;
    if (calculatedData?.chain_id && calculatedData.chain_id !== chain.id) {
      fetchCalculatedTopology(filename, tmAlgorithm, ssAlgorithm, flowType, customUniprotId, chain.id);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [chain?.id]);

  // Prefill the UniProt accession of the uploaded structure (still editable).
  useEffect(() => {
    if (uniprotId && !customUniprotId) setCustomUniprotId(uniprotId);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [uniprotId]);

  const activeTopologyData = topologySource === 'calculated' ? calculatedData : uniprotData;

  // UniProt topology is in UniProt sequence numbering; the uploaded structure (its
  // DSSP/STRIDE result, the 3D selection) is in the file's author numbering. The two
  // only agree by chance, so structure-derived data is combined with the map only
  // when the map itself comes from the structure (Calculated, or the SS fallback).
  // To put UniProt features on the structure use Calculated -> "UniProt API": the
  // backend aligns the UniProt sequence to the chain.
  const numberingMatchesStructure = topologySource === 'calculated' || !activeTopologyData;
  const structureSS = numberingMatchesStructure ? secondaryResult : null;
  const selectResidue = numberingMatchesStructure ? onSelectResidue : undefined;

  // The overlay uses the method picked in the overlay menu (it used to show whatever
  // method was last run on the page, whatever the menu said).
  useEffect(() => {
    if (overlayType === 'none' || !filename) {
      setOverlayResult(null);
      setOverlayError(null);
      return;
    }
    const wanted: 'DSSP' | 'STRIDE' = overlayType === 'dssp' ? 'DSSP' : 'STRIDE';
    if (secondaryResult && secondaryResult.method.toUpperCase() === wanted) {
      setOverlayResult(secondaryResult);
      setOverlayError(null);
      return;
    }
    let cancelled = false;
    setOverlayError(null);
    runSecondaryStructure(filename, wanted)
      .then((result) => {
        if (!cancelled) setOverlayResult(result);
      })
      .catch((err: unknown) => {
        if (cancelled) return;
        setOverlayResult(null);
        setOverlayError(err instanceof Error ? err.message : `${wanted} overlay failed`);
      });
    return () => {
      cancelled = true;
    };
  }, [overlayType, filename, secondaryResult]);
  const overlaySS = numberingMatchesStructure ? overlayResult : null;

  useEffect(() => {
    if (onTopologyDataChange) {
      onTopologyDataChange(activeTopologyData);
    }
  }, [activeTopologyData, onTopologyDataChange]);

  const activeError = topologySource === 'calculated' ? calculatedError : uniprotError;
  const activeLoading = topologySource === 'calculated' ? loadingCalculated : loadingUniProt;

  /* ---------------------------------------------------------------- *
   * Build helices + loops
   * ---------------------------------------------------------------- */

  const { helices, loops, helixCount } = useMemo(
    () =>
      buildTopologyModel({
        chain,
        activeTopologyData,
        secondaryResult,
        structureSS,
        topologySource,
        selectedPaletteKey,
        customHelixColors,
      }),
    [chain, secondaryResult, structureSS, activeTopologyData, topologySource, selectedPaletteKey, customHelixColors]
  );

  const handleExport = useCallback(
    async (format: 'png' | 'jpeg') => {
      if (!svgRef.current) return;
      setExporting(true);
      try {
        const id = activeTopologyData?.uniprot_id ?? 'topology';
        const ext = format === 'jpeg' ? 'jpg' : 'png';
        await exportSvgAsImage(svgRef.current, format, `${id}_TM_topology.${ext}`, 3);
      } catch {
        /* the user can press export again */
      } finally {
        setExporting(false);
      }
    },
    [activeTopologyData]
  );

  const handleAddResidueRule = () => {
    const start = parseInt(resStartInput.trim(), 10);
    const end = resEndInput.trim() ? parseInt(resEndInput.trim(), 10) : start;
    if (isNaN(start)) return;
    setCustomResidueRules((prev) => [
      ...prev,
      {
        id: `rule-${Date.now()}`,
        label: start === end ? `Residue ${start}` : `Residues ${start}–${end}`,
        startRes: Math.min(start, end),
        endRes: Math.max(start, end),
        color: resColorInput,
      },
    ]);
    setResStartInput('');
    setResEndInput('');
  };

  const handleRemoveResidueRule = (id: string) =>
    setCustomResidueRules((prev) => prev.filter((r) => r.id !== id));

  const handleResetColors = () => {
    setSelectedPaletteKey('PAPER_DEFAULT');
    setCustomHelixColors({});
    setCustomResidueRules([]);
  };

  const getCustomResidueColorForHelix = (startRes: number, endRes: number): string | null =>
    customResidueRules.find((rule) => rule.startRes <= endRes && rule.endRes >= startRes)?.color ?? null;

  if (!chain) {
    return <div className="empty-state">Upload a structure to draw the 2D transmembrane map.</div>;
  }

  const isPub = figureTheme === 'publication';
  const slabParamsActive = tmAlgorithm === '3d_slab_geom';
  const gapParamActive = flowType !== 'tm_then_ss' && ssAlgorithm !== 'none';
  const usedParam = (key: string): string => {
    const value = calculatedData?.parameters_used?.[key];
    return value === undefined || value === null ? 'default' : String(value);
  };
  const calcWarnings = topologySource === 'calculated' ? calculatedData?.warnings ?? [] : [];
  const noCrossings = !!activeTopologyData && !activeLoading && helices.length === 0;

  /* ---------------------------------------------------------------- *
   * Geometry
   * ---------------------------------------------------------------- */

  const helixWidth = 40;
  const colSpacing = 72;
  const leftMargin = 96;
  const columnCount = helixCount || 1;
  const canvasWidth = Math.max(960, leftMargin + columnCount * colSpacing + 130);
  const canvasHeight = 440;

  const membraneTopY = 152;
  const membraneBottomY = 258;
  const membraneMidY = (membraneTopY + membraneBottomY) / 2;
  const overhang = 26;
  const breakGap = 9; // half-gap between the two parts of a broken helix

  interface Pos {
    x: number;
    topY: number;
    bottomY: number;
    angle: number;
    nEndY: number;
    cEndY: number;
  }

  const helixPositions: Record<string, Pos> = {};

  helices.forEach((h) => {
    let x = leftMargin + h.column * colSpacing;
    
    // Stagger the halves of a discontinuous helix so they don't perfectly align
    if (h.isSplit) {
      x += h.partIndex === 0 ? -12 : 12;
    }

    const angle = h.column % 2 === 0 ? 4 : -4;

    let topY: number;
    let bottomY: number;

    if (!h.isSplit) {
      topY = membraneTopY - overhang;
      bottomY = membraneBottomY + overhang;
    } else {
      // Part a sits in the half nearest its entry side; part b takes the other half.
      const occupiesUpperHalf = h.partIndex === 0 ? h.entrySide === 'out' : h.exitSide === 'out';
      if (occupiesUpperHalf) {
        topY = membraneTopY - overhang;
        bottomY = membraneMidY - breakGap;
      } else {
        topY = membraneMidY + breakGap;
        bottomY = membraneBottomY + overhang;
      }
    }

    helixPositions[h.id] = {
      x,
      topY,
      bottomY,
      angle,
      nEndY: h.entrySide === 'out' ? topY : bottomY,
      cEndY: h.exitSide === 'out' ? topY : bottomY,
    };
  });

  const firstHelix = helices[0];
  const lastHelix = helices[helices.length - 1];

  return (
    <div 
      className="tm-topology-panel panel"
      style={{
        '--hover-brightness': hoverBrightness,
        '--hover-shadow-opacity': hoverShadow,
      } as React.CSSProperties}
    >
      {/* Header + toolbar */}
      <div className="panel-heading" style={{ flexWrap: 'wrap', gap: '12px' }}>
        <div>
          <span className="section-kicker">
            {topologySource === 'calculated'
              ? 'Calculated topology (DSSP/STRIDE elements + membrane fit)'
              : 'UniProt topology'}
          </span>
          <h2>Transmembrane secondary structure map</h2>
          <p className="panel-subtitle">
            {isCalculatedTopology(activeTopologyData)
              ? `${activeTopologyData.labeler} · chain ${activeTopologyData.chain_id ?? '?'} · membrane score ${
                  (activeTopologyData.membrane_score ?? 0).toFixed(2)
                } · ${helixCount} crossings, ${loops.length} loops`
              : activeTopologyData
              ? `${activeTopologyData.protein_name}${
                  activeTopologyData.gene_name ? ` (${activeTopologyData.gene_name})` : ''
                } · ${activeTopologyData.uniprot_id} · ${helixCount} helices, ${loops.length} loops · UniProt numbering`
              : 'Transmembrane helices, extracellular and cytoplasmic loops, re-entrant segments.'}
          </p>
        </div>

        <div className="tm-toolbar-actions">
          <div className="tm-preset-chip" style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
            <button
              className={`tm-tab-btn ${figureTheme === 'publication' ? 'active' : ''}`}
              onClick={() => setFigureTheme('publication')}
            >
              Publication
            </button>
            <button
              className={`tm-tab-btn ${figureTheme === 'lab' ? 'active' : ''}`}
              onClick={() => setFigureTheme('lab')}
            >
              Dark lab
            </button>
          </div>

          <div className="tm-preset-chip" style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
            <button
              className={`tm-tab-btn ${topologySource === 'uniprot' ? 'active' : ''}`}
              onClick={() => handleTopologySourceChange('uniprot')}
            >
              UniProt
            </button>
            <button
              className={`tm-tab-btn ${topologySource === 'calculated' ? 'active' : ''}`}
              onClick={() => handleTopologySourceChange('calculated')}
              disabled={!filename}
              title={!filename ? 'Upload a structure file to calculate topology' : undefined}
            >
              Calculated (beta)
            </button>
          </div>

          <button className="tm-color-toggle-btn" onClick={() => handleExport('png')} disabled={exporting}>
            <span>{exporting ? 'Exporting…' : 'Export PNG'}</span>
          </button>

          
          <select 
            className="tm-input-field" 
            value={overlayType}
            onChange={(e) => setOverlayType(e.target.value as 'none' | 'dssp' | 'stride')}
            style={{ padding: '4px 8px', width: 'auto' }}
            disabled={!numberingMatchesStructure || !filename}
            title={
              numberingMatchesStructure
                ? 'Secondary Structure Overlay'
                : 'The overlay uses structure numbering; UniProt topology uses UniProt numbering. Use Calculated → UniProt API to map UniProt onto the structure.'
            }
          >
            <option value="none">No Overlay</option>
            <option value="dssp">Overlay DSSP</option>
            <option value="stride">Overlay STRIDE</option>
          </select>
          <button
            className={`tm-color-toggle-btn ${colorDrawerOpen ? 'active' : ''}`}
            onClick={() => setColorDrawerOpen((open) => !open)}
          >
            <span>{colorDrawerOpen ? 'Close colors' : 'Customize colors'}</span>
          </button>

          {topologySource === 'uniprot' ? (
            <div style={{ display: 'flex', alignItems: 'center', gap: '6px', position: 'relative' }}>
              <input
                type="text"
                value={uniprotIdInput}
                onChange={(e) => setUniprotIdInput(e.target.value)}
                placeholder="e.g. P31645"
                className="tm-input-field"
                style={{ width: '96px', fontWeight: 700, textTransform: 'uppercase' }}
                onKeyDown={(e) => {
                  if (e.key === 'Enter') fetchUniProtTopology(uniprotIdInput);
                }}
              />
              <button
                onClick={() => fetchUniProtTopology(uniprotIdInput)}
                disabled={loadingUniProt}
                className="tm-add-btn"
                style={{ background: '#3b82f6' }}
              >
                {loadingUniProt ? '…' : 'Load'}
              </button>

              <button
                onClick={() => setShowUniProtInfo((open) => !open)}
                title="How the UniProt lookup works"
                aria-label="How the UniProt lookup works"
                className="tm-info-btn"
                style={{
                  background: showUniProtInfo ? '#38bdf8' : '#1e293b',
                  color: showUniProtInfo ? '#0f172a' : '#94a3b8',
                }}
              >
                ?
              </button>

              {showUniProtInfo && (
                <div className="tm-info-popover">
                  <strong>How the UniProt lookup works</strong>
                  Enter any UniProt accession (for example <code>P31645</code>) and press Load. The app reads that
                  entry's curated Transmembrane, Topological domain and Intramembrane features and draws the map from
                  them — no structure file needed. For a structure-derived estimate, switch to Calculated (beta).
                  <button onClick={() => setShowUniProtInfo(false)}>Close</button>
                </div>
              )}
            </div>
                    ) : (
            <>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '8px', border: '1px solid #333', padding: '8px', borderRadius: '4px' }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '8px', flexWrap: 'wrap' }}>
                  <select 
                    className="tm-input-field" 
                    value={tmAlgorithm}
                    onChange={(e) => handleTmAlgorithmChange(e.target.value)}
                    style={{ padding: '4px 8px' }}
                  >
                    <option value="3d_slab_geom">3D Slab Geometry (Recommended)</option>
                    <option value="kyte_doolittle_seq">Kyte-Doolittle Sequence</option>
                    <option value="uniprot_api">UniProt API</option>
                  </select>

                  {tmAlgorithm === 'uniprot_api' && (
                    <input 
                      type="text"
                      className="tm-input-field"
                      placeholder="UniProt ID (auto from file)"
                      value={customUniprotId}
                      onChange={(e) => setCustomUniprotId(e.target.value)}
                      style={{ width: '150px', padding: '4px 8px' }}
                    />
                  )}

                  <select 
                    className="tm-input-field" 
                    value={ssAlgorithm}
                    onChange={(e) => setSsAlgorithm(e.target.value)}
                    style={{ padding: '4px 8px' }}
                  >
                    <option value="dssp">DSSP</option>
                    <option value="stride">STRIDE</option>
                    <option value="none">None (Only TM boundaries)</option>
                  </select>

                  <select 
                    className="tm-input-field" 
                    value={flowType}
                    onChange={(e) => setFlowType(e.target.value)}
                    style={{ padding: '4px 8px' }}
                  >
                    <option value="ss_then_tm">Filter TM by SS (Recommended)</option>
                    <option value="tm_then_ss">Filter SS by TM</option>
                    <option value="parallel_merge">Parallel Merge (Strict Intersection)</option>
                  </select>
                </div>
                
                <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                  <button
                    onClick={() =>
                      filename && fetchCalculatedTopology(filename, tmAlgorithm, ssAlgorithm, flowType, customUniprotId, chain?.id)
                    }
                    disabled={loadingCalculated || !filename}
                    className="tm-add-btn"
                  >
                    {loadingCalculated ? 'Computing…' : 'Recalculate'}
                  </button>
                  <button
                    onClick={() => setShowAdvancedParams(!showAdvancedParams)}
                    className="tm-add-btn"
                    style={{ fontSize: '0.8em', opacity: 0.8 }}
                    title="Tune biological thresholds for TM detection"
                  >
                    {showAdvancedParams ? '▲ Parameters' : '▼ Parameters'}
                  </button>
                  {tmAlgorithm === 'uniprot_api' && !customUniprotId && (
                    <span style={{ color: '#94a3b8', fontSize: '12px' }}>
                      Empty = read from the file (DBREF / _struct_ref)
                    </span>
                  )}
                </div>
              </div>
              {showAdvancedParams && (
                <div style={{
                  margin: '8px 0', padding: '10px 14px',
                  background: 'rgba(100,100,140,0.08)', borderRadius: '8px',
                  display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '8px 16px',
                  fontSize: '0.82em',
                }}>
                  <label title="Hydrophobic core thickness (Å). Bacterial IM ~27, eukaryotic PM ~30, ER ~25. (Mitra 2004; OPM database). Used by 3D Slab Geometry." style={{ opacity: slabParamsActive ? 1 : 0.45 }}>
                    Membrane thickness (Å)
                    <input type="number" step="0.5" min="20" max="40" value={tmThickness}
                      placeholder={usedParam('membrane_thickness')} disabled={!slabParamsActive}
                      onChange={e => setTmThickness(e.target.value)}
                      className="tm-input-field" style={{ width: '70px', marginLeft: 4 }} />
                  </label>
                  <label title="Min mean hydrophobicity inside the slab. Below this the protein is treated as soluble. Used by 3D Slab Geometry." style={{ opacity: slabParamsActive ? 1 : 0.45 }}>
                    Min membrane score
                    <input type="number" step="0.1" min="-3" max="3" value={tmMinMembraneScore}
                      placeholder={usedParam('min_membrane_score')} disabled={!slabParamsActive}
                      onChange={e => setTmMinMembraneScore(e.target.value)}
                      className="tm-input-field" style={{ width: '60px', marginLeft: 4 }} />
                  </label>
                  <label title="Kinks/breaks up to this many residues stay inside one crossing when the chain keeps its direction. Used by the SS-driven flows." style={{ opacity: gapParamActive ? 1 : 0.45 }}>
                    Broken gap max
                    <input type="number" step="1" min="0" max="20" value={tmBrokenGapMax}
                      placeholder={usedParam('broken_gap_max')} disabled={!gapParamActive}
                      onChange={e => setTmBrokenGapMax(e.target.value)}
                      className="tm-input-field" style={{ width: '50px', marginLeft: 4 }} />
                  </label>
                  {/* Only read by the SS-element-first predictor (predict_topology_ss_first),
                      not by the TM x SS flows used here - shown for reference, not sent. */}
                  <label title="SS-first predictor only — not used by the TM × SS flows." style={{ opacity: 0.45 }}>
                    Min element in slab
                    <input type="number" value={tmMinElement} disabled
                      placeholder={usedParam('min_tm_element_in_slab')}
                      onChange={e => setTmMinElement(e.target.value)}
                      className="tm-input-field" style={{ width: '50px', marginLeft: 4 }} />
                  </label>
                  <label title="SS-first predictor only — not used by the TM × SS flows." style={{ opacity: 0.45 }}>
                    Min cross span
                    <input type="number" value={tmMinCrossSpan} disabled
                      placeholder={usedParam('min_cross_span_frac')}
                      onChange={e => setTmMinCrossSpan(e.target.value)}
                      className="tm-input-field" style={{ width: '60px', marginLeft: 4 }} />
                  </label>
                  <label title="SS-first predictor only — not used by the TM × SS flows." style={{ opacity: 0.45 }}>
                    Full cross frac
                    <input type="number" value={tmFullCrossFrac} disabled
                      placeholder={usedParam('full_cross_frac')}
                      onChange={e => setTmFullCrossFrac(e.target.value)}
                      className="tm-input-field" style={{ width: '60px', marginLeft: 4 }} />
                  </label>
                  <div style={{ gridColumn: '1 / -1', marginTop: 4, opacity: 0.65, fontSize: '0.9em' }}>
                    Empty fields use the backend defaults (shown greyed, as last used by the server).
                    Defaults are literature-derived (Kyte–Doolittle 1982; Mitra 2004; OPM/PDBTM).
                  </div>
                </div>
              )}
            </>
          )}
        </div>
      </div>

      {topologySource === 'uniprot' && (
        <div className="tm-preset-row" style={{ margin: '4px 0 12px' }}>
          <span className="preset-label">Examples</span>
          {DEFAULT_PRESETS.map((p) => (
            <button
              key={p.id}
              className="tm-preset-chip"
              onClick={() => {
                setUniprotIdInput(p.id);
                fetchUniProtTopology(p.id);
              }}
            >
              {p.label}
            </button>
          ))}
        </div>
      )}

      {/* Color drawer */}
      {colorDrawerOpen && (
        <div className={`tm-color-customizer-drawer ${isPub ? 'publication' : 'lab'}`}>
          <div className="tm-drawer-tabs">
            <button
              className={`tm-tab-btn ${activeTab === 'preset' ? 'active' : ''}`}
              onClick={() => setActiveTab('preset')}
            >
              Palettes
            </button>
            <button
              className={`tm-tab-btn ${activeTab === 'helices' ? 'active' : ''}`}
              onClick={() => setActiveTab('helices')}
            >
              Individual helices
            </button>
            <button
              className={`tm-tab-btn ${activeTab === 'residues' ? 'active' : ''}`}
              onClick={() => setActiveTab('residues')}
            >
              Residue ranges
            </button>
            <button
              className={`tm-tab-btn ${activeTab === 'effects' ? 'active' : ''}`}
              onClick={() => setActiveTab('effects')}
            >
              Effects
            </button>
            <button className="tm-reset-btn" onClick={handleResetColors}>
              Reset
            </button>
          </div>

          {activeTab === 'effects' && (
            <div className="tm-effects-section" style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
                <label style={{ fontSize: '11px', color: isPub ? '#334155' : '#94a3b8', fontWeight: 600, display: 'flex', justifyContent: 'space-between' }}>
                  <span>Hover Brightness</span>
                  <span>{hoverBrightness.toFixed(2)}</span>
                </label>
                <input 
                  type="range" 
                  min="1.0" 
                  max="2.0" 
                  step="0.05" 
                  value={hoverBrightness} 
                  onChange={(e) => setHoverBrightness(parseFloat(e.target.value))} 
                />
              </div>
              <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
                <label style={{ fontSize: '11px', color: isPub ? '#334155' : '#94a3b8', fontWeight: 600, display: 'flex', justifyContent: 'space-between' }}>
                  <span>Hover Shadow Opacity</span>
                  <span>{hoverShadow.toFixed(2)}</span>
                </label>
                <input 
                  type="range" 
                  min="0.0" 
                  max="1.0" 
                  step="0.05" 
                  value={hoverShadow} 
                  onChange={(e) => setHoverShadow(parseFloat(e.target.value))} 
                />
              </div>
            </div>
          )}

          {activeTab === 'preset' && (
            <div className="tm-palette-grid">
              {Object.entries(PALETTES).map(([key, palette]) => (
                <button
                  key={key}
                  className={`tm-palette-btn ${selectedPaletteKey === key ? 'active' : ''}`}
                  onClick={() => {
                    setSelectedPaletteKey(key);
                    setCustomHelixColors({});
                  }}
                >
                  <div className="tm-palette-preview">
                    {buildPalette(palette.colors, Math.max(6, Math.min(12, columnCount))).map((c, i) => (
                      <span key={i} style={{ backgroundColor: c }} />
                    ))}
                  </div>
                  <span>{palette.label}</span>
                </button>
              ))}
            </div>
          )}

          {activeTab === 'helices' && (
            <div className="tm-helix-color-grid">
              {helices.map((h) => (
                <div key={`color-${h.id}`} className="tm-helix-color-item">
                  <span>TM{h.subLabel}</span>
                  <input
                    type="color"
                    className="tm-color-input"
                    value={h.color}
                    onChange={(e) =>
                      setCustomHelixColors((prev) => ({ ...prev, [h.subLabel]: e.target.value }))
                    }
                  />
                </div>
              ))}
            </div>
          )}

          {activeTab === 'residues' && (
            <div className="tm-residue-color-section">
              <div className="tm-residue-form">
                <input
                  type="number"
                  className="tm-input-field"
                  placeholder="Start"
                  value={resStartInput}
                  onChange={(e) => setResStartInput(e.target.value)}
                />
                <span>–</span>
                <input
                  type="number"
                  className="tm-input-field"
                  placeholder="End (optional)"
                  value={resEndInput}
                  onChange={(e) => setResEndInput(e.target.value)}
                />
                <input
                  type="color"
                  className="tm-color-input"
                  value={resColorInput}
                  onChange={(e) => setResColorInput(e.target.value)}
                />
                <button className="tm-add-btn" onClick={handleAddResidueRule}>
                  Add highlight
                </button>
              </div>

              {customResidueRules.length > 0 ? (
                <div className="tm-residue-tag-list">
                  {customResidueRules.map((rule) => (
                    <span key={rule.id} className="tm-residue-tag" style={{ borderLeftColor: rule.color }}>
                      <span className="tm-rule-color-dot" style={{ backgroundColor: rule.color }} />
                      {rule.label}
                      <button onClick={() => handleRemoveResidueRule(rule.id)}>×</button>
                    </span>
                  ))}
                </div>
              ) : (
                <small style={{ color: '#64748b', fontSize: '11px' }}>
                  Enter a residue number or range to highlight positions on the map.
                </small>
              )}
            </div>
          )}
        </div>
      )}

      {activeError && <div className="tm-error-banner">{activeError}</div>}
      {overlayError && <div className="tm-error-banner">Overlay: {overlayError}</div>}
      {calcWarnings.length > 0 && (
        <div className="tm-error-banner" style={{ background: 'rgba(234, 179, 8, 0.12)', borderColor: '#ca8a04', color: isPub ? '#713f12' : '#fde68a' }}>
          <strong>Backend notes</strong>
          <ul style={{ margin: '4px 0 0 16px', padding: 0 }}>
            {calcWarnings.map((w, i) => (
              <li key={i}>{w}</li>
            ))}
          </ul>
        </div>
      )}
      {noCrossings && (
        <div className="tm-error-banner tm-loading-banner">
          No transmembrane segment in this topology
          {activeTopologyData?.protein_name ? ` — ${activeTopologyData.protein_name}` : ''}.
        </div>
      )}
      {activeLoading && !activeTopologyData && (
        <div className="tm-error-banner tm-loading-banner">
          {topologySource === 'calculated' ? 'Computing topology from the structure…' : 'Loading UniProt topology…'}
        </div>
      )}

      {/* Diagram */}
      <div className={`tm-diagram-wrap ${isPub ? 'publication' : 'lab'}`}>
        <svg
          ref={svgRef}
          className="tm-diagram-svg"
          viewBox={`0 0 ${canvasWidth} ${canvasHeight}`}
          role="img"
          aria-label="Transmembrane topology diagram"
        >
          <rect x="0" y="0" width={canvasWidth} height={canvasHeight} fill={isPub ? '#ffffff' : 'transparent'} />

          <defs>
            {/* One barrel gradient per helix: shade, body, highlight, shade. */}
            {helices.map((h) => {
              const c = getCustomResidueColorForHelix(h.startRes, h.endRes) ?? h.color;
              return (
                <linearGradient key={`grad-${h.id}`} id={`cyl-${h.id}`} x1="0%" y1="0%" x2="100%" y2="0%">
                  <stop offset="0%" stopColor={darken(c, 0.4)} />
                  <stop offset="16%" stopColor={c} />
                  <stop offset="38%" stopColor={lighten(c, 0.55)} />
                  <stop offset="60%" stopColor={c} />
                  <stop offset="100%" stopColor={darken(c, 0.32)} />
                </linearGradient>
              );
            })}
            <linearGradient id="short-helix-grad" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stopColor="#94a3b8" />
              <stop offset="35%" stopColor="#e2e8f0" />
              <stop offset="70%" stopColor="#94a3b8" />
              <stop offset="100%" stopColor="#64748b" />
            </linearGradient>

            {/* Break gradients for split helices */}
            {helices.map((h) => {
              if (!h.isSplit || h.partIndex !== 0) return null;
              const partB = helices.find((o) => o.column === h.column && o.partIndex === 1);
              if (!partB) return null;
              const a = helixPositions[h.id];
              const b = helixPositions[partB.id];
              if (!a || !b) return null;
              const isAAbove = a.topY < b.topY;
              const topColor = isAAbove ? h.color : partB.color;
              const bottomColor = isAAbove ? partB.color : h.color;
              return (
                <linearGradient key={`break-grad-${h.id}`} id={`break-grad-${h.id}`} x1="0%" y1="0%" x2="0%" y2="100%">
                  <stop offset="0%" stopColor={topColor} />
                  <stop offset="100%" stopColor={bottomColor} />
                </linearGradient>
              );
            })}
          
            <pattern id="hatch-warning" width="8" height="8" patternTransform="rotate(45 0 0)" patternUnits="userSpaceOnUse">
              <line x1="0" y1="0" x2="0" y2="8" stroke="#ff0000" strokeWidth="3" opacity="0.8" />
            </pattern>

          </defs>

          {/* Lipid bilayer */}
          <g className="membrane-zone">
            <rect
              x="0"
              y={membraneTopY}
              width={canvasWidth}
              height={membraneBottomY - membraneTopY}
              fill={isPub ? '#d9d4cb' : '#1e293b'}
              opacity={isPub ? 0.75 : 0.85}
            />
            <text x="16" y={membraneTopY - 12} className="membrane-label">
              Extracellular
            </text>
            <text x="16" y={membraneBottomY + 24} className="membrane-label">
              Cytoplasmic
            </text>
          </g>

          {/* Loops */}
          <g className="tm-loops">
            {loops.map((loop) => {
              const posPrev = helixPositions[loop.prevHelixId];
              const posNext = helixPositions[loop.nextHelixId];
              if (!posPrev || !posNext) return null;

              const isEL = loop.type === 'EL';
              const startX = posPrev.x + helixWidth / 2;
              const startY = posPrev.cEndY;
              const endX = posNext.x + helixWidth / 2;
              const endY = posNext.nEndY;
              const midX = (startX + endX) / 2;

              let archDepth = 26 + Math.min(46, loop.length * 0.7);
              if (loop.hasExtraFeature) archDepth += 42;

              const apexY = isEL
                ? Math.max(42, Math.min(startY, endY) - archDepth)
                : Math.min(canvasHeight - 34, Math.max(startY, endY) + archDepth);

              const pathD = `M ${startX} ${startY} C ${startX} ${apexY}, ${endX} ${apexY}, ${endX} ${endY}`;
              const loopColorStart = helices.find((h) => h.id === loop.prevHelixId)?.color ?? '#64748b';
              const loopColorEnd = helices.find((h) => h.id === loop.nextHelixId)?.color ?? '#64748b';

              const labelY = isEL
                ? (loop.hasExtraFeature ? apexY + 34 : apexY - 8)
                : (loop.hasExtraFeature ? apexY - 26 : apexY + 16);
              const hideLabel = loop.extraFeatures?.length === 1 && loop.extraFeatures[0].label === loop.label;

              return (
                <g
                  key={loop.id}
                  className="loop-group"
                  onMouseEnter={() =>
                    setHoveredElement({
                      title: `${loop.label} — ${isEL ? 'extracellular' : 'cytoplasmic'} loop`,
                      range: `Residues ${loop.startRes}–${loop.endRes}`,
                      length: loop.length,
                      details: loop.domainName || undefined,
                    })
                  }
                  onMouseLeave={() => setHoveredElement(null)}
                  onClick={() => selectResidue?.(loop.startRes)}
                >
                  <defs>
                    <linearGradient id={`loop-grad-${loop.id}`} x1="0%" y1="0%" x2="100%" y2="0%">
                      <stop offset="0%" stopColor={loopColorStart} />
                      <stop offset="100%" stopColor={loopColorEnd} />
                    </linearGradient>
                  
            <pattern id="hatch-warning" width="8" height="8" patternTransform="rotate(45 0 0)" patternUnits="userSpaceOnUse">
              <line x1="0" y1="0" x2="0" y2="8" stroke="#ff0000" strokeWidth="3" opacity="0.8" />
            </pattern>

          </defs>
                  
                  <path d={pathD} fill="none" stroke={`url(#loop-grad-${loop.id})`} strokeWidth={2.6} strokeLinecap="round" />

                  {!hideLabel && (
                    <text
                      x={midX}
                      y={labelY}
                      className="loop-text-label"
                      textAnchor="middle"
                      style={{ fill: isPub ? '#1e293b' : '#e2e8f0' }}
                    >
                      {loop.label}
                    </text>
                  )}

                  {loop.hasExtraFeature && loop.extraFeatures && (
                    <g className="short-helices">
                      {loop.extraFeatures.map((sh, idx) => {
                        const w = 48; // a bit wider for labels
                        const shX = midX - w / 2 + (sh.offsetFactor ?? 0) * (w + 10);
                        const shY = isEL ? apexY - 6 : apexY - 12;
                        return (
                          <g
                            key={`${loop.id}-sh-${idx}`}
                            transform={`translate(${shX}, ${shY})`}
                            onMouseEnter={(e) => {
                              e.stopPropagation();
                              setHoveredElement({
                                title: `${sh.type}: ${sh.label}`,
                                range: `Residues ${sh.startRes}–${sh.endRes}`,
                                length: sh.endRes - sh.startRes + 1,
                                details: `Inside ${loop.label}`,
                              });
                            }}
                          >
                            {sh.type === 'Strand' || sh.type === 'β-strand' ? (
                              // β-strand outside the membrane: draw as an arrow
                              <path
                                d={`M 0 4 L ${w - 12} 4 L ${w - 12} 0 L ${w} 10.5 L ${w - 12} 21 L ${w - 12} 17 L 0 17 Z`}
                                fill={isPub ? '#e0a64b' : '#c98a2b'}
                                stroke={isPub ? '#475569' : '#0f172a'}
                                strokeWidth="1"
                              />
                            ) : (
                              // α-helix (or UniProt feature) outside the membrane: rounded cylinder
                              <rect
                                x="0"
                                y="0"
                                width={w}
                                height="21"
                                rx={sh.type === 'Helix' || sh.type === 'Intramembrane' ? '10' : '4'}
                                fill="url(#short-helix-grad)"
                                stroke={isPub ? '#475569' : '#0f172a'}
                                strokeWidth="1"
                              />
                            )}
                            <text
                              x={w / 2}
                              y="15"
                              className="short-helix-text"
                              textAnchor="middle"
                              style={{ fill: '#0f172a', fontSize: '10px' }}
                            >
                              {sh.label.substring(0, 7)}
                            </text>
                          </g>
                        );
                      })}
                    </g>
                  )}
                </g>
              );
            })}
          </g>

          
          {/* N-terminus and C-terminus tails */}
          {(() => {
            if (helices.length === 0) return null;
            const first = helices[0];
            const last = helices[helices.length - 1];
            const posFirst = helixPositions[first.id];
            const posLast = helixPositions[last.id];
            if (!posFirst || !posLast) return null;

            // Residue NUMBERS, not sequence length: numbering rarely starts at 1 (the old
            // `sequence.length` bound dropped C-terminal features of e.g. residues 25-144).
            const nFeatures = getExtraFeatures(Number.NEGATIVE_INFINITY, first.startRes - 1, activeTopologyData, structureSS, chain?.id) || [];
            const cFeatures = getExtraFeatures(last.endRes + 1, Number.POSITIVE_INFINITY, activeTopologyData, structureSS, chain?.id) || [];

            return (
              <g className="tm-terminals">
                {/* N-terminus */}
                <path
                  d={`M 20 ${first.entrySide === 'out' ? membraneTopY - 60 : membraneBottomY + 60} Q ${posFirst.x / 2} ${first.entrySide === 'out' ? membraneTopY - 30 : membraneBottomY + 30} ${posFirst.x + helixWidth / 2} ${posFirst.nEndY}`}
                  fill="none"
                  stroke={first.color}
                  strokeWidth="2.6"
                />
                <text x="20" y={first.entrySide === 'out' ? membraneTopY - 70 : membraneBottomY + 75} textAnchor="middle" style={{ fill: isPub ? '#1e293b' : '#e2e8f0', fontSize: '11px', fontWeight: 'bold' }}>N (NH2)</text>
                
                {nFeatures.map((f, i) => {
                   const y = first.entrySide === 'out' ? membraneTopY - 45 : membraneBottomY + 45;
                   const x = 30 + i * 50;
                   return (
                     <g key={`n-${i}`} transform={`translate(${x}, ${y})`} onMouseEnter={(e) => {
                       e.stopPropagation();
                       setHoveredElement({ title: `${f.type}: ${f.label}`, range: `Residues ${f.startRes}-${f.endRes}`, length: f.endRes - f.startRes + 1, details: 'N-terminus' });
                     }} onMouseLeave={() => setHoveredElement(null)}>
                       <rect x="0" y="0" width="45" height="18" rx={f.type === 'Helix' ? 9 : 4} fill="url(#short-helix-grad)" stroke={isPub ? '#475569' : '#0f172a'} />
                       <text x="22.5" y="13" textAnchor="middle" style={{ fill: '#0f172a', fontSize: '9px' }}>{f.label.substring(0, 7)}</text>
                     </g>
                   );
                })}

                {/* C-terminus */}
                <path
                  d={`M ${posLast.x + helixWidth / 2} ${posLast.cEndY} Q ${(posLast.x + canvasWidth) / 2} ${last.exitSide === 'out' ? membraneTopY - 30 : membraneBottomY + 30} ${canvasWidth - 30} ${last.exitSide === 'out' ? membraneTopY - 60 : membraneBottomY + 60}`}
                  fill="none"
                  stroke={last.color}
                  strokeWidth="2.6"
                />
                <text x={canvasWidth - 30} y={last.exitSide === 'out' ? membraneTopY - 70 : membraneBottomY + 75} textAnchor="middle" style={{ fill: isPub ? '#1e293b' : '#e2e8f0', fontSize: '11px', fontWeight: 'bold' }}>C (COOH)</text>
                
                {cFeatures.map((f, i) => {
                   const y = last.exitSide === 'out' ? membraneTopY - 45 : membraneBottomY + 45;
                   const x = canvasWidth - 80 - i * 50;
                   return (
                     <g key={`c-${i}`} transform={`translate(${x}, ${y})`} onMouseEnter={(e) => {
                       e.stopPropagation();
                       setHoveredElement({ title: `${f.type}: ${f.label}`, range: `Residues ${f.startRes}-${f.endRes}`, length: f.endRes - f.startRes + 1, details: 'C-terminus' });
                     }} onMouseLeave={() => setHoveredElement(null)}>
                       <rect x="0" y="0" width="45" height="18" rx={f.type === 'Helix' ? 9 : 4} fill="url(#short-helix-grad)" stroke={isPub ? '#475569' : '#0f172a'} />
                       <text x="22.5" y="13" textAnchor="middle" style={{ fill: '#0f172a', fontSize: '9px' }}>{f.label.substring(0, 7)}</text>
                     </g>
                   );
                })}
              </g>
            );
          })()}

          {/* Break connectors between the two halves of a discontinuous helix */}
          <g className="tm-breaks">
            {helices.map((h) => {
              if (!h.isSplit || h.partIndex !== 0) return null;
              const partB = helices.find((o) => o.column === h.column && o.partIndex === 1);
              if (!partB) return null;
              const a = helixPositions[h.id];
              const b = helixPositions[partB.id];
              if (!a || !b) return null;

              // Helper to get true absolute center of the cap after rotation
              const getCapCenter = (pos: typeof a, isBottom: boolean) => {
                const cx = helixWidth / 2;
                const cylHeight = Math.abs(pos.bottomY - pos.topY);
                const cy = cylHeight / 2;
                const angleRad = (pos.angle * Math.PI) / 180;
                const absX = pos.x + cx + (isBottom ? -cy : cy) * Math.sin(angleRad);
                const absY = pos.topY + cy + (isBottom ? cy : -cy) * Math.cos(angleRad);
                return { x: absX, y: absY };
              };

              const isAAbove = a.topY < b.topY;
              const ptA = getCapCenter(a, isAAbove); // bottom of a, or top of a
              const ptB = getCapCenter(b, !isAAbove); // top of b, or bottom of b
              
              const midY = (ptA.y + ptB.y) / 2;

              return (
                <path
                  key={`break-${h.id}`}
                  d={`M ${ptA.x} ${ptA.y} C ${ptA.x} ${midY}, ${ptB.x} ${midY}, ${ptB.x} ${ptB.y}`}
                  fill="none"
                  stroke={`url(#break-grad-${h.id})`}
                  strokeWidth="3.5"
                  strokeLinecap="round"
                />
              );
            })}
          </g>

          {/* Helices as 3D barrels */}
          <g className="tm-helices">
            {helices.map((h) => {
              const pos = helixPositions[h.id];
              if (!pos) return null;

              const cylHeight = Math.abs(pos.bottomY - pos.topY);
              const isSelected =
                numberingMatchesStructure &&
                selectedResidue != null && selectedResidue >= h.startRes && selectedResidue <= h.endRes;
              const resCustomColor = getCustomResidueColorForHelix(h.startRes, h.endRes);
              const baseColor = resCustomColor ?? h.color;
              const cx = helixWidth / 2;
              const cy = cylHeight / 2;
              const labelColor = getContrastTextColor(baseColor);

              const isBeta = h.isBeta ?? isBetaStrandDescription(h.description);
              const isAlpha = !isBeta;
              
              const pointsDown = pos.nEndY < pos.cEndY;
              const headH = Math.min(20, cylHeight / 2);
              const arrowPath = pointsDown
                ? `M 6 0 L 6 ${cylHeight - headH} L 0 ${cylHeight - headH} L ${cx} ${cylHeight} L ${helixWidth} ${cylHeight - headH} L ${helixWidth - 6} ${cylHeight - headH} L ${helixWidth - 6} 0 Z`
                : `M 6 ${cylHeight} L 6 ${headH} L 0 ${headH} L ${cx} 0 L ${helixWidth} ${headH} L ${helixWidth - 6} ${headH} L ${helixWidth - 6} ${cylHeight} Z`;

              return (
                <g
                  key={h.id}
                  className={`helix-group ${isSelected ? 'selected' : ''}`}
                  transform={`translate(${pos.x}, ${pos.topY}) rotate(${pos.angle}, ${cx}, ${cy})`}
                  onMouseEnter={() =>
                    setHoveredElement({
                      title: `${isBeta ? 'Beta Strand' : 'Helix'} TM${h.subLabel}`,
                      range: `Residues ${h.startRes}–${h.endRes}`,
                      length: h.length,
                      details: h.isSplit
                        ? `Discontinuous segment, part ${h.partIndex === 0 ? '1' : '2'}`
                        : h.description || `Transmembrane ${isBeta ? 'beta strand' : 'alpha helix'}`,
                    })
                  }
                  onMouseLeave={() => setHoveredElement(null)}
                  onClick={() => selectResidue?.(h.startRes)}
                >
                  {isAlpha ? (
                    <>
                      <rect
                        x="0"
                        y="6"
                        width={helixWidth}
                        height={Math.max(0, cylHeight - 12)}
                        fill={`url(#cyl-${h.id})`}
                        className="cylinder-body"
                      />
                      <ellipse
                        cx={cx}
                        cy={cylHeight - 6}
                        rx={helixWidth / 2}
                        ry="6"
                        fill={darken(baseColor, 0.28)}
                      />
                      <ellipse
                        cx={cx}
                        cy="6"
                        rx={helixWidth / 2}
                        ry="6"
                        fill={lighten(baseColor, 0.28)}
                        stroke={darken(baseColor, 0.25)}
                        strokeWidth="0.8"
                      />
                      <rect
                        x="0"
                        y="6"
                        width={helixWidth}
                        height={Math.max(0, cylHeight - 12)}
                        fill="none"
                        stroke={isSelected ? '#ff6f61' : darken(baseColor, 0.3)}
                        strokeWidth={isSelected ? 2.5 : 0.8}
                      />
                    </>
                  ) : (
                    <path
                      d={arrowPath}
                      fill={`url(#cyl-${h.id})`}
                      stroke={isSelected ? '#ff6f61' : darken(baseColor, 0.3)}
                      strokeWidth={isSelected ? 2.5 : 1.2}
                    />
                  )}

                  
                  {/* Disagreement Overlay */}
                  {overlayType !== 'none' && overlaySS && (() => {
                     const ssResidues = overlaySS.residues
                       .filter(r => r.chain_id === chain?.id && r.residue_number >= h.startRes && r.residue_number <= h.endRes)
                       .sort(compareResidues);
                     // The cylinder is drawn top-to-bottom; when the N-terminus sits at the
                     // bottom (entry from the cytoplasm) residue order runs bottom-to-top.
                     const nAtBottom = pos.nEndY > pos.cEndY;
                     if (ssResidues.length === 0) return null;
                     const targetCodes = isBeta ? ['E', 'B'] : ['H', 'G', 'I'];
                     const mismatched = ssResidues.filter(r => !targetCodes.includes(r.code.toUpperCase()));
                     if (mismatched.length === 0) return null;
                     
                     // Draw bands for mismatched regions
                     const bands = [];
                     let startIdx = -1;
                     for (let i = 0; i < ssResidues.length; i++) {
                       const isMismatch = !targetCodes.includes(ssResidues[i].code.toUpperCase());
                       if (isMismatch && startIdx === -1) startIdx = i;
                       if (!isMismatch && startIdx !== -1) {
                         bands.push({ start: startIdx, end: i - 1 });
                         startIdx = -1;
                       }
                     }
                     if (startIdx !== -1) bands.push({ start: startIdx, end: ssResidues.length - 1 });
                     
                     return bands.map((band, idx) => {
                       const f0 = band.start / ssResidues.length;
                       const f1 = (band.end + 1) / ssResidues.length;
                       const top = nAtBottom ? 1 - f1 : f0;
                       const yStart = 6 + top * (cylHeight - 12);
                       const yHeight = Math.max(2, (f1 - f0) * (cylHeight - 12));
                       return (
                         <rect
                           key={`mismatch-${idx}`}
                           x="0"
                           y={yStart}
                           width={helixWidth}
                           height={yHeight}
                           fill="url(#hatch-warning)"
                           style={{ pointerEvents: 'none' }}
                         />
                       );
                     });
                  })()}
                  
                  <text
                    x={cx}
                    y={cy - 2}
                    className="helix-label-text"
                    textAnchor="middle"
                    transform={`rotate(${-pos.angle}, ${cx}, ${cy})`}
                    style={{ fill: labelColor }}
                  >
                    {h.subLabel}
                  </text>
                  <text
                    x={cx}
                    y={cy + 12}
                    className="helix-sub-text"
                    textAnchor="middle"
                    transform={`rotate(${-pos.angle}, ${cx}, ${cy})`}
                    style={{ fill: labelColor, opacity: 0.85 }}
                  >
                    {h.startRes}–{h.endRes}
                  </text>
                </g>
              );
            })}
          </g>

        </svg>
      </div>

      {/* Inspector */}
      <div className="tm-tooltip-bar">
        {hoveredElement ? (
          <div className="tm-tooltip-active">
            <strong className="tm-tooltip-title">{hoveredElement.title}</strong>
            <span>{hoveredElement.range}</span>
            <span className="tm-tooltip-len">{hoveredElement.length} residues</span>
            {hoveredElement.details && <small>{hoveredElement.details}</small>}
          </div>
        ) : (
          <div className="tm-tooltip-placeholder">
            <span>Hover a helix or loop to inspect it; click to jump to that residue in the 3D view.</span>
          </div>
        )}
      </div>
    </div>
  );
}