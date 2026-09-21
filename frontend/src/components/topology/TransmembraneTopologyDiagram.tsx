import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import type { Chain } from '../../types/protein';
import type { SecondaryStructureResult, UniProtTopologyData } from '../../types/secondaryStructure';
import { exportSvgAsImage } from './exportDiagram';
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

const OUT_WORDS = ['extracellular', 'lumenal', 'luminal', 'periplasmic', 'exoplasmic', 'vesicular'];
const IN_WORDS = ['cytoplasmic', 'intracellular', 'cytosolic', 'matrix', 'stromal', 'nuclear'];

/** Map a UniProt topological-domain description to a membrane side, or null if unknown. */
function descriptionToSide(description?: string): MembraneSide | null {
  const d = (description ?? '').toLowerCase();
  if (!d) return null;
  if (OUT_WORDS.some((w) => d.includes(w))) return 'out';
  if (IN_WORDS.some((w) => d.includes(w))) return 'in';
  return null;
}

const opposite = (s: MembraneSide): MembraneSide => (s === 'out' ? 'in' : 'out');

/** Side of the N-terminus, read from the first annotated topological domain. */
function getNTermSide(data: UniProtTopologyData | null): MembraneSide {
  const first = (data?.regions ?? [])
    .filter((r) => r.type === 'Topological domain')
    .sort((a, b) => a.start - b.start)[0];
  return descriptionToSide(first?.description) ?? 'in';
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

/**
 * Two consecutive annotated segments that are really the two halves of one
 * broken helix: nearly touching, each too short to span the bilayer on its own,
 * or explicitly named "…a"/"…b" or "part 1"/"part 2".
 */
function looksLikeHalfPair(prev: RawTM, curr: RawTM): boolean {
  const gap = curr.start - prev.end - 1;
  if (gap < 0 || gap > 7) return false;

  const prevText = `${prev.name ?? ''} ${prev.description ?? ''}`.toLowerCase();
  const currText = `${curr.name ?? ''} ${curr.description ?? ''}`.toLowerCase();

  const prevNum = prevText.match(/(\d+)\s*a\b/);
  const currNum = currText.match(/(\d+)\s*b\b/);
  if (prevNum && currNum && prevNum[1] === currNum[1]) return true;
  if (prevText.includes('part 1') && currText.includes('part 2')) return true;

  // β-strands are legitimately short and pack close together in a barrel — never
  // fuse two of them into "the two halves of one broken helix" (that heuristic is
  // only meant for a single α-helix that UniProt split across the bilayer).
  if (isBetaStrandDescription(prev.description) || isBetaStrandDescription(curr.description)) return false;

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
      .sort((a, b) => a.residue_number - b.residue_number);

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
}: Props) {
  const [uniprotIdInput, setUniprotIdInput] = useState<string>('P31645');
  const [uniprotData, setUniprotData] = useState<UniProtTopologyData | null>(null);
  const [loadingUniProt, setLoadingUniProt] = useState<boolean>(false);
  const [uniprotError, setUniprotError] = useState<string | null>(null);

  const [topologySource, setTopologySource] = useState<TopologySource>('uniprot');
  const [calcAlgorithm, setCalcAlgorithm] = useState<
    'dssp_ss' | 'stride_ss' | 'dssp_slab' | 'stride_slab' | 'kd_slab' | 'tmhmm_seq'
  >('dssp_ss');
  const [overlayType, setOverlayType] = useState<'none' | 'dssp' | 'stride'>('none');
  const [calculatedData, setCalculatedData] = useState<UniProtTopologyData | null>(null);
  const [loadingCalculated, setLoadingCalculated] = useState<boolean>(false);
  const [calculatedError, setCalculatedError] = useState<string | null>(null);
  const [showUniProtInfo, setShowUniProtInfo] = useState<boolean>(false);

  // --- Advanced TM parameters (user-tunable biological thresholds) ---
  const [showAdvancedParams, setShowAdvancedParams] = useState(false);
  const [tmThickness, setTmThickness] = useState<string>('30');
  const [tmMinElement, setTmMinElement] = useState<string>('4');
  const [tmMinCrossSpan, setTmMinCrossSpan] = useState<string>('0.45');
  const [tmFullCrossFrac, setTmFullCrossFrac] = useState<string>('0.66');
  const [tmBrokenGapMax, setTmBrokenGapMax] = useState<string>('9');
  const [tmMinMembraneScore, setTmMinMembraneScore] = useState<string>('0.5');

  const [colorDrawerOpen, setColorDrawerOpen] = useState<boolean>(false);
  const [activeTab, setActiveTab] = useState<'preset' | 'helices' | 'residues'>('preset');
  const [selectedPaletteKey, setSelectedPaletteKey] = useState<string>('PAPER_DEFAULT');
  const [customHelixColors, setCustomHelixColors] = useState<Record<string, string>>({});
  const [customResidueRules, setCustomResidueRules] = useState<CustomResidueColorRule[]>([]);

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
      const response = await fetch(`http://localhost:8000/api/secondary-structure/uniprot/${uniprotIdToFetch.trim()}`);
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

  const fetchCalculatedTopology = useCallback(async (filenameToFetch: string, algorithm: string) => {
    if (!filenameToFetch.trim()) return;
    setLoadingCalculated(true);
    setCalculatedError(null);
    try {
      // Build query string with algorithm + any non-default advanced params
      const qp = new URLSearchParams({ algorithm });
      const th = parseFloat(tmThickness);    if (!isNaN(th) && th !== 30)    qp.set('thickness', String(th));
      const me = parseInt(tmMinElement);      if (!isNaN(me) && me !== 4)     qp.set('min_tm_element', String(me));
      const cs = parseFloat(tmMinCrossSpan);  if (!isNaN(cs) && cs !== 0.45)  qp.set('min_cross_span', String(cs));
      const fc = parseFloat(tmFullCrossFrac); if (!isNaN(fc) && fc !== 0.66)  qp.set('full_cross_frac', String(fc));
      const bg = parseInt(tmBrokenGapMax);    if (!isNaN(bg) && bg !== 9)     qp.set('broken_gap_max', String(bg));
      const ms = parseFloat(tmMinMembraneScore); if (!isNaN(ms) && ms !== 0.5) qp.set('min_membrane_score', String(ms));

      const response = await fetch(
        `http://localhost:8000/api/secondary-structure/predict-topology/${encodeURIComponent(filenameToFetch.trim())}?${qp.toString()}`
      );
      if (!response.ok) {
        const errData = await response.json().catch(() => ({}));
        throw new Error(errData.detail || `Could not compute topology for ${filenameToFetch}`);
      }
      setCalculatedData((await response.json()) as UniProtTopologyData);
    } catch (err: any) {
      setCalculatedError(err.message || 'Could not compute topology');
      setCalculatedData(null);
    } finally {
      setLoadingCalculated(false);
    }
  }, [tmThickness, tmMinElement, tmMinCrossSpan, tmFullCrossFrac, tmBrokenGapMax, tmMinMembraneScore]);

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
    setCalculatedData(null);
    setCalculatedError(null);
    if (topologySource === 'calculated' && filename) fetchCalculatedTopology(filename, calcAlgorithm);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [filename, calcAlgorithm]);

  useEffect(() => {
    if (topologySource === 'calculated' && filename && !calculatedData && !loadingCalculated) {
      fetchCalculatedTopology(filename, calcAlgorithm);
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [topologySource]);

  const activeTopologyData = topologySource === 'calculated' ? calculatedData : uniprotData;
  const activeError = topologySource === 'calculated' ? calculatedError : uniprotError;
  const activeLoading = topologySource === 'calculated' ? loadingCalculated : loadingUniProt;

  /* ---------------------------------------------------------------- *
   * Build helices + loops
   * ---------------------------------------------------------------- */

  const { helices, loops, helixCount } = useMemo(() => {
    const empty = { helices: [] as TMHelix[], loops: [] as TMLoop[], helixCount: 0 };
    if (!chain) return empty;

    /* --- Step 1: collect raw transmembrane segments --- */
    let rawTMs: RawTM[] = [];
    let domainRegions: { start: number; end: number; description: string }[] = [];
    let usingAnnotation = false;

    if (activeTopologyData?.regions?.length) {
      usingAnnotation = true;
      rawTMs = activeTopologyData.regions
        .filter((r) => r.type === 'Transmembrane')
        .map((r) => ({ start: r.start, end: r.end, name: r.name, description: r.description }))
        .sort((a, b) => a.start - b.start);
      domainRegions = activeTopologyData.regions
        .filter((r) => r.type === 'Topological domain')
        .map((r) => ({ start: r.start, end: r.end, description: r.description ?? '' }));
    } else {
      // Fallback: derive candidate TM segments from the assigned secondary structure.
      const ssResidues = (secondaryResult?.residues ?? [])
        .filter((r) => r.chain_id === chain.id)
        .sort((a, b) => a.residue_number - b.residue_number);

      let run: RawTM | null = null;
      const flush = () => {
        if (run && run.end - run.start + 1 >= 15) rawTMs.push(run);
        run = null;
      };
      for (const r of ssResidues) {
        if (['H', 'G', 'I'].includes(r.code.toUpperCase())) {
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
    const groups: RawTM[][] = [];
    for (const tm of rawTMs) {
      const last = groups[groups.length - 1];
      if (last && last.length === 1 && looksLikeHalfPair(last[0], tm)) last.push(tm);
      else groups.push([tm]);
    }

    /* --- Step 3: palette sized to the actual number of helices --- */
    const base = PALETTES[selectedPaletteKey]?.colors ?? PALETTES.PAPER_DEFAULT.colors;
    const palette = buildPalette(base, groups.length);
    const colorFor = (subLabel: string, hNum: number) =>
      customHelixColors[subLabel] ?? customHelixColors[`TM${hNum}`] ?? palette[hNum - 1] ?? base[0];

    /* --- Step 4: which side does each helix start on? --- */
    const nTermSide = usingAnnotation ? getNTermSide(activeTopologyData) : 'in';

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
      const shouldSplit = preSplit || mentionsDiscontinuity(first) || spanLen >= LONG_TM_RESIDUES;

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
        secondaryResult,
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
  }, [chain, secondaryResult, activeTopologyData, selectedPaletteKey, customHelixColors]);

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
    <div className="tm-topology-panel panel">
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
            {activeTopologyData
              ? `${activeTopologyData.protein_name}${
                  activeTopologyData.gene_name ? ` (${activeTopologyData.gene_name})` : ''
                } · ${activeTopologyData.uniprot_id} · ${helixCount} helices, ${loops.length} loops`
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
              onClick={() => setTopologySource('uniprot')}
            >
              UniProt
            </button>
            <button
              className={`tm-tab-btn ${topologySource === 'calculated' ? 'active' : ''}`}
              onClick={() => setTopologySource('calculated')}
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
            onChange={(e) => setOverlayType(e.target.value as any)}
            style={{ padding: '4px 8px', width: 'auto' }}
            title="Secondary Structure Overlay"
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
              <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                <select 
                  className="tm-input-field" 
                  value={calcAlgorithm}
                  onChange={(e) => setCalcAlgorithm(e.target.value as any)}
                  style={{ padding: '4px 8px' }}
                >
                  <option value="dssp_ss">DSSP elements + membrane (recommended)</option>
                  <option value="stride_ss">STRIDE elements + membrane</option>
                  <option value="dssp_slab">DSSP + slab (legacy)</option>
                  <option value="stride_slab">STRIDE + slab (legacy)</option>
                  <option value="kd_slab">Geometry only (no SS)</option>
                  <option value="tmhmm_seq">Sequence only (no 3D)</option>
                </select>
                <button
                  onClick={() => filename && fetchCalculatedTopology(filename, calcAlgorithm)}
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
              </div>
              {showAdvancedParams && (
                <div style={{
                  margin: '8px 0', padding: '10px 14px',
                  background: 'rgba(100,100,140,0.08)', borderRadius: '8px',
                  display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '8px 16px',
                  fontSize: '0.82em',
                }}>
                  <label title="Hydrophobic core thickness (Å). Bacterial IM ~27, eukaryotic PM ~30, ER ~25. (Mitra 2004; OPM database)">
                    Membrane thickness (Å)
                    <input type="number" step="0.5" min="20" max="40" value={tmThickness}
                      onChange={e => setTmThickness(e.target.value)}
                      className="tm-input-field" style={{ width: '70px', marginLeft: 4 }} />
                  </label>
                  <label title="Min residues of a helix/strand inside the slab to count as TM. Lower = more sensitive.">
                    Min element in slab
                    <input type="number" step="1" min="2" max="15" value={tmMinElement}
                      onChange={e => setTmMinElement(e.target.value)}
                      className="tm-input-field" style={{ width: '50px', marginLeft: 4 }} />
                  </label>
                  <label title="Min fraction of thickness a crossing must span. Lower admits shallower crossings.">
                    Min cross span
                    <input type="number" step="0.05" min="0.2" max="0.9" value={tmMinCrossSpan}
                      onChange={e => setTmMinCrossSpan(e.target.value)}
                      className="tm-input-field" style={{ width: '60px', marginLeft: 4 }} />
                  </label>
                  <label title="Fraction of thickness above which a single element is a full crossing (no fusion).">
                    Full cross frac
                    <input type="number" step="0.05" min="0.3" max="1.0" value={tmFullCrossFrac}
                      onChange={e => setTmFullCrossFrac(e.target.value)}
                      className="tm-input-field" style={{ width: '60px', marginLeft: 4 }} />
                  </label>
                  <label title="Max gap (residues) to fuse two partial helices into one crossing (broken helix detection).">
                    Broken gap max
                    <input type="number" step="1" min="3" max="20" value={tmBrokenGapMax}
                      onChange={e => setTmBrokenGapMax(e.target.value)}
                      className="tm-input-field" style={{ width: '50px', marginLeft: 4 }} />
                  </label>
                  <label title="Min mean hydrophobicity inside the slab. Below this the protein is treated as soluble.">
                    Min membrane score
                    <input type="number" step="0.1" min="0" max="3" value={tmMinMembraneScore}
                      onChange={e => setTmMinMembraneScore(e.target.value)}
                      className="tm-input-field" style={{ width: '60px', marginLeft: 4 }} />
                  </label>
                  <div style={{ gridColumn: '1 / -1', marginTop: 4, opacity: 0.65, fontSize: '0.9em' }}>
                    Defaults are literature-derived (Kyte–Doolittle 1982; Mitra 2004; OPM/PDBTM).
                    Hover each label for details.
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
            <button className="tm-reset-btn" onClick={handleResetColors}>
              Reset
            </button>
          </div>

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
                  onClick={() => onSelectResidue?.(loop.startRes)}
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

            const nFeatures = getExtraFeatures(1, first.startRes - 1, activeTopologyData, secondaryResult, chain?.id) || [];
            const lastResNum = chain?.sequence?.length || 10000;
            const cFeatures = getExtraFeatures(last.endRes + 1, lastResNum, activeTopologyData, secondaryResult, chain?.id) || [];

            return (
              <g className="tm-terminals">
                {/* N-terminus */}
                <path
                  d={`M 20 ${first.entrySide === 'out' ? membraneTopY - 60 : membraneBottomY + 60} Q ${posFirst.x / 2} ${first.entrySide === 'out' ? membraneTopY - 30 : membraneBottomY + 30} ${posFirst.x + helixWidth / 2} ${posFirst.nEndY}`}
                  fill="none"
                  stroke={first.color}
                  strokeWidth="2.6"
                />
                <text x="20" y={first.entrySide === 'out' ? membraneTopY - 70 : membraneBottomY + 75} textAnchor="middle" style={{ fill: isPub ? '#1e293b' : '#e2e8f0', fontSize: '11px', fontWeight: 'bold' }}>NH2</text>
                
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
                <text x={canvasWidth - 30} y={last.exitSide === 'out' ? membraneTopY - 70 : membraneBottomY + 75} textAnchor="middle" style={{ fill: isPub ? '#1e293b' : '#e2e8f0', fontSize: '11px', fontWeight: 'bold' }}>COOH</text>
                
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
                selectedResidue != null && selectedResidue >= h.startRes && selectedResidue <= h.endRes;
              const resCustomColor = getCustomResidueColorForHelix(h.startRes, h.endRes);
              const baseColor = resCustomColor ?? h.color;
              const cx = helixWidth / 2;
              const cy = cylHeight / 2;
              const labelColor = getContrastTextColor(baseColor);

              const isBeta = isBetaStrandDescription(h.description);
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
                  onClick={() => onSelectResidue?.(h.startRes)}
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
                  {overlayType !== 'none' && secondaryResult && (() => {
                     const ssResidues = secondaryResult.residues.filter(r => r.chain_id === chain?.id && r.residue_number >= h.startRes && r.residue_number <= h.endRes);
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
                       const yStart = 6 + (band.start / ssResidues.length) * (cylHeight - 12);
                       const yHeight = Math.max(2, ((band.end - band.start + 1) / ssResidues.length) * (cylHeight - 12));
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

          {/* N and C termini */}
          {firstHelix && lastHelix && (
            <g className="terminals">
              {(() => {
                const p = helixPositions[firstHelix.id];
                if (!p) return null;
                const y = p.nEndY;
                const down = firstHelix.entrySide === 'in';
                const tipY = down ? y + 42 : y - 42;
                return (
                  <>
                    <path
                      d={`M ${p.x + helixWidth / 2} ${y} C ${p.x - 4} ${(y + tipY) / 2}, ${p.x - 30} ${
                        (y + tipY) / 2
                      }, ${p.x - 34} ${tipY}`}
                      fill="none"
                      stroke={firstHelix.color}
                      strokeWidth="2.4"
                      strokeLinecap="round"
                    />
                    <text
                      x={p.x - 48}
                      y={tipY + 5}
                      className="terminal-text"
                      textAnchor="middle"
                      style={{ fill: isPub ? '#0f172a' : '#e2e8f0' }}
                    >
                      N
                    </text>
                  </>
                );
              })()}

              {(() => {
                const p = helixPositions[lastHelix.id];
                if (!p) return null;
                const y = p.cEndY;
                const down = lastHelix.exitSide === 'in';
                const tipY = down ? y + 42 : y - 42;
                return (
                  <>
                    <path
                      d={`M ${p.x + helixWidth / 2} ${y} C ${p.x + helixWidth + 6} ${(y + tipY) / 2}, ${
                        p.x + helixWidth + 26
                      } ${(y + tipY) / 2}, ${p.x + helixWidth + 32} ${tipY}`}
                      fill="none"
                      stroke={lastHelix.color}
                      strokeWidth="2.4"
                      strokeLinecap="round"
                    />
                    <text
                      x={p.x + helixWidth + 48}
                      y={tipY + 5}
                      className="terminal-text"
                      textAnchor="middle"
                      style={{ fill: isPub ? '#0f172a' : '#e2e8f0' }}
                    >
                      C
                    </text>
                  </>
                );
              })()}
            </g>
          )}
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