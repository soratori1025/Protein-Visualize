import React, { useCallback, useEffect, useMemo, useRef, useState } from 'react';
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
import { ConsensusAnalysisMap } from './ConsensusAnalysisMap';

export type FigureTheme = 'publication' | 'lab';
export type TopologySource = 'uniprot' | 'calculated';
export type MembraneSide = 'out' | 'in'; // out = extracellular/lumenal, in = cytoplasmic

/* ------------------------------------------------------------------ *
 * Color utilities
 * ------------------------------------------------------------------ */

export function hexToRgb(hex: string): [number, number, number] {
  const h = hex.replace('#', '').trim();
  const full = h.length === 3 ? h.split('').map((c) => c + c).join('') : h;
  return [
    parseInt(full.slice(0, 2), 16) || 0,
    parseInt(full.slice(2, 4), 16) || 0,
    parseInt(full.slice(4, 6), 16) || 0,
  ];
}

export function rgbToHex(r: number, g: number, b: number): string {
  const c = (v: number) => Math.max(0, Math.min(255, Math.round(v))).toString(16).padStart(2, '0');
  return `#${c(r)}${c(g)}${c(b)}`;
}

/** Linear blend between two hex colors. amount = 0 → a, 1 → b. */
export function mixHex(a: string, b: string, amount: number): string {
  const [r1, g1, b1] = hexToRgb(a);
  const [r2, g2, b2] = hexToRgb(b);
  const t = Math.max(0, Math.min(1, amount));
  return rgbToHex(r1 + (r2 - r1) * t, g1 + (g2 - g1) * t, b1 + (b2 - b1) * t);
}

export const lighten = (hex: string, amount: number) => mixHex(hex, '#ffffff', amount);
export const darken = (hex: string, amount: number) => mixHex(hex, '#000000', amount);

export const renderRibbon = (width: number, hSpan: number, baseColor: string, id: string, yStart: number) => {
  const pitch = 18; // Smaller pitch for more turns
  const turns = Math.max(1, Math.floor(hSpan / pitch));
  const actualPitch = hSpan / turns;
  
  const backFaces = [];
  const frontFaces = [];
  const thickness = width * 0.22; // Thinner ribbon!
  
  // To keep the rounded caps strictly inside the visual bounds (optional but looks cleaner)
  const pad = thickness / 2;
  const w = width - pad * 2;
  const xOffset = pad;
  
  for (let i = 0; i < turns; i++) {
    const y0 = yStart + i * actualPitch;
    const yMid = y0 + actualPitch * 0.5;
    const yEnd = y0 + actualPitch;
    
    // Zero-derivative at the edges creates a mathematically perfect 2D projection of a 3D helix
    backFaces.push(
      <path key={`back-${i}`} d={`M ${xOffset + w} ${y0} C ${xOffset + w * 0.5} ${y0}, ${xOffset + w * 0.5} ${yMid}, ${xOffset} ${yMid}`} 
            fill="none" stroke={darken(baseColor, 0.55)} strokeWidth={thickness} strokeLinecap="round" />
    );
    
    frontFaces.push(
      <path key={`front-${i}`} d={`M ${xOffset} ${yMid} C ${xOffset + w * 0.5} ${yMid}, ${xOffset + w * 0.5} ${yEnd}, ${xOffset + w} ${yEnd}`} 
            fill="none" stroke={`url(#cyl-${id})`} strokeWidth={thickness} strokeLinecap="round" />
    );
  }
  
  return (
    <React.Fragment key={`ribbon-${yStart}`}>
      <g className="ribbon-back">{backFaces}</g>
      <g className="ribbon-front">{frontFaces}</g>
    </React.Fragment>
  );
};

/** WCAG relative luminance — decides whether a label should be white or near-black. */
export function relativeLuminance(hex: string): number {
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
export function getContrastTextColor(bgHex: string): string {
  return relativeLuminance(bgHex) > 0.45 ? '#0f172a' : '#ffffff';
}

/**
 * Resample a base palette to exactly `count` colors by interpolating between
 * neighbouring stops. Works for 3 TMs or 24 TMs — no modulo wrap-around, so two
 * helices never share the same color.
 */
export function buildPalette(base: string[], count: number): string[] {
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
export const OUT_WORDS = [
  'extracellular', 'lumenal', 'luminal', 'periplasm', 'exoplasmic', 'vesicular',
  'intermembrane', 'perinuclear', 'vacuolar', 'peroxisomal', 'thylakoid', 'virion surface',
];
export const IN_WORDS = ['cytoplasm', 'cytosol', 'intracellular', 'mitochondrial matrix', 'matrix', 'stromal', 'intravirion'];

/** Map a UniProt topological-domain description to a membrane side, or null if unknown. */
export function descriptionToSide(description?: string): MembraneSide | null {
  const d = (description ?? '').toLowerCase();
  if (!d) return null;
  if (OUT_WORDS.some((w) => d.includes(w))) return 'out';
  if (IN_WORDS.some((w) => d.includes(w)) || d.trim() === 'nuclear') return 'in';
  return null;
}

/** Side of a region: the backend's explicit `side` first, the description otherwise. */
export function regionSide(region: Pick<UniProtTopologyRegion, 'side' | 'description'>): MembraneSide | null {
  if (region.side === 'Cytoplasmic') return 'in';
  if (region.side === 'Extracellular') return 'out';
  return descriptionToSide(region.description);
}

export const opposite = (s: MembraneSide): MembraneSide => (s === 'out' ? 'in' : 'out');

/**
 * Side on which the FIRST crossing starts. Taken from the first topological domain
 * whose side is known, flipped once per crossing that lies before it. (The old
 * version used that domain's side directly, which draws the whole map upside down
 * whenever the first annotated domain comes after TM1 — e.g. a chain that starts
 * inside the membrane, or a UniProt entry without an N-terminal domain.)
 */
export function getNTermSide(data: UniProtTopologyData | null, crossingStarts: number[]): MembraneSide {
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

export interface RawTM {
  start: number;
  end: number;
  name?: string;
  description?: string;
  /** 'Helix' | 'Strand' | 'Irregular' from the calculated endpoint. */
  ss?: string | null;
  confidence?: string | null;
  /** Calculated endpoint: crossing number + 'a'/'b' part of a broken crossing. */
  crossing?: number | null;
  part?: string | null;
  observedSpans?: { start: number; end: number }[];
}

export function mentionsDiscontinuity(tm: RawTM): boolean {
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
export function isBetaStrandDescription(description?: string): boolean {
  const d = (description ?? '').toLowerCase();
  if (!d) return false;
  const hasBeta = d.includes('beta') || d.includes('β');
  return hasBeta && (d.includes('strand') || d.includes('barrel') || d.includes('sheet'));
}

/** β segment: explicit `ss` from the calculated endpoint, else the description. */
export function isBetaSegment(tm: { ss?: string | null; description?: string }): boolean {
  return tm.ss === 'Strand' || isBetaStrandDescription(tm.description);
}

/**
 * Two consecutive annotated segments that are really the two halves of one
 * broken helix: nearly touching, each too short to span the bilayer on its own,
 * or explicitly named "…a"/"…b" or "part 1"/"part 2".
 */
export function looksLikeHalfPair(prev: RawTM, curr: RawTM, trustSegments = false): boolean {
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

export interface Props {
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
  distinguishTurns?: boolean;
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
  /** Drawn dashed: irregular membrane segment or low-confidence call. */
  uncertain?: boolean;
  confidence?: string | null;
  observedSpans?: { start: number; end: number }[];
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

export const PALETTES: Record<string, { label: string; colors: string[] }> = {
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

export const DEFAULT_PRESETS = [
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
  chainId: string | undefined,
  distinguishTurns?: boolean
): ExtraFeature[] | undefined {
  if (lEnd < lStart) return undefined;

  const found: ExtraFeature[] = [];

  // 1. UniProt features
  for (const region of topologyData?.regions ?? []) {
    if (
      ['Intramembrane', 'Signal', 'Domain', 'Region', 'Motif', 'Propeptide', 'Repeat', 'Turn'].includes(region.type) &&
      region.start >= lStart &&
      region.end <= lEnd
    ) {
      found.push({
        startRes: region.start,
        endRes: region.end,
        type: region.type,
        label: region.topology_label || region.name || region.description || region.type,
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
        // backend topology label (EL2, EL3a, IL1 …) names the element by its loop,
        // like published topology figures; SS type stays in `type`
        label:
          region.topology_label ||
          (region.description.includes('Interfacial') ? 'interf. α' : ss === 'Helix' ? 'α-helix' : 'β-strand'),
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

  if (topologyData?.consensus_map && distinguishTurns) {
    let turnStart = null;
    for (let i = lStart; i <= lEnd; i++) {
       const res = topologyData.consensus_map.find(r => r.residue_number === i);
       if (res && res.label && res.label.startsWith('Turn_')) {
          if (turnStart === null) turnStart = i;
       } else {
          if (turnStart !== null) {
             found.push({
                label: 'Turn',
                startRes: turnStart,
                endRes: i - 1,
                type: 'turn'
             });
             turnStart = null;
          }
       }
    }
    if (turnStart !== null) {
       found.push({
          label: 'Turn',
          startRes: turnStart,
          endRes: lEnd,
          type: 'turn'
       });
    }
  }

  const filteredFound = !(isCalculatedTopology(topologyData) && topologyData.residues && topologyData.membrane)
    ? found
    : found.filter((feature) => {
        // Temporarily hide RE and IL/EL segments that are outside the lipid membrane
        const isReOrLoopHelix = feature.type === 'Intramembrane' || feature.type === 'Helix' || feature.type === 'Strand';
        if (isReOrLoopHelix) {
          const insideCount = topologyData.residues!.filter(
            (res) => res.residue_number >= feature.startRes && res.residue_number <= feature.endRes && (res.zone === 'CORE' || res.zone === 'EDGE')
          ).length;
          return insideCount >= 5;
        }
        return true;
      });

  if (filteredFound.length === 0) return undefined;

  // Keep the longest 3 to avoid cluttering
  const kept = filteredFound
    .sort((a, b) => b.endRes - b.startRes - (a.endRes - a.startRes))
    .slice(0, 3)
    .sort((a, b) => a.startRes - b.startRes);

  return kept.map((sh, idx) => ({
    ...sh,
    offsetFactor: kept.length > 1 ? idx - (kept.length - 1) / 2 : 0,
  }));
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
  distinguishTurns?: boolean;
}

export function buildTopologyModel({
  chain,
  activeTopologyData,
  secondaryResult = null,
  structureSS = null,
  topologySource,
  selectedPaletteKey = 'PAPER_DEFAULT',
  customHelixColors = {},
  distinguishTurns = false,
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

    if (isCalculatedTopology(activeTopologyData) && activeTopologyData.consensus_map) {
      // Synchronize strictly with the Consensus Merge Analysis map
      const tmInMap = new Map<number, typeof activeTopologyData.consensus_map>();
      activeTopologyData.consensus_map.forEach((r) => {
        if (r.label === 'TM_in' && r.tm_segment != null) {
          if (!tmInMap.has(r.tm_segment)) tmInMap.set(r.tm_segment, []);
          tmInMap.get(r.tm_segment)!.push(r);
        }
      });

      const segments: (typeof activeTopologyData.consensus_map)[] = [];
      const observedSpansMap = new Map<number, { start: number; end: number }[]>();

      for (const residues of tmInMap.values()) {
        if (residues.length === 0) continue;
        
        let currentSubSeg = [residues[0]];
        const observedSpans: { start: number; end: number }[] = [];
        
        for (let i = 1; i < residues.length; i++) {
          if (residues[i].index === residues[i - 1].index + 1) {
            currentSubSeg.push(residues[i]);
          } else {
            observedSpans.push({
              start: currentSubSeg[0].residue_number,
              end: currentSubSeg[currentSubSeg.length - 1].residue_number,
            });
            currentSubSeg = [residues[i]];
          }
        }
        observedSpans.push({
          start: currentSubSeg[0].residue_number,
          end: currentSubSeg[currentSubSeg.length - 1].residue_number,
        });

        segments.push([...residues]);
        observedSpansMap.set(residues[0].index, observedSpans);
      }

      segments.sort((a, b) => a[0].index - b[0].index);
      
      const crossingCounts = new Map<number, number>();
      for (const seg of segments) {
        if (seg[0].crossing != null) {
          crossingCounts.set(seg[0].crossing, (crossingCounts.get(seg[0].crossing) || 0) + 1);
        }
      }

      const crossingSeen = new Map<number, number>();
      rawTMs = segments.map((seg) => {
        const start = seg[0].residue_number;
        const end = seg[seg.length - 1].residue_number;
        const crossing = seg[0].crossing;
        const observedSpans = observedSpansMap.get(seg[0].index);
        
        let part: string | null = null;
        if (crossing != null && crossingCounts.get(crossing)! > 1) {
          const seen = crossingSeen.get(crossing) || 0;
          part = String.fromCharCode(97 + seen); // 'a', 'b', 'c', etc.
          crossingSeen.set(crossing, seen + 1);
        }

        return {
          start,
          end,
          name: part ? `TM${crossing || ''}${part}` : `TM${crossing || ''}`,
          description: part ? 'Transmembrane Unwound' : 'Transmembrane', // Helps shouldSplit
          ss: seg[0].ss_raw === 'E' ? 'Strand' : 'Helix',
          confidence: 'high',
          crossing,
          part,
          observedSpans,
        };
      });

      domainRegions = activeTopologyData.regions
        .filter((r) => r.type === 'Topological domain')
        .map((r) => ({ start: r.start, end: r.end, description: r.description ?? '' }));
    } else {
      rawTMs = (activeTopologyData.regions ?? [])
        // the unwound stretch between the two halves of a broken crossing is drawn as
        // the connector between them, not as a helix of its own
        .filter((r) => r.type === 'Transmembrane' && !(r.description ?? '').includes('Unwound'))
        .map((r) => ({
          start: r.start,
          end: r.end,
          name: r.name,
          description: r.description,
          ss: r.ss,
          confidence: r.confidence,
          crossing: r.crossing,
          part: r.part,
        }))
        .sort((a, b) => a.start - b.start);
      domainRegions = activeTopologyData.regions
        .filter((r) => r.type === 'Topological domain')
        .map((r) => ({ start: r.start, end: r.end, description: r.description ?? '' }));
    }
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
  const hasCrossingIds = trustSegments && rawTMs.some((tm) => tm.crossing != null);
  for (const tm of rawTMs) {
    const last = groups[groups.length - 1];
    if (hasCrossingIds) {
      // The backend decided which pieces form one crossing (broken helix = same
      // crossing number, parts 'a' and 'b'); draw exactly that.
      if (last && tm.crossing != null && last[0].crossing === tm.crossing) last.push(tm);
      else groups.push([tm]);
    } else if (last && last.length === 1 && looksLikeHalfPair(last[0], tm, trustSegments)) {
      last.push(tm);
    } else {
      groups.push([tm]);
    }
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
    const uncertain = first.ss === 'Irregular' || first.confidence === 'low';

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
        uncertain,
        confidence: first.confidence,
        entrySide,
        exitSide,
        column: gi,
        observedSpans: first.observedSpans,
      });
      return;
    }

    const halves = preSplit
      ? [
          { start: first.start, end: first.end, description: first.description, observedSpans: first.observedSpans },
          { start: last.start, end: last.end, description: last.description, observedSpans: last.observedSpans },
        ]
      : (() => {
          const mid = Math.floor((first.start + last.end) / 2);
          return [
            { start: first.start, end: mid, description: first.description, observedSpans: first.observedSpans },
            { start: mid + 1, end: last.end, description: first.description, observedSpans: first.observedSpans },
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
        uncertain,
        confidence: first.confidence,
        // Part a runs from the entry side inwards; part b carries on to the exit side.
        entrySide: idx === 0 ? entrySide : entrySide,
        exitSide: idx === 0 ? entrySide : exitSide,
        column: gi,
        observedSpans: half.observedSpans,
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
      chain.id,
      distinguishTurns
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

