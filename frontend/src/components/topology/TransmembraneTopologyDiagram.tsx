import { useCallback, useEffect, useMemo, useRef, useState } from 'react';
import type { Chain } from '../../types/protein';
import type { SecondaryStructureResult, UniProtTopologyData } from '../../types/secondaryStructure';
import { exportSvgAsImage } from './exportDiagram';
import './TransmembraneTopologyDiagram.css';

type FigureTheme = 'publication' | 'lab';

function getNTermInside(uniprotData: UniProtTopologyData | null): boolean {
  if (!uniprotData?.regions?.length) return true;
  const firstDomain = uniprotData.regions
    .filter((r) => r.type === 'Topological domain')
    .sort((a, b) => a.start - b.start)[0];
  if (!firstDomain) return true;
  return !firstDomain.description.toLowerCase().includes('extracellular');
}

function classifyLoopType(domainDesc: string, loopIndex: number, nTermInside: boolean): 'EL' | 'IL' {
  const desc = domainDesc.toLowerCase();
  if (desc.includes('extracellular')) return 'EL';
  if (desc.includes('cytoplasmic') || desc.includes('intracellular')) return 'IL';
  const startWithEL = nTermInside;
  return loopIndex % 2 === 0 ? (startWithEL ? 'EL' : 'IL') : startWithEL ? 'IL' : 'EL';
}


interface Props {
  chain: Chain | undefined;
  secondaryResult?: SecondaryStructureResult | null;
  selectedResidue?: number | null;
  onSelectResidue?: (residueNumber: number) => void;
  uniprotId?: string | null;
}

export interface TMHelix {
  id: string;
  helixNumber: number;
  subLabel: string; // e.g. "1a", "1b", "2", "6a"
  startRes: number;
  endRes: number;
  length: number;
  color: string;
  isSplit?: boolean;
  partIndex?: number; // 0 for a, 1 for b
  description?: string;
}

export interface TMLoop {
  id: string;
  type: 'EL' | 'IL'; // Extracellular or Intracellular
  label: string; // e.g. "EL1", "IL1"
  startRes: number;
  endRes: number;
  length: number;
  prevHelixId: string;
  nextHelixId: string;
  hasShortHelix?: boolean;
  shortHelices?: { label: string; offsetFactor?: number }[];
  domainName?: string;
}

export interface CustomResidueColorRule {
  id: string;
  label: string;
  startRes: number;
  endRes: number;
  color: string;
}

// Preset Palettes
const PALETTES: Record<string, { label: string; colors: string[] }> = {
  PAPER_DEFAULT: {
    label: 'Figure (A) Paper Default',
    colors: [
      '#ef4444', // Red (TM1)
      '#f97316', // Orange (TM2)
      '#eab308', // Yellow (TM3)
      '#84cc16', // Lime (TM4)
      '#10b981', // Emerald Green (TM5)
      '#06b6d4', // Cyan (TM6)
      '#2563eb', // Royal Blue (TM7)
      '#6366f1', // Indigo (TM8)
      '#9333ea', // Purple (TM9)
      '#ec4899', // Pink (TM10)
      '#64748b', // Slate (TM11)
      '#334155', // Charcoal (TM12)
    ],
  },
  RAINBOW: {
    label: 'Rainbow Spectrum',
    colors: ['#ef4444', '#f97316', '#eab308', '#10b981', '#06b6d4', '#3b82f6', '#8b5cf6', '#d946ef', '#f43f5e', '#fb923c', '#a3e635', '#38bdf8'],
  },
  HYDROPHOBIC: {
    label: 'Hydrophobicity & Charge',
    colors: ['#3b82f6', '#60a5fa', '#93c5fd', '#f97316', '#fb923c', '#fdba74', '#10b981', '#34d399', '#6ee7b7', '#ec4899', '#f472b6', '#fbcfe8'],
  },
  NEON: {
    label: 'Vivid Neon',
    colors: ['#ff0055', '#ff6600', '#ffcc00', '#00ff66', '#00ffff', '#0099ff', '#9900ff', '#ff00ff', '#ff3399', '#00ffaa', '#ffff00', '#ff3300'],
  },
  PASTEL: {
    label: 'Pastel Bio',
    colors: ['#fca5a5', '#fdba74', '#fde047', '#bef264', '#86efac', '#67e8f9', '#93c5fd', '#c084fc', '#f472b6', '#fda4af', '#cbd5e1', '#94a3b8'],
  },
};

const DEFAULT_PRESETS = [
  { id: 'P31645', label: 'P31645 (hSERT)' },
  { id: 'P23975', label: 'P23975 (hNET)' },
  { id: 'P08183', label: 'P08183 (P-gp)' },
  { id: 'P00533', label: 'P00533 (EGFR)' },
];

export function TransmembraneTopologyDiagram({ chain, secondaryResult, selectedResidue, onSelectResidue, uniprotId }: Props) {
  const [uniprotIdInput, setUniprotIdInput] = useState<string>('P31645');
  const [uniprotData, setUniprotData] = useState<UniProtTopologyData | null>(null);
  const [loadingUniProt, setLoadingUniProt] = useState<boolean>(false);
  const [uniprotError, setUniprotError] = useState<string | null>(null);

  // Color Customizer state
  const [colorDrawerOpen, setColorDrawerOpen] = useState<boolean>(false);
  const [activeTab, setActiveTab] = useState<'preset' | 'helices' | 'residues'>('preset');
  const [selectedPaletteKey, setSelectedPaletteKey] = useState<string>('PAPER_DEFAULT');
  const [customHelixColors, setCustomHelixColors] = useState<Record<string, string>>({});
  const [customResidueRules, setCustomResidueRules] = useState<CustomResidueColorRule[]>([]);

  // Residue range input state
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

  // Fetch UniProt topology data from backend endpoint
  const fetchUniProtTopology = async (uniprotIdToFetch: string) => {
    if (!uniprotIdToFetch.trim()) return;
    setLoadingUniProt(true);
    setUniprotError(null);
    try {
      const response = await fetch(`http://localhost:8000/api/secondary-structure/uniprot/${uniprotIdToFetch.trim()}`);
      if (!response.ok) {
        const errData = await response.json().catch(() => ({}));
        throw new Error(errData.detail || `Failed to fetch UniProt data for ${uniprotIdToFetch}`);
      }
      const data: UniProtTopologyData = await response.json();
      setUniprotData(data);
    } catch (err: any) {
      setUniprotError(err.message || 'Error fetching UniProt data');
      setUniprotData(null);
    } finally {
      setLoadingUniProt(false);
    }
  };

  useEffect(() => {
    if (uniprotId) {
      setUniprotIdInput(uniprotId);
      fetchUniProtTopology(uniprotId);
    } else {
      fetchUniProtTopology('P31645');
    }
  }, [uniprotId]);

  // Compute active base colors array from preset or custom
  const currentPaletteColors = PALETTES[selectedPaletteKey]?.colors || PALETTES.PAPER_DEFAULT.colors;

  const getHelixColor = (subLabel: string, hNum: number): string => {
    if (customHelixColors[subLabel]) return customHelixColors[subLabel];
    if (customHelixColors[`TM${hNum}`]) return customHelixColors[`TM${hNum}`];
    return currentPaletteColors[(hNum - 1) % currentPaletteColors.length];
  };

  // Compute helices and loops
  const { helices, loops } = useMemo(() => {
    if (!chain) return { helices: [], loops: [] };

    // --- CASE A: UniProt Topology Data ---
    if (uniprotData && uniprotData.regions && uniprotData.regions.length > 0) {
      const tmRegions = uniprotData.regions.filter((r) => r.type === 'Transmembrane');
      const domainRegions = uniprotData.regions.filter((r) => r.type === 'Topological domain');
      const intraRegions = uniprotData.regions.filter((r) => r.type === 'Intramembrane');

      const tmHelices: TMHelix[] = [];
      let hNum = 1;

      tmRegions.forEach((tm) => {
        const len = tm.end - tm.start + 1;
        const rawName = tm.name || `${hNum}`;

        // Split helices: TM1 (1a/1b), TM6 (6a/6b), TM12 (12a/12b) or len >= 32
        const isSplitTarget = hNum === 1 || hNum === 6 || hNum === 12 || len >= 32;

        if (isSplitTarget && !rawName.includes('a') && !rawName.includes('b')) {
          const mid = Math.floor((tm.start + tm.end) / 2);
          tmHelices.push({
            id: `tm-${hNum}a`,
            helixNumber: hNum,
            subLabel: `${hNum}a`,
            startRes: tm.start,
            endRes: mid,
            length: mid - tm.start + 1,
            color: getHelixColor(`${hNum}a`, hNum),
            isSplit: true,
            partIndex: 0,
            description: tm.description,
          });
          tmHelices.push({
            id: `tm-${hNum}b`,
            helixNumber: hNum,
            subLabel: `${hNum}b`,
            startRes: mid + 1,
            endRes: tm.end,
            length: tm.end - mid,
            color: getHelixColor(`${hNum}b`, hNum),
            isSplit: true,
            partIndex: 1,
            description: tm.description,
          });
        } else {
          const subL = rawName.includes('a') || rawName.includes('b') ? rawName : `${hNum}`;
          tmHelices.push({
            id: `tm-${subL}`,
            helixNumber: hNum,
            subLabel: subL,
            startRes: tm.start,
            endRes: tm.end,
            length: len,
            color: getHelixColor(subL, hNum),
            isSplit: subL.endsWith('a') || subL.endsWith('b'),
            description: tm.description,
          });
        }
        hNum++;
      });

      // Build loops — alternate EL/IL based on UniProt topological domain
      const tmLoops: TMLoop[] = [];
      let elCount = 1;
      let ilCount = 1;
      let loopIndex = 0;
      const nTermInside = getNTermInside(uniprotData);

      for (let i = 0; i < tmHelices.length - 1; i++) {
        const hCurr = tmHelices[i];
        const hNext = tmHelices[i + 1];

        if (hCurr.isSplit && hNext.isSplit && hCurr.helixNumber === hNext.helixNumber) {
          continue;
        }

        const lStart = hCurr.endRes + 1;
        const lEnd = hNext.startRes - 1;
        const lLen = Math.max(0, lEnd - lStart + 1);
        const midRes = Math.floor((lStart + lEnd) / 2);

        const matchedDomain = domainRegions.find((d) => midRes >= d.start && midRes <= d.end);
        const domainDesc = matchedDomain ? matchedDomain.description : '';
        const loopSide = classifyLoopType(domainDesc, loopIndex, nTermInside);
        loopIndex++;

        const label = loopSide === 'EL' ? `EL${elCount++}` : `IL${ilCount++}`;

        // Mini helices definition based on biological paper Figure A:
        // EL2: "EL2", EL3: "3a", "3b", EL4: "4a", "4b", EL6: "EL6", IL1: "IL1", IL5: "IL5"
        let shortHelices: { label: string; offsetFactor?: number }[] | undefined = undefined;

        if (label === 'EL2') {
          shortHelices = [{ label: 'EL2' }];
        } else if (label === 'EL3') {
          shortHelices = [{ label: '3a', offsetFactor: -0.25 }, { label: '3b', offsetFactor: 0.25 }];
        } else if (label === 'EL4') {
          shortHelices = [{ label: '4a', offsetFactor: -0.25 }, { label: '4b', offsetFactor: 0.25 }];
        } else if (label === 'EL6') {
          shortHelices = [{ label: 'EL6' }];
        } else if (label === 'IL1') {
          shortHelices = [{ label: 'IL1' }];
        } else if (label === 'IL5') {
          shortHelices = [{ label: 'IL5' }];
        } else if (lLen > 16) {
          shortHelices = [{ label: label }];
        }

        tmLoops.push({
          id: `loop-${label}`,
          type: loopSide,
          label,
          startRes: lStart,
          endRes: lEnd,
          length: lLen,
          prevHelixId: hCurr.id,
          nextHelixId: hNext.id,
          hasShortHelix: !!shortHelices && shortHelices.length > 0,
          shortHelices,
          domainName: domainDesc,
        });
      }

      return { helices: tmHelices, loops: tmLoops };
    }

    // --- CASE B: DSSP / STRIDE Fallback ---
    const ssResidues = secondaryResult?.residues?.filter((r) => r.chain_id === chain.id) ?? [];
    const rawHelices: { start: number; end: number; code: string }[] = [];
    let currentHelix: { start: number; end: number; code: string } | null = null;

    if (ssResidues.length) {
      for (const r of ssResidues) {
        const code = r.code.toUpperCase();
        const isHelix = ['H', 'G', 'I'].includes(code);

        if (isHelix) {
          if (!currentHelix) {
            currentHelix = { start: r.residue_number, end: r.residue_number, code };
          } else {
            currentHelix.end = r.residue_number;
          }
        } else if (currentHelix) {
          if (currentHelix.end - currentHelix.start + 1 >= 5) {
            rawHelices.push(currentHelix);
          }
          currentHelix = null;
        }
      }
      if (currentHelix && currentHelix.end - currentHelix.start + 1 >= 5) {
        rawHelices.push(currentHelix);
      }
    } else {
      let curr = 15;
      const total = chain.residue_count;
      while (curr + 18 <= total) {
        rawHelices.push({ start: curr, end: Math.min(curr + 22, total - 10), code: 'H' });
        curr += 40;
      }
    }

    const tmHelices: TMHelix[] = [];
    let hNum = 1;

    rawHelices.forEach((h) => {
      const len = h.end - h.start + 1;
      const isSplitTarget = hNum === 1 || hNum === 6 || hNum === 12 || len >= 32;

      if (isSplitTarget) {
        const mid = Math.floor((h.start + h.end) / 2);
        tmHelices.push({
          id: `tm-${hNum}a`,
          helixNumber: hNum,
          subLabel: `${hNum}a`,
          startRes: h.start,
          endRes: mid,
          length: mid - h.start + 1,
          color: getHelixColor(`${hNum}a`, hNum),
          isSplit: true,
          partIndex: 0,
        });
        tmHelices.push({
          id: `tm-${hNum}b`,
          helixNumber: hNum,
          subLabel: `${hNum}b`,
          startRes: mid + 1,
          endRes: h.end,
          length: h.end - mid,
          color: getHelixColor(`${hNum}b`, hNum),
          isSplit: true,
          partIndex: 1,
        });
      } else {
        tmHelices.push({
          id: `tm-${hNum}`,
          helixNumber: hNum,
          subLabel: `${hNum}`,
          startRes: h.start,
          endRes: h.end,
          length: len,
          color: getHelixColor(`${hNum}`, hNum),
          isSplit: false,
        });
      }
      hNum++;
    });

    const tmLoops: TMLoop[] = [];
    let elCount = 1;
    let ilCount = 1;
    let loopIndex = 0;

    for (let i = 0; i < tmHelices.length - 1; i++) {
      const hCurr = tmHelices[i];
      const hNext = tmHelices[i + 1];

      if (hCurr.isSplit && hNext.isSplit && hCurr.helixNumber === hNext.helixNumber) {
        continue;
      }

      const lStart = hCurr.endRes + 1;
      const lEnd = hNext.startRes - 1;
      const lLen = Math.max(0, lEnd - lStart + 1);

      const loopSide = classifyLoopType('', loopIndex, true);
      loopIndex++;
      const label = loopSide === 'EL' ? `EL${elCount++}` : `IL${ilCount++}`;

      let shortHelices: { label: string; offsetFactor?: number }[] | undefined = undefined;
      if (label === 'EL2') shortHelices = [{ label: 'EL2' }];
      else if (label === 'EL3') shortHelices = [{ label: '3a', offsetFactor: -0.25 }, { label: '3b', offsetFactor: 0.25 }];
      else if (label === 'EL4') shortHelices = [{ label: '4a', offsetFactor: -0.25 }, { label: '4b', offsetFactor: 0.25 }];
      else if (label === 'EL6') shortHelices = [{ label: 'EL6' }];
      else if (label === 'IL1') shortHelices = [{ label: 'IL1' }];
      else if (label === 'IL5') shortHelices = [{ label: 'IL5' }];

      tmLoops.push({
        id: `loop-${label}`,
        type: loopSide,
        label,
        startRes: lStart,
        endRes: lEnd,
        length: lLen,
        prevHelixId: hCurr.id,
        nextHelixId: hNext.id,
        hasShortHelix: !!shortHelices,
        shortHelices,
      });
    }

    return { helices: tmHelices, loops: tmLoops };
  }, [chain, secondaryResult, uniprotData, selectedPaletteKey, customHelixColors]);

  const handleExport = useCallback(
    async (format: 'png' | 'jpeg') => {
      if (!svgRef.current) return;
      setExporting(true);
      try {
        const id = uniprotData?.uniprot_id ?? 'topology';
        const ext = format === 'jpeg' ? 'jpg' : 'png';
        await exportSvgAsImage(svgRef.current, format, `${id}_TM_topology.${ext}`, 3);
      } catch {
        /* export failed silently — user can retry */
      } finally {
        setExporting(false);
      }
    },
    [uniprotData?.uniprot_id],
  );

  // Handle adding custom residue color
  const handleAddResidueRule = () => {
    const start = parseInt(resStartInput.trim(), 10);
    const end = resEndInput.trim() ? parseInt(resEndInput.trim(), 10) : start;
    if (isNaN(start)) return;

    const newRule: CustomResidueColorRule = {
      id: `rule-${Date.now()}`,
      label: start === end ? `Residue ${start}` : `Residues ${start}–${end}`,
      startRes: Math.min(start, end),
      endRes: Math.max(start, end),
      color: resColorInput,
    };

    setCustomResidueRules((prev) => [...prev, newRule]);
    setResStartInput('');
    setResEndInput('');
  };

  const handleRemoveResidueRule = (id: string) => {
    setCustomResidueRules((prev) => prev.filter((r) => r.id !== id));
  };

  const handleResetColors = () => {
    setSelectedPaletteKey('PAPER_DEFAULT');
    setCustomHelixColors({});
    setCustomResidueRules([]);
  };

  // Check if a helix or residue range has a custom residue color assigned
  const getCustomResidueColorForHelix = (startRes: number, endRes: number): string | null => {
    const matched = customResidueRules.find((rule) => rule.startRes <= endRes && rule.endRes >= startRes);
    return matched ? matched.color : null;
  };

  if (!chain) {
    return <div className="empty-state">Upload a structure to view the 2D Transmembrane Secondary Structure Map.</div>;
  }

  const isPub = figureTheme === 'publication';

  // Layout geometry calculations
  const helixWidth = 44;
  const colSpacing = 68;
  const canvasWidth = Math.max(960, helices.length * colSpacing + 160);
  const membraneTopY = 155;
  const membraneBottomY = 255;
  const membraneHeight = membraneBottomY - membraneTopY;

  const helixPositions: Record<string, { x: number; topY: number; bottomY: number; angle: number }> = {};
  let currentX = 100;

  helices.forEach((h) => {
    const isEvenHelix = h.helixNumber % 2 === 0;
    const tiltAngle = isEvenHelix ? -5 : 5;

    if (h.isSplit) {
      if (h.partIndex === 0) {
        // Part A: Intracellular side or lower half
        helixPositions[h.id] = {
          x: currentX,
          topY: h.subLabel.endsWith('a') && (h.helixNumber === 1 || h.helixNumber === 6) ? membraneTopY + membraneHeight / 2 + 5 : membraneTopY - 14,
          bottomY: h.subLabel.endsWith('a') && (h.helixNumber === 1 || h.helixNumber === 6) ? membraneBottomY + 16 : membraneTopY + membraneHeight / 2 - 5,
          angle: tiltAngle + 3,
        };
      } else {
        // Part B: Extracellular side or upper half
        helixPositions[h.id] = {
          x: currentX,
          topY: h.subLabel.endsWith('b') && (h.helixNumber === 1 || h.helixNumber === 6) ? membraneTopY - 16 : membraneTopY + membraneHeight / 2 + 5,
          bottomY: h.subLabel.endsWith('b') && (h.helixNumber === 1 || h.helixNumber === 6) ? membraneTopY + membraneHeight / 2 - 5 : membraneBottomY + 16,
          angle: tiltAngle - 3,
        };
        currentX += colSpacing;
      }
    } else {
      helixPositions[h.id] = {
        x: currentX,
        topY: membraneTopY - 20,
        bottomY: membraneBottomY + 20,
        angle: tiltAngle,
      };
      currentX += colSpacing;
    }
  });

  return (
    <div className="tm-topology-panel panel">
      {/* Panel Header & Toolbar */}
      <div className="panel-heading" style={{ flexWrap: 'wrap', gap: '12px' }}>
        <div>
          <span className="section-kicker">UNIPROT TOPOLOGY · 2D STRUCTURE MAP</span>
          <h2>Transmembrane Secondary Structure Map</h2>
          <p className="panel-subtitle">
            {uniprotData
              ? `${uniprotData.protein_name}${uniprotData.gene_name ? ` (${uniprotData.gene_name})` : ''} · ${uniprotData.uniprot_id}`
              : 'Transmembrane helices, extracellular/intracellular loops, and re-entrant segments.'}
          </p>
        </div>

        <div className="tm-toolbar-actions">
          {/* Theme Selector */}
          <div className="tm-preset-chip" style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
            <button
              className={`tm-tab-btn ${figureTheme === 'publication' ? 'active' : ''}`}
              onClick={() => setFigureTheme('publication')}
            >
              Publication (Paper)
            </button>
            <button className={`tm-tab-btn ${figureTheme === 'lab' ? 'active' : ''}`} onClick={() => setFigureTheme('lab')}>
              Dark Lab
            </button>
          </div>

          {/* Export PNG/JPEG Buttons */}
          <button className="tm-color-toggle-btn" onClick={() => handleExport('png')} disabled={exporting}>
            <span>📷</span>
            <span>{exporting ? 'Exporting...' : 'Export PNG'}</span>
          </button>

          {/* Color Customizer Drawer Toggle Button */}
          <button className={`tm-color-toggle-btn ${colorDrawerOpen ? 'active' : ''}`} onClick={() => setColorDrawerOpen((open) => !open)}>
            <span>🎨</span>
            <span>{colorDrawerOpen ? 'Close Color Picker' : 'Customize Colors'}</span>
          </button>

          {/* UniProt Search Input */}
          <div style={{ display: 'flex', alignItems: 'center', gap: '6px' }}>
            <input
              type="text"
              value={uniprotIdInput}
              onChange={(e) => setUniprotIdInput(e.target.value)}
              placeholder="e.g. P31645"
              style={{
                width: '90px',
                padding: '4px 8px',
                borderRadius: '5px',
                border: '1px solid #334155',
                background: '#0f172a',
                color: '#f8fafc',
                fontSize: '12px',
                fontWeight: 700,
                textTransform: 'uppercase',
              }}
              onKeyDown={(e) => {
                if (e.key === 'Enter') fetchUniProtTopology(uniprotIdInput);
              }}
            />
            <button
              onClick={() => fetchUniProtTopology(uniprotIdInput)}
              disabled={loadingUniProt}
              style={{
                padding: '4px 10px',
                borderRadius: '5px',
                border: 'none',
                background: '#3b82f6',
                color: '#ffffff',
                fontSize: '11px',
                fontWeight: 700,
                cursor: 'pointer',
              }}
            >
              {loadingUniProt ? '...' : 'Fetch'}
            </button>
          </div>
        </div>
      </div>

      {/* Quick Presets row */}
      <div style={{ display: 'flex', alignItems: 'center', justifyContent: 'space-between', gap: '8px', margin: '4px 0 12px', flexWrap: 'wrap' }}>
        <div style={{ display: 'flex', gap: '6px', alignItems: 'center' }}>
          <small style={{ fontSize: '11px', color: '#94a3b8', fontWeight: 600 }}>UniProt Presets:</small>
          {DEFAULT_PRESETS.map((preset) => (
            <button
              key={preset.id}
              onClick={() => {
                setUniprotIdInput(preset.id);
                fetchUniProtTopology(preset.id);
              }}
              style={{
                padding: '2px 8px',
                borderRadius: '4px',
                border: '1px solid #334155',
                background: uniprotData?.uniprot_id === preset.id ? '#1e293b' : '#090d16',
                color: uniprotData?.uniprot_id === preset.id ? '#38bdf8' : '#94a3b8',
                fontSize: '10px',
                cursor: 'pointer',
                fontWeight: 600,
              }}
            >
              {preset.label}
            </button>
          ))}
        </div>

        <span className="tag">
          {helices.length} HELICES · {loops.length} LOOPS
        </span>
      </div>

      {/* COLOR CUSTOMIZER DRAWER */}
      {colorDrawerOpen && (
        <div className="tm-color-customizer-drawer">
          <div className="tm-drawer-header">
            <h4>
              <span>🖌️ Residue & Helix Color Customizer</span>
            </h4>
            <div className="tm-drawer-tabs">
              <button className={`tm-tab-btn ${activeTab === 'preset' ? 'active' : ''}`} onClick={() => setActiveTab('preset')}>
                Presets & Themes
              </button>
              <button className={`tm-tab-btn ${activeTab === 'helices' ? 'active' : ''}`} onClick={() => setActiveTab('helices')}>
                Helix Colors
              </button>
              <button className={`tm-tab-btn ${activeTab === 'residues' ? 'active' : ''}`} onClick={() => setActiveTab('residues')}>
                Custom Residues ({customResidueRules.length})
              </button>
            </div>
          </div>

          {/* TAB 1: PRESETS */}
          {activeTab === 'preset' && (
            <div className="tm-preset-row">
              <label>Choose Color Theme:</label>
              {Object.keys(PALETTES).map((key) => (
                <button
                  key={key}
                  className={`tm-preset-chip ${selectedPaletteKey === key ? 'active' : ''}`}
                  onClick={() => {
                    setSelectedPaletteKey(key);
                    setCustomHelixColors({});
                  }}
                >
                  {PALETTES[key].label}
                </button>
              ))}
              <button className="tm-reset-btn" onClick={handleResetColors} style={{ marginLeft: 'auto' }}>
                Reset All Colors
              </button>
            </div>
          )}

          {/* TAB 2: HELIX COLORS */}
          {activeTab === 'helices' && (
            <div>
              <div className="tm-helix-color-grid">
                {helices.map((h) => (
                  <div key={h.id} className="tm-helix-color-item">
                    <span>{`Helix ${h.subLabel}`}</span>
                    <input
                      type="color"
                      value={h.color}
                      className="tm-color-input"
                      onChange={(e) => {
                        setCustomHelixColors((prev) => ({ ...prev, [h.subLabel]: e.target.value }));
                      }}
                    />
                  </div>
                ))}
              </div>
            </div>
          )}

          {/* TAB 3: CUSTOM RESIDUE / RANGE COLORS */}
          {activeTab === 'residues' && (
            <div className="tm-residue-color-section">
              <div className="tm-residue-form">
                <label style={{ fontSize: '11px', color: '#94a3b8' }}>Residue # or Range:</label>
                <input
                  type="number"
                  placeholder="Start (e.g. 45)"
                  value={resStartInput}
                  onChange={(e) => setResStartInput(e.target.value)}
                  className="tm-input-field"
                />
                <span style={{ color: '#64748b' }}>to</span>
                <input
                  type="number"
                  placeholder="End (e.g. 60)"
                  value={resEndInput}
                  onChange={(e) => setResEndInput(e.target.value)}
                  className="tm-input-field"
                />
                <input type="color" value={resColorInput} className="tm-color-input" onChange={(e) => setResColorInput(e.target.value)} />
                <button className="tm-add-btn" onClick={handleAddResidueRule}>
                  + Add Residue Color
                </button>

                {selectedResidue && (
                  <button
                    className="tm-add-btn"
                    style={{ background: '#3b82f6' }}
                    onClick={() => {
                      setResStartInput(selectedResidue.toString());
                      setResEndInput(selectedResidue.toString());
                    }}
                  >
                    Color Selected Residue ({selectedResidue})
                  </button>
                )}
              </div>

              {/* Active Rules List */}
              {customResidueRules.length > 0 ? (
                <div className="tm-residue-tag-list">
                  {customResidueRules.map((rule) => (
                    <span key={rule.id} className="tm-residue-tag">
                      <span style={{ width: '10px', height: '10px', borderRadius: '50%', background: rule.color, display: 'inline-block' }} />
                      <strong>{rule.label}</strong>
                      <button onClick={() => handleRemoveResidueRule(rule.id)}>✕</button>
                    </span>
                  ))}
                </div>
              ) : (
                <small style={{ color: '#64748b', fontSize: '11px' }}>
                  No custom residue colors added yet. Enter a residue number or range above to highlight specific positions on the map.
                </small>
              )}
            </div>
          )}
        </div>
      )}

      {uniprotError && (
        <div className="tm-error-banner">{uniprotError}</div>
      )}

      {/* 2D SVG DIAGRAM */}
      <div className={`tm-diagram-wrap ${isPub ? 'publication' : 'lab'}`}>
        <svg
          ref={svgRef}
          className="tm-diagram-svg"
          viewBox={`0 0 ${canvasWidth} 440`}
          role="img"
          aria-label="Transmembrane 2D Topology Diagram"
        >
          {isPub && <rect x="0" y="0" width={canvasWidth} height="440" fill="#ffffff" />}
          <defs>
            {/* 3D Opaque Cylinder Linear Gradients for Helices */}
            {helices.map((h) => {
              const resColor = getCustomResidueColorForHelix(h.startRes, h.endRes);
              const mainColor = resColor || h.color;
              return (
                <linearGradient key={`grad-${h.id}`} id={`cyl-grad-${h.id}`} x1="0%" y1="0%" x2="100%" y2="0%">
                  <stop offset="0%" stopColor={mainColor} stopOpacity="1" />
                  <stop offset="32%" stopColor="#ffffff" stopOpacity="0.6" />
                  <stop offset="68%" stopColor={mainColor} stopOpacity="1" />
                  <stop offset="100%" stopColor={isPub ? '#1e293b' : '#0b131b'} stopOpacity="1" />
                </linearGradient>
              );
            })}

            {/* Short Helices Metallic Gradient */}
            <linearGradient id="short-helix-grad" x1="0%" y1="0%" x2="100%" y2="0%">
              <stop offset="0%" stopColor="#64748b" stopOpacity="1" />
              <stop offset="40%" stopColor="#ffffff" stopOpacity="0.65" />
              <stop offset="75%" stopColor="#475569" stopOpacity="1" />
              <stop offset="100%" stopColor="#1e293b" stopOpacity="1" />
            </linearGradient>

            {/* Membrane Lipid Bilayer Gradient */}
            <linearGradient id="membrane-grad" x1="0%" y1="0%" x2="0%" y2="100%">
              <stop offset="0%" stopColor="#334155" stopOpacity="0.35" />
              <stop offset="50%" stopColor="#1e293b" stopOpacity="0.2" />
              <stop offset="100%" stopColor="#334155" stopOpacity="0.35" />
            </linearGradient>
          </defs>

          {isPub && (
            <text x="28" y="32" className="panel-letter">(A)</text>
          )}

          {/* Central Lipid Bilayer Band */}
          <g className="membrane-zone">
            <rect
              x="20"
              y={membraneTopY}
              width={canvasWidth - 40}
              height={membraneHeight}
              fill={isPub ? '#e2e8f0' : 'url(#membrane-grad)'}
              rx="4"
              stroke={isPub ? '#94a3b8' : 'none'}
              strokeWidth={isPub ? 1 : 0}
            />
            <line
              x1="20"
              y1={membraneTopY}
              x2={canvasWidth - 20}
              y2={membraneTopY}
              stroke={isPub ? '#64748b' : '#64748b'}
              strokeWidth={isPub ? 1.5 : 2}
              strokeDasharray={isPub ? 'none' : '6 4'}
              opacity="0.85"
            />
            <line
              x1="20"
              y1={membraneBottomY}
              x2={canvasWidth - 20}
              y2={membraneBottomY}
              stroke={isPub ? '#64748b' : '#64748b'}
              strokeWidth={isPub ? 1.5 : 2}
              strokeDasharray={isPub ? 'none' : '6 4'}
              opacity="0.85"
            />
            <text x="18" y={membraneTopY - 10} className="membrane-label" textAnchor="start">
              Extracellular
            </text>
            <text x="18" y={membraneBottomY + 22} className="membrane-label" textAnchor="start">
              Cytoplasmic
            </text>
          </g>

          {/* Loop Curves & Short Helices */}
          <g className="tm-loops">
            {loops.map((loop) => {
              const posPrev = helixPositions[loop.prevHelixId];
              const posNext = helixPositions[loop.nextHelixId];
              if (!posPrev || !posNext) return null;

              const isEL = loop.type === 'EL';
              const startX = posPrev.x + helixWidth / 2;
              const startY = isEL ? posPrev.topY : posPrev.bottomY;
              const endX = posNext.x + helixWidth / 2;
              const endY = isEL ? posNext.topY : posNext.bottomY;

              const midX = (startX + endX) / 2;
              const archDepth = 30 + Math.min(38, loop.length * 0.8);
              const apexY = isEL
                ? Math.max(48, Math.min(startY, endY) - archDepth)
                : Math.min(360, Math.max(startY, endY) + archDepth);

              const pathD = `M ${startX} ${startY} C ${startX} ${apexY}, ${endX} ${apexY}, ${endX} ${endY}`;
              const loopStroke = isPub ? '#334155' : (helices.find((h) => h.id === loop.prevHelixId)?.color ?? '#64748b');

              // Position loop title badge vs short helix mini-cylinders so they never overlap
              const badgeY = isEL
                ? (loop.hasShortHelix ? apexY - 30 : apexY - 20)
                : (loop.hasShortHelix ? apexY + 18 : apexY + 8);
              const badgeTextY = isEL
                ? (loop.hasShortHelix ? apexY - 18 : apexY - 8)
                : (loop.hasShortHelix ? apexY + 30 : apexY + 20);

              return (
                <g
                  key={loop.id}
                  className="loop-group"
                  onMouseEnter={() =>
                    setHoveredElement({
                      title: `${loop.label} (${loop.type === 'EL' ? 'Extracellular' : 'Intracellular'} Loop)`,
                      range: `Residues ${loop.startRes}–${loop.endRes}`,
                      length: loop.length,
                      details: `Connects ${loop.prevHelixId} → ${loop.nextHelixId}${loop.domainName ? ` · ${loop.domainName}` : ''}`,
                    })
                  }
                  onMouseLeave={() => setHoveredElement(null)}
                  onClick={() => onSelectResidue?.(loop.startRes)}
                >
                  <path d={pathD} fill="none" stroke={loopStroke} strokeWidth={isPub ? 2.2 : 2.8} strokeLinecap="round" opacity={isPub ? 0.85 : 0.9} />

                  {/* Loop title badge pill */}
                  <rect
                    x={midX - 19}
                    y={badgeY}
                    width="38"
                    height="16"
                    rx="8"
                    fill={isPub ? '#ffffff' : '#1e293b'}
                    stroke={isPub ? '#94a3b8' : '#475569'}
                    strokeWidth="1.2"
                    opacity="0.98"
                  />
                  <text x={midX} y={badgeTextY} className="loop-text-label" textAnchor="middle">
                    {loop.label}
                  </text>

                  {/* Short loop helices (re-entrant mini-cylinders) */}
                  {loop.shortHelices &&
                    loop.shortHelices.map((sh, idx) => {
                      const offset = (sh.offsetFactor || 0) * 42;
                      const shX = midX + offset - 18;
                      const shY = apexY - 6;

                      return (
                        <g key={`${loop.id}-sh-${idx}`} transform={`translate(${shX}, ${shY})`}>
                          <rect x="0" y="0" width="36" height="15" rx="7" fill="url(#short-helix-grad)" stroke={isPub ? '#475569' : '#ffffff'} strokeWidth="1" />
                          <text x="18" y="11" className="short-helix-text" textAnchor="middle">
                            {sh.label}
                          </text>
                        </g>
                      );
                    })}
                </g>
              );
            })}
          </g>

          {/* Transmembrane 3D Cylinders */}
          <g className="tm-helices">
            {helices.map((h) => {
              const pos = helixPositions[h.id];
              if (!pos) return null;

              const cylHeight = Math.abs(pos.bottomY - pos.topY);
              const isSelected = selectedResidue !== null && selectedResidue !== undefined && selectedResidue >= h.startRes && selectedResidue <= h.endRes;
              const resCustomColor = getCustomResidueColorForHelix(h.startRes, h.endRes);

              return (
                <g
                  key={h.id}
                  className={`helix-group ${isSelected ? 'selected' : ''}`}
                  transform={`translate(${pos.x}, ${pos.topY}) rotate(${pos.angle}, ${helixWidth / 2}, ${cylHeight / 2})`}
                  onMouseEnter={() =>
                    setHoveredElement({
                      title: `Helix ${h.subLabel} (TM${h.helixNumber})`,
                      range: `Residues ${h.startRes}–${h.endRes}`,
                      length: h.length,
                      details: h.description ? `UniProt: ${h.description}` : `Secondary Structure: Transmembrane Alpha Helix`,
                    })
                  }
                  onMouseLeave={() => setHoveredElement(null)}
                  onClick={() => onSelectResidue?.(h.startRes)}
                >
                  {/* 3D Opaque Cylinder Body */}
                  <rect
                    x="0"
                    y="0"
                    width={helixWidth}
                    height={cylHeight}
                    rx="14"
                    fill={`url(#cyl-grad-${h.id})`}
                    stroke={isSelected ? '#ff6f61' : resCustomColor ? resCustomColor : (isPub ? '#334155' : '#0f172a')}
                    strokeWidth={isSelected || resCustomColor ? '2.5' : '1.2'}
                    className="cylinder-body"
                  />

                  {/* Top Cap */}
                  <ellipse
                    cx={helixWidth / 2}
                    cy="7"
                    rx={helixWidth / 2 - 1}
                    ry="5"
                    fill={resCustomColor || h.color}
                    opacity="1.0"
                    stroke="#ffffff"
                    strokeWidth="1"
                  />

                  {/* Helix Label Inside Cylinder */}
                  <text
                    x={helixWidth / 2}
                    y={cylHeight <= 60 ? cylHeight / 2 - 3 : cylHeight / 2 - 4}
                    className="helix-label-text"
                    textAnchor="middle"
                  >
                    {h.subLabel}
                  </text>

                  {/* Residue Range Text */}
                  <text
                    x={helixWidth / 2}
                    y={cylHeight <= 60 ? cylHeight / 2 + 11 : cylHeight / 2 + 13}
                    className="helix-sub-text"
                    textAnchor="middle"
                  >
                    {h.startRes}-{h.endRes}
                  </text>

                  {/* Custom Residue Highlight Indicator */}
                  {resCustomColor && (
                    <circle cx={helixWidth / 2} cy="16" r="4" fill={resCustomColor} stroke="#ffffff" strokeWidth="1" />
                  )}
                </g>
              );
            })}
          </g>

          {/* Termini (N & C) */}
          {helices.length > 0 && (
            <g className="terminals">
              <text x={(helixPositions[helices[0].id]?.x ?? 40) - 20} y={membraneBottomY + 45} className="terminal-text" textAnchor="middle">
                N
              </text>
              <g transform={`translate(${(helixPositions[helices[helices.length - 1].id]?.x ?? canvasWidth - 60) + 40}, ${membraneBottomY + 30})`}>
                <rect x="0" y="0" width="36" height="15" rx="7" fill="url(#short-helix-grad)" stroke={isPub ? '#475569' : '#ffffff'} strokeWidth="1" />
                <text x="18" y="11" className="short-helix-text" textAnchor="middle">
                  C-term
                </text>
                <text x="48" y="14" className="terminal-text" textAnchor="middle">
                  C
                </text>
              </g>
            </g>
          )}

          {/* Bottom 1D domain track */}
          <g className="bottom-domain-track" transform="translate(40, 395)">
            <line x1="20" y1="12" x2={canvasWidth - 100} y2="12" stroke={isPub ? '#cbd5e1' : '#334155'} strokeWidth="2.5" />
            <text x="5" y="16" className="domain-track-label">
              N
            </text>

            {helices.map((h) => {
              const pos = helixPositions[h.id];
              if (!pos) return null;
              const resColor = getCustomResidueColorForHelix(h.startRes, h.endRes);
              const barColor = resColor || h.color;

              return (
                <g key={`track-${h.id}`} transform={`translate(${pos.x - 20}, 3)`} onClick={() => onSelectResidue?.(h.startRes)}>
                  <rect x="0" y="0" width={h.isSplit ? 26 : 38} height="18" rx="4" fill={barColor} stroke="#0f172a" strokeWidth="1" />
                  <text x={h.isSplit ? 13 : 19} y="13" className="domain-bar-text" textAnchor="middle">
                    {h.subLabel}
                  </text>
                </g>
              );
            })}
            <text x={canvasWidth - 85} y="16" className="domain-track-label">
              C
            </text>
          </g>
        </svg>
      </div>

      {/* Inspector / Tooltip Bar */}
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
            <span>Hover or click any transmembrane helix cylinder, loop curve, or 1D track block to inspect coordinates & navigate.</span>
          </div>
        )}
      </div>
    </div>
  );
}
