import { useMemo, useState } from 'react';
import type { SecondaryStructureResidue, SecondaryStructureResult } from '../../types/secondaryStructure';

const labels: Record<string, string> = {
  H: 'Alpha helix',
  G: '3-10 helix',
  I: 'Pi helix',
  E: 'Beta strand',
  B: 'Beta bridge',
  T: 'Turn',
  S: 'Bend',
  '-': 'Loop / coil',
};

const colors: Record<string, string> = {
  H: '#e76f51',
  G: '#f4a261',
  I: '#e9c46a',
  E: '#4ea8de',
  B: '#72b7d6',
  T: '#9b8afb',
  S: '#8aa17b',
  '-': '#496673',
};

interface Props {
  result: SecondaryStructureResult;
  selectedResidue?: number | null;
  onSelectResidue?: (residueNumber: number) => void;
}

interface SSEBlock {
  id: string;
  type: 'helix' | 'strand' | 'turn' | 'coil';
  code: string;
  label: string;
  start: number;
  end: number;
  length: number;
}

export function SecondaryStructureTrack({ result, selectedResidue, onSelectResidue }: Props) {
  const [hoveredResidue, setHoveredResidue] = useState<SecondaryStructureResidue | null>(null);
  const [activeFilter, setActiveFilter] = useState<'ALL' | 'HELIX' | 'STRAND' | 'TURN' | 'COIL'>('ALL');
  const [hoveredSSEBlock, setHoveredSSEBlock] = useState<SSEBlock | null>(null);

  const totalResidues = result.residues.length;

  // Category counts & percentages
  const categories = useMemo(() => {
    let helix = 0;
    let strand = 0;
    let turn = 0;
    let coil = 0;

    for (const r of result.residues) {
      const code = r.code.toUpperCase();
      if (['H', 'G', 'I'].includes(code)) helix++;
      else if (['E', 'B'].includes(code)) strand++;
      else if (['T', 'S'].includes(code)) turn++;
      else coil++;
    }

    return {
      helix: { count: helix, pct: totalResidues ? ((helix / totalResidues) * 100).toFixed(1) : '0' },
      strand: { count: strand, pct: totalResidues ? ((strand / totalResidues) * 100).toFixed(1) : '0' },
      turn: { count: turn, pct: totalResidues ? ((turn / totalResidues) * 100).toFixed(1) : '0' },
      coil: { count: coil, pct: totalResidues ? ((coil / totalResidues) * 100).toFixed(1) : '0' },
    };
  }, [result.residues, totalResidues]);

  // Group contiguous residues into Secondary Structure Elements (SSEs)
  const sseBlocks = useMemo(() => {
    const blocks: SSEBlock[] = [];
    if (!result.residues.length) return blocks;

    let current: { type: 'helix' | 'strand' | 'turn' | 'coil'; code: string; start: number; end: number } | null = null;
    let hCount = 1;
    let sCount = 1;
    let tCount = 1;

    for (const r of result.residues) {
      const code = r.code.toUpperCase();
      let type: 'helix' | 'strand' | 'turn' | 'coil' = 'coil';
      if (['H', 'G', 'I'].includes(code)) type = 'helix';
      else if (['E', 'B'].includes(code)) type = 'strand';
      else if (['T', 'S'].includes(code)) type = 'turn';

      if (!current) {
        current = { type, code, start: r.residue_number, end: r.residue_number };
      } else if (current.type === type) {
        current.end = r.residue_number;
      } else {
        const len = current.end - current.start + 1;
        let label = 'Loop';
        if (current.type === 'helix') label = `α${hCount++}`;
        else if (current.type === 'strand') label = `β${sCount++}`;
        else if (current.type === 'turn') label = `Turn ${tCount++}`;

        blocks.push({
          id: `${current.type}-${current.start}`,
          type: current.type,
          code: current.code,
          label,
          start: current.start,
          end: current.end,
          length: len,
        });

        current = { type, code, start: r.residue_number, end: r.residue_number };
      }
    }

    if (current) {
      const len = current.end - current.start + 1;
      let label = 'Loop';
      if (current.type === 'helix') label = `α${hCount++}`;
      else if (current.type === 'strand') label = `β${sCount++}`;
      else if (current.type === 'turn') label = `Turn ${tCount++}`;

      blocks.push({
        id: `${current.type}-${current.start}`,
        type: current.type,
        code: current.code,
        label,
        start: current.start,
        end: current.end,
        length: len,
      });
    }

    return blocks;
  }, [result.residues]);

  // Filtered SSE blocks
  const filteredBlocks = useMemo(() => {
    if (activeFilter === 'HELIX') return sseBlocks.filter((b) => b.type === 'helix');
    if (activeFilter === 'STRAND') return sseBlocks.filter((b) => b.type === 'strand');
    if (activeFilter === 'TURN') return sseBlocks.filter((b) => b.type === 'turn');
    if (activeFilter === 'COIL') return sseBlocks.filter((b) => b.type === 'coil');
    return sseBlocks;
  }, [sseBlocks, activeFilter]);

  // Generate ruler markers (every 50 residues or at start/end)
  const rulerMarkers = useMemo(() => {
    if (!result.residues.length) return [];
    const minRes = result.residues[0].residue_number;
    const maxRes = result.residues[result.residues.length - 1].residue_number;

    const markers: { number: number; index: number }[] = [];
    result.residues.forEach((r, index) => {
      if (r.residue_number === minRes || r.residue_number === maxRes || r.residue_number % 50 === 0) {
        markers.push({ number: r.residue_number, index });
      }
    });
    return markers;
  }, [result.residues]);

  const activeResidue = hoveredResidue ?? result.residues.find((r) => r.residue_number === selectedResidue);

  return (
    <div className="secondary-track-workspace">
      {/* Header Banner */}
      <div className="secondary-track-header">
        <div className="secondary-title-group">
          <span className="analysis-label">ASSIGNMENT MAP</span>
          <strong className="secondary-method-tag">{result.method} ASSIGNMENT</strong>
        </div>
        <div className="secondary-badge">
          <span>{totalResidues} RESIDUES</span>
        </div>
      </div>

      {/* Composition Proportional Bar */}
      <div className="secondary-composition-bar-wrap">
        <div className="secondary-composition-bar" role="img" aria-label="Secondary structure composition">
          <div className="comp-segment helix-seg" style={{ width: `${categories.helix.pct}%` }} title={`Helices: ${categories.helix.count} res (${categories.helix.pct}%)`} />
          <div className="comp-segment strand-seg" style={{ width: `${categories.strand.pct}%` }} title={`Strands: ${categories.strand.count} res (${categories.strand.pct}%)`} />
          <div className="comp-segment turn-seg" style={{ width: `${categories.turn.pct}%` }} title={`Turns: ${categories.turn.count} res (${categories.turn.pct}%)`} />
          <div className="comp-segment coil-seg" style={{ width: `${categories.coil.pct}%` }} title={`Coils: ${categories.coil.count} res (${categories.coil.pct}%)`} />
        </div>
      </div>

      {/* Category Summary Cards / Quick Filters */}
      <div className="secondary-category-cards">
        <button
          className={`category-card helix-card ${activeFilter === 'HELIX' ? 'active' : ''}`}
          onClick={() => setActiveFilter(activeFilter === 'HELIX' ? 'ALL' : 'HELIX')}
        >
          <div className="cat-icon-dot helix-dot" />
          <div className="cat-info">
            <span className="cat-name">Alpha Helix (α)</span>
            <strong className="cat-val">{categories.helix.count} <small>({categories.helix.pct}%)</small></strong>
          </div>
        </button>

        <button
          className={`category-card strand-card ${activeFilter === 'STRAND' ? 'active' : ''}`}
          onClick={() => setActiveFilter(activeFilter === 'STRAND' ? 'ALL' : 'STRAND')}
        >
          <div className="cat-icon-dot strand-dot" />
          <div className="cat-info">
            <span className="cat-name">Beta Strand (β)</span>
            <strong className="cat-val">{categories.strand.count} <small>({categories.strand.pct}%)</small></strong>
          </div>
        </button>

        <button
          className={`category-card turn-card ${activeFilter === 'TURN' ? 'active' : ''}`}
          onClick={() => setActiveFilter(activeFilter === 'TURN' ? 'ALL' : 'TURN')}
        >
          <div className="cat-icon-dot turn-dot" />
          <div className="cat-info">
            <span className="cat-name">Turn / Bend</span>
            <strong className="cat-val">{categories.turn.count} <small>({categories.turn.pct}%)</small></strong>
          </div>
        </button>

        <button
          className={`category-card coil-card ${activeFilter === 'COIL' ? 'active' : ''}`}
          onClick={() => setActiveFilter(activeFilter === 'COIL' ? 'ALL' : 'COIL')}
        >
          <div className="cat-icon-dot coil-dot" />
          <div className="cat-info">
            <span className="cat-name">Loop / Coil</span>
            <strong className="cat-val">{categories.coil.count} <small>({categories.coil.pct}%)</small></strong>
          </div>
        </button>
      </div>

      {/* Hover & Active Inspector Card */}
      <div className={`secondary-inspector-card ${activeResidue ? 'active' : ''}`}>
        {activeResidue ? (
          <div className="inspector-content">
            <div className="inspector-main">
              <span className="inspector-res-chip" style={{ backgroundColor: colors[activeResidue.code] ?? colors['-'] }}>
                {activeResidue.code}
              </span>
              <div>
                <strong className="inspector-res-title">
                  Chain {activeResidue.chain_id} : Residue {activeResidue.residue_number} · {activeResidue.residue_name}
                </strong>
                <span className="inspector-ss-label">{labels[activeResidue.code] ?? activeResidue.code}</span>
              </div>
            </div>
            <div className="inspector-angles">
              {activeResidue.phi !== null && activeResidue.phi !== undefined && (
                <span><b>Φ</b> {activeResidue.phi.toFixed(1)}°</span>
              )}
              {activeResidue.psi !== null && activeResidue.psi !== undefined && (
                <span><b>Ψ</b> {activeResidue.psi.toFixed(1)}°</span>
              )}
              {activeResidue.asa !== null && activeResidue.asa !== undefined && (
                <span><b>ASA</b> {activeResidue.asa.toFixed(1)}</span>
              )}
            </div>
          </div>
        ) : (
          <div className="inspector-placeholder">
            <span>Hover or click any residue block below to inspect secondary structure & dihedral angles.</span>
          </div>
        )}
      </div>

      {/* Scrollable Linear Sequence Track with Ruler */}
      <div className="secondary-seq-track-container">
        <div className="secondary-ruler">
          {rulerMarkers.map((marker) => {
            const leftPct = (marker.index / Math.max(totalResidues - 1, 1)) * 100;
            return (
              <span key={marker.number} className="ruler-tick" style={{ left: `${leftPct}%` }}>
                <i />
                <b>{marker.number}</b>
              </span>
            );
          })}
        </div>

        <div className="secondary-interactive-track" role="img" aria-label={`${result.method} secondary structure track`}>
          {result.residues.map((residue, index) => {
            const isSelected = selectedResidue === residue.residue_number;
            const isHovered = hoveredResidue?.residue_number === residue.residue_number;
            const isInHoveredBlock =
              hoveredSSEBlock &&
              residue.residue_number >= hoveredSSEBlock.start &&
              residue.residue_number <= hoveredSSEBlock.end;

            return (
              <button
                key={`${residue.chain_id}-${residue.residue_number}-${index}`}
                type="button"
                className={`secondary-seq-block ${isSelected ? 'selected' : ''} ${isHovered ? 'hovered' : ''} ${isInHoveredBlock ? 'block-highlight' : ''}`}
                style={{ backgroundColor: colors[residue.code] ?? colors['-'] }}
                onMouseEnter={() => setHoveredResidue(residue)}
                onMouseLeave={() => setHoveredResidue(null)}
                onClick={() => onSelectResidue?.(residue.residue_number)}
                title={`${residue.chain_id}:${residue.residue_number} ${residue.residue_name} · ${labels[residue.code] ?? residue.code}`}
              />
            );
          })}
        </div>
      </div>

      {/* Secondary Structure Elements (SSE) Block Explorer */}
      <div className="secondary-sse-explorer">
        <div className="sse-explorer-header">
          <span className="analysis-label">STRUCTURAL ELEMENTS (SSE)</span>
          <small>{filteredBlocks.length} elements {activeFilter !== 'ALL' ? `(${activeFilter})` : ''}</small>
        </div>

        <div className="sse-block-list">
          {filteredBlocks.map((block) => (
            <div
              key={block.id}
              className={`sse-pill-chip sse-${block.type}`}
              onMouseEnter={() => setHoveredSSEBlock(block)}
              onMouseLeave={() => setHoveredSSEBlock(null)}
              onClick={() => onSelectResidue?.(block.start)}
              title={`Click to navigate to residue ${block.start}`}
            >
              <span className="sse-chip-badge">{block.label}</span>
              <span className="sse-chip-range">Res {block.start}–{block.end}</span>
              <small className="sse-chip-len">{block.length} aa</small>
            </div>
          ))}
        </div>
      </div>
    </div>
  );
}