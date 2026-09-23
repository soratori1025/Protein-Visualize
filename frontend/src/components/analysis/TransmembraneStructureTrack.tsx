import { useMemo, useState } from 'react';
import type { UniProtTopologyData } from '../../types/secondaryStructure';
import type { Chain } from '../../types/protein';

const regionColors: Record<string, string> = {
  'Transmembrane': '#e76f51', // Orange/Red
  'Extracellular': '#4ea8de', // Blue
  'Cytoplasmic': '#8aa17b', // Green
  'Intramembrane': '#e9c46a', // Yellow
  'Unknown': '#496673', // Gray/Dark Blue
};

interface Props {
  chain: Chain;
  data: UniProtTopologyData;
  selectedResidue?: number | null;
  onSelectResidue?: (residueNumber: number) => void;
}

export function TransmembraneStructureTrack({ chain, data, selectedResidue, onSelectResidue }: Props) {
  const [hoveredRes, setHoveredRes] = useState<number | null>(null);

  const totalResidues = chain.sequence?.length || 0;

  // Build an array of length N, where each element is the type of region for that residue
  const resArray = useMemo(() => {
    const arr: { number: number; type: string; description: string }[] = [];
    for (let i = 1; i <= totalResidues; i++) {
      arr.push({ number: i, type: 'Unknown', description: 'Unknown domain' });
    }

    // Apply domains (Extracellular / Cytoplasmic)
    const domains = (data?.regions || []).filter(r => r.type === 'Topological domain');
    domains.forEach(d => {
      let type = 'Unknown';
      const desc = (d.description || '').toLowerCase();
      if (['extracellular', 'lumenal', 'luminal', 'periplasmic', 'exoplasmic', 'vesicular'].some(w => desc.includes(w))) {
        type = 'Extracellular';
      } else if (['cytoplasmic', 'intracellular', 'cytosolic', 'matrix', 'stromal', 'nuclear'].some(w => desc.includes(w))) {
        type = 'Cytoplasmic';
      }

      for (let i = Math.max(1, d.start); i <= Math.min(totalResidues, d.end); i++) {
        arr[i - 1].type = type;
        arr[i - 1].description = d.description || type;
      }
    });

    // Apply TM / Intramembrane over top
    const tms = (data?.regions || []).filter(r => ['Transmembrane', 'Intramembrane'].includes(r.type));
    tms.forEach(tm => {
      for (let i = Math.max(1, tm.start); i <= Math.min(totalResidues, tm.end); i++) {
        arr[i - 1].type = tm.type;
        arr[i - 1].description = tm.description || tm.type;
      }
    });

    return arr;
  }, [data, totalResidues]);

  const categories = useMemo(() => {
    let tm = 0, extra = 0, cyto = 0, intra = 0, unk = 0;
    resArray.forEach(r => {
      if (r.type === 'Transmembrane') tm++;
      else if (r.type === 'Extracellular') extra++;
      else if (r.type === 'Cytoplasmic') cyto++;
      else if (r.type === 'Intramembrane') intra++;
      else unk++;
    });
    return {
      tm: { count: tm, pct: totalResidues ? ((tm / totalResidues) * 100).toFixed(1) : '0' },
      extra: { count: extra, pct: totalResidues ? ((extra / totalResidues) * 100).toFixed(1) : '0' },
      cyto: { count: cyto, pct: totalResidues ? ((cyto / totalResidues) * 100).toFixed(1) : '0' },
      intra: { count: intra, pct: totalResidues ? ((intra / totalResidues) * 100).toFixed(1) : '0' },
    };
  }, [resArray, totalResidues]);

  const rulerMarkers = useMemo(() => {
    const markers: { number: number; index: number }[] = [];
    if (totalResidues === 0) return markers;
    for (let i = 1; i <= totalResidues; i++) {
      if (i === 1 || i === totalResidues || i % 50 === 0) {
        markers.push({ number: i, index: i - 1 });
      }
    }
    return markers;
  }, [totalResidues]);

  const activeResIndex = hoveredRes ?? selectedResidue;
  const activeRes = activeResIndex ? resArray[activeResIndex - 1] : null;

  return (
    <div className="secondary-track-workspace">
      <div className="secondary-track-header">
        <div className="secondary-title-group">
          <span className="analysis-label">TRANSMEMBRANE MAP</span>
          <strong className="secondary-method-tag">TOPOLOGY ASSIGNMENT</strong>
        </div>
        <div className="secondary-badge">
          <span>{totalResidues} RESIDUES</span>
        </div>
      </div>

      <div className="secondary-composition-bar-wrap">
        <div className="secondary-composition-bar">
          <div className="comp-segment" style={{ width: `${categories.tm.pct}%`, backgroundColor: regionColors['Transmembrane'] }} title={`Transmembrane: ${categories.tm.count}`} />
          <div className="comp-segment" style={{ width: `${categories.extra.pct}%`, backgroundColor: regionColors['Extracellular'] }} title={`Extracellular: ${categories.extra.count}`} />
          <div className="comp-segment" style={{ width: `${categories.cyto.pct}%`, backgroundColor: regionColors['Cytoplasmic'] }} title={`Cytoplasmic: ${categories.cyto.count}`} />
        </div>
      </div>

      <div className="secondary-category-cards">
        <div className="category-card">
          <div className="cat-icon-dot" style={{ backgroundColor: regionColors['Transmembrane'] }} />
          <div className="cat-info">
            <span className="cat-name">Transmembrane</span>
            <strong className="cat-val">{categories.tm.count} <small>({categories.tm.pct}%)</small></strong>
          </div>
        </div>
        <div className="category-card">
          <div className="cat-icon-dot" style={{ backgroundColor: regionColors['Extracellular'] }} />
          <div className="cat-info">
            <span className="cat-name">Extracellular</span>
            <strong className="cat-val">{categories.extra.count} <small>({categories.extra.pct}%)</small></strong>
          </div>
        </div>
        <div className="category-card">
          <div className="cat-icon-dot" style={{ backgroundColor: regionColors['Cytoplasmic'] }} />
          <div className="cat-info">
            <span className="cat-name">Cytoplasmic</span>
            <strong className="cat-val">{categories.cyto.count} <small>({categories.cyto.pct}%)</small></strong>
          </div>
        </div>
      </div>

      <div className={`secondary-inspector-card ${activeRes ? 'active' : ''}`}>
        {activeRes ? (
          <div className="inspector-content">
            <div className="inspector-main">
              <span className="inspector-res-chip" style={{ backgroundColor: regionColors[activeRes.type] ?? regionColors['Unknown'] }}>
                {activeRes.type.substring(0, 2)}
              </span>
              <div>
                <strong className="inspector-res-title">
                  Chain {chain.id} : Residue {activeRes.number}
                </strong>
                <span className="inspector-ss-label">{activeRes.type} — {activeRes.description}</span>
              </div>
            </div>
          </div>
        ) : (
          <div className="inspector-placeholder">
            <span>Hover or click any residue block below to inspect topology.</span>
          </div>
        )}
      </div>

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

        <div className="secondary-interactive-track">
          {resArray.map((r, index) => {
            const isSelected = selectedResidue === r.number;
            const isHovered = hoveredRes === r.number;
            return (
              <button
                key={`${r.number}-${index}`}
                type="button"
                className={`secondary-seq-block ${isSelected ? 'selected' : ''} ${isHovered ? 'hovered' : ''}`}
                style={{ backgroundColor: regionColors[r.type] ?? regionColors['Unknown'] }}
                onMouseEnter={() => setHoveredRes(r.number)}
                onMouseLeave={() => setHoveredRes(null)}
                onClick={() => onSelectResidue?.(r.number)}
                title={`Residue ${r.number} · ${r.type}`}
              />
            );
          })}
        </div>
      </div>
    </div>
  );
}
