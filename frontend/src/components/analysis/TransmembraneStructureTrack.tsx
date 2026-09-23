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

  // Build an array where each element is the type of region for that residue
  const resArray = useMemo(() => {
    const arr: { number: number; type: string; description: string }[] = [];
    
    const residues = chain.residues || [];
    const resNumbers = residues.length > 0 
      ? residues.map(r => r.id) 
      : Array.from({length: chain.sequence?.length || 0}, (_, i) => i + 1);

    resNumbers.forEach(num => {
      arr.push({ number: num, type: 'Unknown', description: 'Unknown domain' });
    });

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

      arr.forEach(item => {
        if (item.number >= d.start && item.number <= d.end) {
          item.type = type;
          item.description = d.description || type;
        }
      });
    });

    // Apply TM / Intramembrane over top
    const tms = (data?.regions || []).filter(r => ['Transmembrane', 'Intramembrane'].includes(r.type));
    tms.forEach(tm => {
      arr.forEach(item => {
        if (item.number >= tm.start && item.number <= tm.end) {
          item.type = tm.type;
          item.description = tm.description || tm.type;
        }
      });
    });

    return arr;
  }, [chain.residues, chain.sequence, data]);

  const categories = useMemo(() => {
    let tm = 0, extra = 0, cyto = 0, intra = 0, unk = 0;
    resArray.forEach(r => {
      if (r.type === 'Transmembrane') tm++;
      else if (r.type === 'Extracellular') extra++;
      else if (r.type === 'Cytoplasmic') cyto++;
      else if (r.type === 'Intramembrane') intra++;
      else unk++;
    });
    const total = resArray.length;
    return {
      tm: { count: tm, pct: total ? ((tm / total) * 100).toFixed(1) : '0' },
      extra: { count: extra, pct: total ? ((extra / total) * 100).toFixed(1) : '0' },
      cyto: { count: cyto, pct: total ? ((cyto / total) * 100).toFixed(1) : '0' },
      intra: { count: intra, pct: total ? ((intra / total) * 100).toFixed(1) : '0' },
    };
  }, [resArray]);

  const rulerMarkers = useMemo(() => {
    const markers: { number: number; index: number }[] = [];
    if (resArray.length === 0) return markers;
    const minRes = resArray[0].number;
    const maxRes = resArray[resArray.length - 1].number;
    
    resArray.forEach((r, index) => {
      if (r.number === minRes || r.number === maxRes || r.number % 50 === 0) {
        markers.push({ number: r.number, index });
      }
    });
    return markers;
  }, [resArray]);

  const activeResIndex = hoveredRes ?? selectedResidue;
  const activeRes = activeResIndex ? resArray.find(r => r.number === activeResIndex) : null;

  return (
    <div className="secondary-track-workspace">
      <div className="secondary-track-header">
        <div className="secondary-title-group">
          <span className="analysis-label">TRANSMEMBRANE MAP</span>
          <strong className="secondary-method-tag">TOPOLOGY ASSIGNMENT</strong>
        </div>
        <div className="secondary-badge">
          <span>{resArray.length} RESIDUES</span>
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
            const leftPct = (marker.index / Math.max(resArray.length - 1, 1)) * 100;
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
