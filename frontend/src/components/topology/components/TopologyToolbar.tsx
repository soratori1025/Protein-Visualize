import React from 'react';
import type { TopologySource, FigureTheme } from '../topologyUtils';
import type { CalculatedTopologyData } from '../../../types/secondaryStructure';

interface Props {
  figureTheme: FigureTheme;
  setFigureTheme: (t: FigureTheme) => void;
  topologySource: TopologySource;
  handleTopologySourceChange: (s: TopologySource) => void;
  filename?: string | null;
  handleExport: (format: 'png' | 'jpeg') => void;
  exporting: boolean;
  colorDrawerOpen: boolean;
  setColorDrawerOpen: (v: boolean | ((prev: boolean) => boolean)) => void;
  uniprotIdInput: string;
  setUniprotIdInput: (v: string) => void;
  fetchUniProtTopology: (id: string) => void;
  loadingUniProt: boolean;
  showUniProtInfo: boolean;
  setShowUniProtInfo: (v: boolean | ((prev: boolean) => boolean)) => void;
  tmAlgorithm: string;
  handleTmAlgorithmChange: (v: string) => void;
  customUniprotId: string;
  setCustomUniprotId: (v: string) => void;
  ssAlgorithm: string;
  setSsAlgorithm: (v: string) => void;
  fetchCalculatedTopology: () => void;
  loadingCalculated: boolean;
  showAdvancedParams: boolean;
  setShowAdvancedParams: (v: boolean) => void;
  tmThickness: string;
  setTmThickness: (v: string) => void;
  tmMinMembraneScore: string;
  setTmMinMembraneScore: (v: string) => void;
  tmMinCrossSpan: string;
  setTmMinCrossSpan: (v: string) => void;
  tmFullCrossFrac: string;
  setTmFullCrossFrac: (v: string) => void;
  tmTreatTurnAsHelix: boolean;
  setTmTreatTurnAsHelix: (v: boolean) => void;
  calculatedData: CalculatedTopologyData | null;
  slabParamsActive: boolean;
  crossParamsActive: boolean;
}

export function TopologyToolbar({
  figureTheme, setFigureTheme, topologySource, handleTopologySourceChange, filename,
  handleExport, exporting, colorDrawerOpen, setColorDrawerOpen,
  uniprotIdInput, setUniprotIdInput, fetchUniProtTopology, loadingUniProt, showUniProtInfo, setShowUniProtInfo,
  tmAlgorithm, handleTmAlgorithmChange, customUniprotId, setCustomUniprotId, ssAlgorithm, setSsAlgorithm,
  fetchCalculatedTopology, loadingCalculated, showAdvancedParams, setShowAdvancedParams,
  tmThickness, setTmThickness, tmMinMembraneScore, setTmMinMembraneScore,
  tmMinCrossSpan, setTmMinCrossSpan, tmFullCrossFrac, setTmFullCrossFrac,
  tmTreatTurnAsHelix, setTmTreatTurnAsHelix, calculatedData, slabParamsActive, crossParamsActive
}: Props) {
  const usedParam = (key: string): string => {
    const value = calculatedData?.parameters_used?.[key];
    return value === undefined || value === null ? 'default' : String(value);
  };

  return (
    <div className="tm-toolbar-actions">
      <div className="tm-preset-chip" style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
        <button className={`tm-tab-btn ${figureTheme === 'publication' ? 'active' : ''}`} onClick={() => setFigureTheme('publication')}>Light Mode</button>
        <button className={`tm-tab-btn ${figureTheme === 'lab' ? 'active' : ''}`} onClick={() => setFigureTheme('lab')}>Dark Mode</button>
      </div>

      <div className="tm-preset-chip" style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
        <button className={`tm-tab-btn ${topologySource === 'calculated' ? 'active' : ''}`} onClick={() => handleTopologySourceChange('calculated')} disabled={!filename} title={!filename ? 'Upload a structure file to calculate topology' : undefined}>Calculated (beta)</button>
        <button className={`tm-tab-btn ${topologySource === 'uniprot' ? 'active' : ''}`} onClick={() => handleTopologySourceChange('uniprot')}>UniProt</button>
      </div>

      <button className="tm-color-toggle-btn" onClick={() => handleExport('png')} disabled={exporting}>
        <span>{exporting ? 'Exporting…' : 'Export PNG'}</span>
      </button>

      <button className={`tm-color-toggle-btn ${colorDrawerOpen ? 'active' : ''}`} onClick={() => setColorDrawerOpen(open => !open)}>
        <span>{colorDrawerOpen ? 'Close customize bar' : 'Customize styles'}</span>
      </button>

      {topologySource === 'uniprot' ? (
        <div style={{ display: 'flex', alignItems: 'center', gap: '6px', position: 'relative' }}>
          <input type="text" value={uniprotIdInput} onChange={(e) => setUniprotIdInput(e.target.value)} placeholder="e.g. P31645" className="tm-input-field" style={{ width: '96px', fontWeight: 700, textTransform: 'uppercase' }} onKeyDown={(e) => { if (e.key === 'Enter') fetchUniProtTopology(uniprotIdInput); }} />
          <button onClick={() => fetchUniProtTopology(uniprotIdInput)} disabled={loadingUniProt} className="tm-add-btn" style={{ background: '#3b82f6' }}>{loadingUniProt ? 'Loading...' : 'Load UniProt'}</button>
          <button onClick={() => setShowUniProtInfo(open => !open)} title="How the UniProt lookup works" className="tm-info-btn" style={{ background: showUniProtInfo ? '#38bdf8' : '#1e293b', color: showUniProtInfo ? '#0f172a' : '#94a3b8' }}>?</button>
          {showUniProtInfo && (
            <div className="tm-info-popover">
              <strong>How the UniProt lookup works</strong>
              Enter any UniProt accession (for example <code>P31645</code>) and press Load. The app reads that entry's curated Transmembrane, Topological domain and Intramembrane features and draws the map from them — no structure file needed. For a structure-derived estimate, switch to Calculated (beta).
              <button onClick={() => setShowUniProtInfo(false)}>Close</button>
            </div>
          )}
        </div>
      ) : (
        <>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '8px', border: '1px solid #333', padding: '8px', borderRadius: '4px' }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px', flexWrap: 'wrap' }}>
              <select className="tm-input-field" value={tmAlgorithm} onChange={(e) => handleTmAlgorithmChange(e.target.value)} style={{ padding: '4px 8px' }}>
                <option value="uniprot_api">UniProt API</option>
              </select>
              {tmAlgorithm === 'uniprot_api' && (
                <input type="text" className="tm-input-field" placeholder="UniProt ID (auto from file)" value={customUniprotId} onChange={(e) => setCustomUniprotId(e.target.value)} style={{ width: '150px', padding: '4px 8px' }} />
              )}
              <select className="tm-input-field" value={ssAlgorithm} onChange={(e) => setSsAlgorithm(e.target.value)} style={{ padding: '4px 8px' }}>
                <option value="dssp">DSSP</option>
                <option value="stride">STRIDE</option>
                <option value="none">None (Only TM boundaries)</option>
              </select>
            </div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <button onClick={fetchCalculatedTopology} disabled={loadingCalculated || !filename} className="tm-add-btn">{loadingCalculated ? 'Computing…' : 'Recalculate'}</button>
              <button onClick={() => setShowAdvancedParams(!showAdvancedParams)} className="tm-add-btn" style={{ fontSize: '0.8em', opacity: 0.8 }} title="Tune biological thresholds for TM detection">{showAdvancedParams ? '▲ Parameters' : '▼ Parameters'}</button>
              {tmAlgorithm === 'uniprot_api' && !customUniprotId && (
                <span style={{ color: '#94a3b8', fontSize: '12px' }}>Empty = read from the file (DBREF / _struct_ref)</span>
              )}
            </div>
          </div>
          {showAdvancedParams && (
            <div style={{ margin: '8px 0', padding: '10px 14px', background: 'rgba(100,100,140,0.08)', borderRadius: '8px', display: 'grid', gridTemplateColumns: '1fr 1fr 1fr', gap: '8px 16px', fontSize: '0.82em' }}>
              <label title="Hydrophobic core thickness (Å)." style={{ opacity: slabParamsActive ? 1 : 0.45 }}>
                Membrane thickness (Å)
                <input type="number" step="0.5" min="20" max="40" value={tmThickness} placeholder={usedParam('membrane_thickness')} disabled={!slabParamsActive} onChange={e => setTmThickness(e.target.value)} className="tm-input-field" style={{ width: '70px', marginLeft: 4 }} />
              </label>
              <label title="Min mean hydrophobicity inside the slab." style={{ opacity: slabParamsActive ? 1 : 0.45 }}>
                Min membrane score
                <input type="number" step="0.1" min="-3" max="3" value={tmMinMembraneScore} placeholder={usedParam('min_membrane_score')} disabled={!slabParamsActive} onChange={e => setTmMinMembraneScore(e.target.value)} className="tm-input-field" style={{ width: '60px', marginLeft: 4 }} />
              </label>
              <label title="Hairpin guard: fraction of the thickness two TM_in runs must span together to count as one traversal of the bilayer." style={{ opacity: crossParamsActive ? 1 : 0.45 }}>
                Min cross span
                <input type="number" step="0.05" min="0.05" max="1" value={tmMinCrossSpan} disabled={!crossParamsActive} placeholder={usedParam('min_cross_span_frac')} onChange={e => setTmMinCrossSpan(e.target.value)} className="tm-input-field" style={{ width: '60px', marginLeft: 4 }} />
              </label>
              <label title="Hairpin guard: a TM_in run spanning this fraction of the thickness crosses the bilayer on its own." style={{ opacity: crossParamsActive ? 1 : 0.45 }}>
                Full cross frac
                <input type="number" step="0.05" min="0.1" max="1.5" value={tmFullCrossFrac} disabled={!crossParamsActive} placeholder={usedParam('full_cross_frac')} onChange={e => setTmFullCrossFrac(e.target.value)} className="tm-input-field" style={{ width: '60px', marginLeft: 4 }} />
              </label>
              <label title="Turn (T) / Bend (S) with a helix on both sides becomes helix." style={{ opacity: 1, display: 'flex', alignItems: 'center', gap: '6px' }}>
                <input type="checkbox" checked={tmTreatTurnAsHelix} onChange={e => setTmTreatTurnAsHelix(e.target.checked)} /> Treat Turn/Bend as Helix
              </label>
            </div>
          )}
        </>
      )}
    </div>
  );
}
