import React from 'react';
import type { CalculatedTopologyData } from '../../../types/secondaryStructure';

interface TopologySummaryPanelProps {
  data: CalculatedTopologyData;
}

export function TopologySummaryPanel({ data }: TopologySummaryPanelProps) {
  const numCrossings = data.regions.filter((r) => r.type === 'Transmembrane' || r.membrane_role === 'TM_CROSSING').length;
  
  return (
    <div className="topology-summary-panel" style={{
      background: '#0a1622',
      border: '1px solid #1e293b',
      borderRadius: '8px',
      padding: '16px',
      marginTop: '24px',
      color: '#cbd5e1',
      fontSize: '14px',
      fontFamily: 'system-ui, -apple-system, sans-serif'
    }}>
      <h3 style={{ margin: '0 0 12px 0', fontSize: '16px', color: '#fff' }}>Transmembrane Analysis Summary</h3>
      
      <div style={{ display: 'grid', gridTemplateColumns: '1fr 1fr', gap: '16px', marginBottom: '16px' }}>
        <div>
          <div style={{ color: '#94a3b8', fontSize: '12px', textTransform: 'uppercase', marginBottom: '4px' }}>Algorithm</div>
          <div style={{ fontWeight: 500, color: '#38bdf8' }}>{data.labeler || 'Unknown'}</div>
        </div>
        <div>
          <div style={{ color: '#94a3b8', fontSize: '12px', textTransform: 'uppercase', marginBottom: '4px' }}>Chain Analyzed</div>
          <div style={{ fontWeight: 500 }}>Chain {data.chain_id || 'A'}</div>
        </div>
        <div>
          <div style={{ color: '#94a3b8', fontSize: '12px', textTransform: 'uppercase', marginBottom: '4px' }}>Topology</div>
          <div>{numCrossings} TM Crossings {data.domain_type ? `· ${data.domain_type.replace('_', ' ')}` : ''}</div>
        </div>
        <div>
          <div style={{ color: '#94a3b8', fontSize: '12px', textTransform: 'uppercase', marginBottom: '4px' }}>Membrane Geometry</div>
          <div>
            {data.membrane ? `Thickness: ${(data.membrane.half_thickness * 2).toFixed(1)} Å` : 'N/A'}
            {data.membrane_score ? ` (Score: ${data.membrane_score.toFixed(2)})` : ''}
          </div>
        </div>
      </div>

      {data.warnings && data.warnings.length > 0 && (
        <div style={{ marginTop: '16px', paddingTop: '16px', borderTop: '1px solid #1e293b' }}>
          <div style={{ color: '#fbbf24', fontSize: '12px', textTransform: 'uppercase', marginBottom: '8px', fontWeight: 600 }}>Backend Notes & Validation</div>
          <ul style={{ margin: 0, paddingLeft: '20px', color: '#94a3b8', fontSize: '13px' }}>
            {data.warnings.map((warn, i) => (
              <li key={i} style={{ marginBottom: '4px' }}>{warn}</li>
            ))}
          </ul>
        </div>
      )}

      {/* Chemical Hotspots Layer */}
      {(() => {
        if (!data.residues) return null;
        const chargedAAs = ['D', 'E', 'R', 'K', 'H', 'ASP', 'GLU', 'ARG', 'LYS', 'HIS'];
        const polarAAs = ['S', 'T', 'N', 'Q', 'Y', 'SER', 'THR', 'ASN', 'GLN', 'TYR'];

        const hotspots = data.residues.filter((r) => r.zone === 'CORE' && (chargedAAs.includes(r.aa.toUpperCase()) || polarAAs.includes(r.aa.toUpperCase())));
        if (hotspots.length === 0) return null;

        const charged = hotspots.filter(r => chargedAAs.includes(r.aa.toUpperCase()));
        const polar = hotspots.filter(r => polarAAs.includes(r.aa.toUpperCase()));

        return (
          <div style={{ marginTop: '16px', paddingTop: '16px', borderTop: '1px solid #1e293b' }}>
            <div style={{ color: '#f43f5e', fontSize: '12px', textTransform: 'uppercase', marginBottom: '8px', fontWeight: 600 }}>
              Chemical Hotspots (Membrane Core)
            </div>
            <div style={{ color: '#94a3b8', fontSize: '13px', marginBottom: '8px' }}>
              Unusual charged or polar residues buried deep in the hydrophobic membrane core. These often correlate with ion transport, proton transfer, or binding sites.
            </div>
            
            <div style={{ display: 'flex', gap: '24px', flexWrap: 'wrap' }}>
              {charged.length > 0 && (
                <div>
                  <strong style={{ color: '#fb7185', fontSize: '12px' }}>Charged in Core:</strong>
                  <div style={{ display: 'flex', gap: '6px', flexWrap: 'wrap', marginTop: '4px' }}>
                    {charged.map((r) => (
                      <span key={r.index} style={{ background: 'rgba(244, 63, 94, 0.15)', color: '#fda4af', padding: '2px 6px', borderRadius: '4px', fontSize: '12px' }}>
                        {r.aa}{r.residue_number} ({r.label}, depth: {r.depth?.toFixed(1)}Å)
                      </span>
                    ))}
                  </div>
                </div>
              )}
              {polar.length > 0 && (
                <div>
                  <strong style={{ color: '#38bdf8', fontSize: '12px' }}>Polar in Core:</strong>
                  <div style={{ display: 'flex', gap: '6px', flexWrap: 'wrap', marginTop: '4px' }}>
                    {polar.map((r) => (
                      <span key={r.index} style={{ background: 'rgba(56, 189, 248, 0.15)', color: '#bae6fd', padding: '2px 6px', borderRadius: '4px', fontSize: '12px' }}>
                        {r.aa}{r.residue_number} ({r.label}, depth: {r.depth?.toFixed(1)}Å)
                      </span>
                    ))}
                  </div>
                </div>
              )}
            </div>
          </div>
        );
      })()}
    </div>
  );
}
