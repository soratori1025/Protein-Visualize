import React, { useMemo, useState } from 'react';
import { ConsensusResidue } from '../../../types/secondaryStructure';
import { ProteinViewer } from '../../viewer/ProteinViewer';
import { useProtein } from '../../../contexts/ProteinContext';
import './ConsensusAnalysisMap.css';
import './ConsensusAnalysisMap.css';

const labelColors: Record<string, string> = {
  TM_E: '#ff9f43', // Extracellular (Orange)
  TM_in: '#00d2d3', // Membrane (Teal)
  TM_C: '#5f27cd', // Cytoplasmic (Purple)
};

const labelY: Record<string, number> = {
  TM_E: 40,
  TM_in: 140,
  TM_C: 240,
};

interface ConsensusAnalysisMapProps {
  consensusMap: ConsensusResidue[];
}

export const ConsensusAnalysisMap: React.FC<ConsensusAnalysisMapProps> = ({ consensusMap }) => {
  const [hoveredRes, setHoveredRes] = useState<ConsensusResidue | null>(null);
  
  const { protein, chainId, selectedResidue, setSelectedResidue } = useProtein();
  const chain = protein?.models[0]?.chains.find((item) => item.id === chainId) ?? protein?.models[0]?.chains[0];
  const chains = protein?.models[0]?.chains ?? [];

  // Layout calculation
  const nodeSpacing = 24;
  const gapSpacing = 64;

  const nodes = useMemo(() => {
    if (!consensusMap || consensusMap.length === 0) return [];
    const result = [];
    let currentX = 40;
    
    for (let i = 0; i < consensusMap.length; i++) {
      const res = consensusMap[i];
      if (i > 0) {
        const prev = consensusMap[i - 1];
        if (res.index === prev.index + 1) {
          currentX += nodeSpacing;
        } else {
          currentX += gapSpacing;
        }
      }
      
      result.push({
        ...res,
        cx: currentX,
        cy: labelY[res.label] || 140,
      });
    }
    return result;
  }, [consensusMap]);

  if (!consensusMap || consensusMap.length === 0) {
    return (
      <div className="consensus-map-empty">
        <p>No consensus data available to map.</p>
      </div>
    );
  }

  const canvasWidth = Math.max(nodes[nodes.length - 1].cx + 60, 800);
  const canvasHeight = 280;

  return (
    <div className="consensus-map-container" style={{ marginTop: '24px' }}>
      <div className="consensus-map-header">
        <h3 className="consensus-map-title">Consensus Merge Analysis</h3>
        <p className="consensus-map-subtitle">
          Complete topology including extra-membrane helices (TM_C, TM_E) and intra-membrane (TM_in).
          Showing exactly <b>{consensusMap.length} residues</b> as distinct nodes.
        </p>
      </div>

      <div style={{ background: '#0f172a', borderRadius: '8px', border: '1px solid #1e293b' }}>
        {/* The SVG scrollable container */}
        <div style={{ overflowX: 'auto' }}>
          <svg
            width={canvasWidth}
            height={canvasHeight}
            viewBox={`0 0 ${canvasWidth} ${canvasHeight}`}
            style={{ minWidth: '100%', display: 'block', background: '#1e293b', borderRadius: '8px 8px 0 0' }}
          >
            {/* Background Bands */}
            <rect x="0" y="0" width={canvasWidth} height="80" fill="#1e293b" />
            <text x="16" y="24" fill="#64748b" fontSize="12" fontWeight="600" letterSpacing="1">EXTRACELLULAR</text>

            <rect x="0" y="80" width={canvasWidth} height="120" fill="#0f172a" />
            <line x1="0" y1="80" x2={canvasWidth} y2="80" stroke="#334155" strokeWidth="1" strokeDasharray="4 4" />
            <line x1="0" y1="200" x2={canvasWidth} y2="200" stroke="#334155" strokeWidth="1" strokeDasharray="4 4" />
            <text x="16" y="104" fill="#00d2d3" fontSize="12" fontWeight="600" letterSpacing="1" opacity="0.7">LIQUID MEMBRANE</text>

            <rect x="0" y="200" width={canvasWidth} height="80" fill="#1e293b" />
            <text x="16" y="270" fill="#64748b" fontSize="12" fontWeight="600" letterSpacing="1">CYTOPLASMIC</text>

            {/* Lines between nodes */}
            {nodes.map((node, i) => {
              if (i === 0) return null;
              const prev = nodes[i - 1];
              const isGap = node.index > prev.index + 1;
              return (
                <line
                  key={`link-${i}`}
                  x1={prev.cx}
                  y1={prev.cy}
                  x2={node.cx}
                  y2={node.cy}
                  stroke="#475569"
                  strokeWidth="2"
                  strokeDasharray={isGap ? "4 4" : "none"}
                />
              );
            })}

            {/* Nodes */}
            {nodes.map((node, i) => {
              const isHovered = hoveredRes === node;
              const color = labelColors[node.label] || '#94a3b8';
              return (
                <g 
                  key={`node-${i}`} 
                  transform={`translate(${node.cx}, ${node.cy})`}
                  onMouseEnter={() => setHoveredRes(node)}
                  onMouseLeave={() => setHoveredRes(null)}
                  style={{ cursor: 'crosshair' }}
                >
                  <rect
                    x="-12"
                    y="-12"
                    width="24"
                    height="24"
                    rx="4"
                    fill={color}
                    stroke={isHovered ? '#ffffff' : color}
                    strokeWidth={isHovered ? 2 : 0}
                  />
                  <text
                    x="0"
                    y="4"
                    fill="#ffffff"
                    fontSize="11"
                    fontWeight="bold"
                    textAnchor="middle"
                    style={{ pointerEvents: 'none' }}
                  >
                    {node.aa}
                  </text>
                  
                  {/* Always show index for the first node, last node, or if hovered */}
                  {(isHovered || i === 0 || i === nodes.length - 1 || node.index % 10 === 0) && (
                    <text
                      x="0"
                      y={node.label === 'TM_E' ? -18 : 26}
                      fill={isHovered ? '#f8fafc' : '#94a3b8'}
                      fontSize="10"
                      fontWeight={isHovered ? 'bold' : 'normal'}
                      textAnchor="middle"
                      style={{ pointerEvents: 'none' }}
                    >
                      {node.index}
                    </text>
                  )}
                </g>
              );
            })}
          </svg>
        </div>

        {/* Legend & Inspector */}
        <div style={{ padding: '16px', display: 'flex', gap: '24px', alignItems: 'center', flexWrap: 'wrap', borderTop: '1px solid #1e293b' }}>
          <div style={{ display: 'flex', gap: '16px' }}>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <div style={{ width: '12px', height: '12px', borderRadius: '2px', background: labelColors.TM_E }} />
              <span style={{ fontSize: '14px', color: '#cbd5e1' }}>Extracellular (TM_E)</span>
            </div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <div style={{ width: '12px', height: '12px', borderRadius: '2px', background: labelColors.TM_in }} />
              <span style={{ fontSize: '14px', color: '#cbd5e1' }}>Transmembrane (TM_in)</span>
            </div>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <div style={{ width: '12px', height: '12px', borderRadius: '2px', background: labelColors.TM_C }} />
              <span style={{ fontSize: '14px', color: '#cbd5e1' }}>Cytoplasmic (TM_C)</span>
            </div>
          </div>

          <div style={{ flex: 1, borderLeft: '1px solid #334155', paddingLeft: '24px', minHeight: '32px', display: 'flex', alignItems: 'center' }}>
            {hoveredRes ? (
              <div style={{ display: 'flex', gap: '20px', alignItems: 'center' }}>
                <span style={{ fontSize: '18px', fontWeight: 'bold', color: labelColors[hoveredRes.label] }}>
                  {hoveredRes.aa}
                </span>
                <span style={{ fontSize: '14px', color: '#f8fafc' }}>
                  Residue: <b>{hoveredRes.residue_number}</b> <span style={{ color: '#94a3b8' }}>(Idx: {hoveredRes.index})</span>
                </span>
                <span style={{ fontSize: '14px', color: '#cbd5e1' }}>
                  Label: <b style={{ color: labelColors[hoveredRes.label] }}>{hoveredRes.label}</b>
                </span>
                <span style={{ fontSize: '14px', color: '#cbd5e1' }}>
                  SS Raw: <b>{hoveredRes.ss_raw}</b>
                </span>
                {hoveredRes.tm_segment && (
                  <span style={{ fontSize: '14px', color: '#94a3b8', border: '1px solid #475569', padding: '2px 6px', borderRadius: '4px' }}>
                    TM Segment {hoveredRes.tm_segment}
                  </span>
                )}
                {hoveredRes.crossing && (
                  <span style={{ fontSize: '14px', color: '#94a3b8', border: '1px solid #475569', padding: '2px 6px', borderRadius: '4px' }}>
                    Crossing {hoveredRes.crossing}
                  </span>
                )}
              </div>
            ) : (
              <span style={{ fontSize: '14px', color: '#64748b' }}>Hover over a residue node to see details</span>
            )}
          </div>
        </div>
      </div>
      
      {chain && (
        <div style={{ marginTop: '24px', background: '#0f172a', borderRadius: '8px', border: '1px solid #1e293b', overflow: 'hidden' }}>
          <div className="consensus-map-header" style={{ padding: '16px', borderBottom: '1px solid #1e293b' }}>
            <h3 className="consensus-map-title">3D Consensus Validation</h3>
            <p className="consensus-map-subtitle">
              Verify the structural separation of Extracellular, Intramembrane, and Cytoplasmic segments in 3D space.
            </p>
          </div>
          <div style={{ height: '400px', position: 'relative' }}>
            <ProteinViewer 
              chain={chain} 
              chains={chains} 
              filename={protein?.filename} 
              variant="interactive" 
              focusChainId={chain.id} 
              selectedResidue={hoveredRes ? hoveredRes.residue_number : selectedResidue} 
              onSelectResidue={setSelectedResidue}
              consensusMap={consensusMap}
              defaultColorScheme="consensus"
            />
          </div>
        </div>
      )}
    </div>
  );
};
