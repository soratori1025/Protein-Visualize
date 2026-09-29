import React, { useMemo, useState } from 'react';
import { ConsensusResidue } from '../../../types/secondaryStructure';
import { ProteinViewer } from '../../viewer/ProteinViewer';
import { useProtein } from '../../../contexts/ProteinContext';
import './ConsensusAnalysisMap.css';

const labelColors: Record<string, string> = {
  TM_E: '#ff9f43',   // Extracellular (Orange)
  TM_in: '#00d2d3',  // Membrane (Teal)
  TM_C: '#5f27cd',   // Cytoplasmic (Purple)
  Turn_E: '#ff6b6b', // Turn Extracellular (Coral)
  Turn_in: '#f368e0', // Turn Membrane (Pink)
  Turn_C: '#c44569', // Turn Cytoplasmic (Magenta)
};

const labelY: Record<string, number> = {
  TM_E: 40,
  TM_in: 140,
  TM_C: 240,
  Turn_E: 40,
  Turn_in: 140,
  Turn_C: 240,
};

// Turn labels use a diamond shape to visually distinguish from TM helix squares
const TURN_LABELS = new Set(['Turn_in', 'Turn_E', 'Turn_C']);

interface ConsensusAnalysisMapProps {
  consensusMap: ConsensusResidue[];
  customHelixColors?: Record<string, string>;
  customRegionColors?: Record<string, string>;
  figureTheme?: 'lab' | 'publication';
}

export const ConsensusAnalysisMap: React.FC<ConsensusAnalysisMapProps> = ({ consensusMap, customHelixColors, customRegionColors, figureTheme }) => {
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
    <div className="consensus-map-container">
      {chain && (
        <div>
            <ProteinViewer 
              chain={chain} 
              chains={chains} 
              filename={protein?.filename} 
              variant="interactive" 
              focusChainId={chain.id} 
              selectedResidue={hoveredRes ? hoveredRes.residue_number : selectedResidue} 
              onSelectResidue={setSelectedResidue}
              consensusMap={consensusMap}
              defaultColorScheme="helices"
              customHelixColors={customHelixColors}
              customRegionColors={customRegionColors}
              figureTheme={figureTheme}
            />
        </div>
      )}
    </div>
  );
};
