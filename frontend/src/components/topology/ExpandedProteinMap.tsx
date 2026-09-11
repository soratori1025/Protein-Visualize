import type { Chain } from '../../types/protein';
import type { SecondaryStructureResult } from '../../types/secondaryStructure';

export interface SSElement {
  id: string;
  type: 'helix' | 'strand' | 'coil';
  code: string;
  label: string;
  startRes: number;
  endRes: number;
  length: number;
}

interface Props {
  chains: Chain[];
  selectedChain: string | undefined;
  secondaryResult?: SecondaryStructureResult | null;
  onSelectChain: (chainId: string) => void;
}

export function extractSSElements(
  chainId: string,
  chainResiduesCount: number,
  secondaryResult?: SecondaryStructureResult | null
): SSElement[] {
  const chainResidues = secondaryResult?.residues?.filter((r) => r.chain_id === chainId) ?? [];

  if (!chainResidues.length) {
    const elements: SSElement[] = [];
    let current = 1;
    let helixIdx = 1;
    let strandIdx = 1;

    while (current <= chainResiduesCount) {
      const loopLen = Math.min(8, chainResiduesCount - current + 1);
      if (loopLen > 0) {
        elements.push({
          id: `loop-${current}`,
          type: 'coil',
          code: 'C',
          label: 'Loop',
          startRes: current,
          endRes: current + loopLen - 1,
          length: loopLen,
        });
        current += loopLen;
      }
      if (current > chainResiduesCount) break;

      const isHelix = elements.length % 4 !== 3;
      if (isHelix) {
        const hLen = Math.min(20, chainResiduesCount - current + 1);
        if (hLen >= 4) {
          elements.push({
            id: `helix-${current}`,
            type: 'helix',
            code: 'H',
            label: `α${helixIdx++}`,
            startRes: current,
            endRes: current + hLen - 1,
            length: hLen,
          });
          current += hLen;
        }
      } else {
        const sLen = Math.min(8, chainResiduesCount - current + 1);
        if (sLen >= 3) {
          elements.push({
            id: `strand-${current}`,
            type: 'strand',
            code: 'E',
            label: `β${strandIdx++}`,
            startRes: current,
            endRes: current + sLen - 1,
            length: sLen,
          });
          current += sLen;
        }
      }
    }
    return elements;
  }

  const elements: SSElement[] = [];
  let currentGroup: { type: 'helix' | 'strand' | 'coil'; code: string; start: number; end: number } | null = null;
  let helixCount = 1;
  let strandCount = 1;

  for (const res of chainResidues) {
    const code = res.code.toUpperCase();
    let type: 'helix' | 'strand' | 'coil' = 'coil';
    if (['H', 'G', 'I'].includes(code)) type = 'helix';
    else if (['E', 'B'].includes(code)) type = 'strand';

    if (!currentGroup) {
      currentGroup = { type, code, start: res.residue_number, end: res.residue_number };
    } else if (currentGroup.type === type) {
      currentGroup.end = res.residue_number;
    } else {
      const len = currentGroup.end - currentGroup.start + 1;
      let label = 'Loop';
      if (currentGroup.type === 'helix') label = `α${helixCount++}`;
      else if (currentGroup.type === 'strand') label = `β${strandCount++}`;

      elements.push({
        id: `${currentGroup.type}-${currentGroup.start}`,
        type: currentGroup.type,
        code: currentGroup.code,
        label,
        startRes: currentGroup.start,
        endRes: currentGroup.end,
        length: len,
      });
      currentGroup = { type, code, start: res.residue_number, end: res.residue_number };
    }
  }

  if (currentGroup) {
    const len = currentGroup.end - currentGroup.start + 1;
    let label = 'Loop';
    if (currentGroup.type === 'helix') label = `α${helixCount++}`;
    else if (currentGroup.type === 'strand') label = `β${strandCount++}`;

    elements.push({
      id: `${currentGroup.type}-${currentGroup.start}`,
      type: currentGroup.type,
      code: currentGroup.code,
      label,
      startRes: currentGroup.start,
      endRes: currentGroup.end,
      length: len,
    });
  }

  return elements;
}

const typeColors = {
  helix: '#e76f51',
  strand: '#4ea8de',
  coil: '#496673',
};

export function ExpandedProteinMap({ chains, selectedChain, secondaryResult, onSelectChain }: Props) {
  if (!chains.length) return <div className="empty-state spread-empty">Upload a structure to expand its secondary structure layout.</div>;

  const width = Math.max(920, chains.length * 190);
  const laneWidth = 140;
  const gap = (width - 60 - chains.length * laneWidth) / Math.max(chains.length - 1, 1);
  const laneX = (index: number) => 30 + index * (laneWidth + gap);

  return (
    <div className="spread-map-wrap">
      <svg className="spread-map" viewBox={`0 0 ${width} 330`} role="img" aria-label="Expanded protein secondary structure map">
        {chains.length > 1 && (
          <>
            <path className="outer-connector" d={`M ${laneX(0) + 40} 32 C ${width * 0.25} 2, ${width * 0.68} 2, ${laneX(chains.length - 1) + 100} 32`} />
            <path className="outer-connector faint" d={`M ${laneX(0) + 100} 296 C ${width * 0.42} 326, ${width * 0.68} 326, ${laneX(chains.length - 1) + 40} 296`} />
          </>
        )}

        {chains.map((chain, chainIndex) => {
          const x = laneX(chainIndex);
          const elements = extractSSElements(chain.id, chain.residue_count, secondaryResult);
          const helices = elements.filter((e) => e.type === 'helix');
          const strands = elements.filter((e) => e.type === 'strand');
          const displayedElements = elements.slice(0, 8);

          return (
            <g
              key={chain.id}
              className={selectedChain === chain.id ? 'spread-lane active' : 'spread-lane'}
              onClick={() => onSelectChain(chain.id)}
            >
              <rect className="compound-card" x={x} y="44" width={laneWidth} height="246" rx="12" />
              <text className="compound-label" x={x + laneWidth / 2} y="68" textAnchor="middle">CHAIN {chain.id}</text>
              <text className="sse-summary" x={x + laneWidth / 2} y="84" textAnchor="middle">
                {helices.length}α · {strands.length}β · {chain.residue_count} res
              </text>
              <line className="inner-strip" x1={x + 18} y1="96" x2={x + laneWidth - 18} y2="96" />

              <g transform={`translate(${x + 12}, 106)`}>
                {displayedElements.map((elem, idx) => {
                  const elemY = idx * 20;
                  const color = typeColors[elem.type];
                  return (
                    <g key={elem.id} className="sse-block-group">
                      <title>{`${elem.label} (${elem.startRes}-${elem.endRes}): ${elem.length} residues`}</title>
                      {elem.type === 'helix' ? (
                        <rect className="sse-helix" x="0" y={elemY} width={laneWidth - 24} height="15" rx="7.5" fill={color} />
                      ) : elem.type === 'strand' ? (
                        <rect className="sse-strand" x="0" y={elemY} width={laneWidth - 24} height="15" rx="3" fill={color} />
                      ) : (
                        <rect className="sse-coil" x="12" y={elemY + 4} width={laneWidth - 48} height="7" rx="3.5" fill={color} opacity={0.6} />
                      )}
                      {elem.type !== 'coil' && (
                        <text className="sse-block-text" x={(laneWidth - 24) / 2} y={elemY + 11} textAnchor="middle">
                          {elem.label} ({elem.startRes}-{elem.endRes})
                        </text>
                      )}
                    </g>
                  );
                })}
                {elements.length > 8 && (
                  <text className="sse-more-text" x={(laneWidth - 24) / 2} y={166} textAnchor="middle">
                    + {elements.length - 8} more SSEs
                  </text>
                )}
              </g>

              <text className="residue-count" x={x + laneWidth / 2} y="278" textAnchor="middle">
                Secondary Structure
              </text>
            </g>
          );
        })}
      </svg>
      <div className="spread-caption">
        <span className="spread-key helix-key" style={{ background: '#e76f51' }} /> Alpha Helix (α)
        <span className="spread-key strand-key" style={{ background: '#4ea8de', marginLeft: '12px' }} /> Beta Strand (β)
        <span className="spread-key coil-key" style={{ background: '#496673', marginLeft: '12px' }} /> Loop / Coil
      </div>
    </div>
  );
}