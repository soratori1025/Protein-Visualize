import type { Chain } from '../../types/protein';
import type { SecondaryStructureResult } from '../../types/secondaryStructure';
import { residueLabel } from '../../types/secondaryStructure';

export interface SSElement {
  id: string;
  type: 'helix' | 'strand' | 'coil';
  code: string;
  label: string;
  startRes: number;
  endRes: number;
  /** "100A"-style labels (insertion codes kept). */
  startLabel: string;
  endLabel: string;
  /** Number of assigned residues in the element (not end - start + 1). */
  length: number;
}

interface Props {
  chains: Chain[];
  selectedChain: string | undefined;
  secondaryResult?: SecondaryStructureResult | null;
  onSelectChain: (chainId: string) => void;
}

/**
 * Group consecutive residues of one chain into helix / strand / coil elements.
 *
 * A new element starts whenever the class changes OR the residue numbering jumps
 * (unresolved residues): residues on both sides of a gap are not covalently
 * connected, so helix 10–30 and helix 41–60 around a missing loop are two helices,
 * not one 10–60 helix. With no assignment the result is empty — nothing is invented.
 */
export function extractSSElements(
  chainId: string,
  _chainResiduesCount: number,
  secondaryResult?: SecondaryStructureResult | null
): SSElement[] {
  const chainResidues = secondaryResult?.residues?.filter((r) => r.chain_id === chainId) ?? [];
  if (!chainResidues.length) return [];

  type Group = {
    type: 'helix' | 'strand' | 'coil';
    code: string;
    start: number;
    startIcode: string;
    end: number;
    endIcode: string;
    count: number;
  };
  const elements: SSElement[] = [];
  let current: Group | null = null;
  let helixCount = 1;
  let strandCount = 1;

  const flush = (g: Group) => {
    let label = 'Loop';
    if (g.type === 'helix') label = `α${helixCount++}`;
    else if (g.type === 'strand') label = `β${strandCount++}`;
    elements.push({
      id: `${g.type}-${residueLabel(g.start, g.startIcode)}`,
      type: g.type,
      code: g.code,
      label,
      startRes: g.start,
      endRes: g.end,
      startLabel: residueLabel(g.start, g.startIcode),
      endLabel: residueLabel(g.end, g.endIcode),
      length: g.count,
    });
  };

  for (const res of chainResidues) {
    const code = (res.code ?? '').toUpperCase();
    let type: Group['type'] = 'coil';
    if (['H', 'G', 'I'].includes(code)) type = 'helix';
    else if (['E', 'B'].includes(code)) type = 'strand';
    const icode = res.insertion_code ?? '';

    // same number (insertion code) or +1 = contiguous; anything else is a gap
    const contiguous =
      current !== null &&
      (res.residue_number === current.end || res.residue_number === current.end + 1);

    if (current && current.type === type && contiguous) {
      current.end = res.residue_number;
      current.endIcode = icode;
      current.count += 1;
    } else {
      if (current) flush(current);
      current = { type, code, start: res.residue_number, startIcode: icode, end: res.residue_number, endIcode: icode, count: 1 };
    }
  }
  if (current) flush(current);
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
                {elements.length === 0 && (
                  <text className="sse-more-text" x={(laneWidth - 24) / 2} y={20} textAnchor="middle">
                    Run DSSP / STRIDE
                  </text>
                )}
                {displayedElements.map((elem, idx) => {
                  const elemY = idx * 20;
                  const color = typeColors[elem.type];
                  return (
                    <g key={elem.id} className="sse-block-group">
                      <title>{`${elem.label} (${elem.startLabel}-${elem.endLabel}): ${elem.length} residues`}</title>
                      {elem.type === 'helix' ? (
                        <rect className="sse-helix" x="0" y={elemY} width={laneWidth - 24} height="15" rx="7.5" fill={color} />
                      ) : elem.type === 'strand' ? (
                        <rect className="sse-strand" x="0" y={elemY} width={laneWidth - 24} height="15" rx="3" fill={color} />
                      ) : (
                        <rect className="sse-coil" x="12" y={elemY + 4} width={laneWidth - 48} height="7" rx="3.5" fill={color} opacity={0.6} />
                      )}
                      {elem.type !== 'coil' && (
                        <text className="sse-block-text" x={(laneWidth - 24) / 2} y={elemY + 11} textAnchor="middle">
                          {elem.label} ({elem.startLabel}-{elem.endLabel})
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