import type { Chain } from '../../types/protein';

export function TopologyView({ chain, selectedResidue, onSelectResidue }: { chain: Chain | undefined; selectedResidue: number | null; onSelectResidue: (id: number) => void }) {
  if (!chain) return <div className="empty-state">Topology will appear after upload.</div>;
  const width = Math.max(760, chain.residues.length * 8);
  const step = (width - 52) / Math.max(chain.residues.length - 1, 1);
  return (
    <div className="topology-scroll polished-topology">
      <svg viewBox={`0 0 ${width} 150`} role="img" aria-label="Protein topology">
        <rect className="topology-track" x="12" y="48" width={width - 24} height="42" rx="21" />
        <line className="topology-axis" x1="26" y1="69" x2={width - 26} y2="69" />
        {chain.residues.map((residue, index) => {
          const x = 26 + index * step;
          const isSelected = selectedResidue === residue.id;
          const showLabel = index % 10 === 0 || isSelected;
          return <g key={residue.id} className={isSelected ? 'topology-node selected' : 'topology-node'} onClick={() => onSelectResidue(residue.id)}>
            <circle className="topology-residue" cx={x} cy="69" r={isSelected ? 8 : 4} />
            {showLabel && <text className="topology-label" x={x} y="31" textAnchor="middle">{residue.id}</text>}
            {showLabel && <text className="topology-aa" x={x} y="120" textAnchor="middle">{chain.sequence[index] ?? 'X'}</text>}
          </g>;
        })}
      </svg>
      <div className="topology-footer"><span>CHAIN {chain.id}</span><span>{chain.residue_count} residues</span><span>Click a node to link it to the 3D isolate</span></div>
    </div>
  );
}