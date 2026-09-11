import type { SecondaryStructureResult } from '../../types/secondaryStructure';

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
  H: '#e76f51', G: '#f4a261', I: '#f6bd60', E: '#4ea8de', B: '#72b7d6',
  T: '#9b8afb', S: '#8aa17b', '-': '#496673',
};

export function SecondaryStructureTrack({ result }: { result: SecondaryStructureResult }) {
  const counts = result.residues.reduce<Record<string, number>>((summary, residue) => {
    summary[residue.code] = (summary[residue.code] ?? 0) + 1;
    return summary;
  }, {});

  return (
    <div className="secondary-track-card">
      <div className="secondary-track-header"><span className="analysis-label">ASSIGNMENT MAP</span><span>{result.residues.length} residues · {result.method}</span></div>
      <div className="secondary-track" role="img" aria-label={`${result.method} secondary structure assignment`}>
        {result.residues.map((residue, index) => <span key={`${residue.chain_id}-${residue.residue_number}-${index}`} className="secondary-block" style={{ backgroundColor: colors[residue.code] ?? colors['-'] }} title={`${residue.chain_id}:${residue.residue_number} ${residue.residue_name} · ${labels[residue.code] ?? residue.code}`} />)}
      </div>
      <div className="secondary-legend">{Object.entries(counts).map(([code, count]) => <span key={code}><i style={{ backgroundColor: colors[code] ?? colors['-'] }} /><b>{code}</b>{labels[code] ?? 'Other'} · {count}</span>)}</div>
    </div>
  );
}