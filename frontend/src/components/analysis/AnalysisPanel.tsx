import type { ChainAnalysis } from '../../types/analysis';

interface Props {
  analysis: ChainAnalysis | null;
  loading: boolean;
  onRun: () => void;
}

export function AnalysisPanel({ analysis, loading, onRun }: Props) {
  const profile = analysis?.result.physicochemical;
  const sequence = analysis?.result.sequence;
  const geometry = analysis?.result.geometry;
  return (
    <div className="analysis-toolbox">
      <div className="analysis-actions"><div><span className="section-kicker">RESEARCH TOOLBOX</span><h2>Chain analysis</h2><p className="panel-subtitle">Reproducible sequence, physicochemical, geometry, and contact metrics.</p></div><button className="run-button analysis-run" onClick={onRun} disabled={loading}>{loading ? 'Running...' : 'Run analysis'}</button></div>
      {analysis ? <>
        <div className="analysis-metrics">
          <div><span>Length</span><strong>{sequence?.length}</strong><small>residues</small></div>
          <div><span>Molecular weight</span><strong>{profile?.molecular_weight?.toLocaleString() ?? '-'} </strong><small>Da</small></div>
          <div><span>pI</span><strong>{profile?.isoelectric_point ?? '-'}</strong><small>isoelectric point</small></div>
          <div><span>Hydropathy</span><strong>{profile?.mean_hydropathy ?? '-'}</strong><small>Kyte-Doolittle mean</small></div>
          <div><span>Contacts</span><strong>{geometry?.contacts}</strong><small>CA pairs under {analysis.result.parameters.contact_cutoff} A</small></div>
        </div>
        <div className="analysis-lower"><div><span className="analysis-label">AMINO ACID COMPOSITION</span><div className="composition-bars">{Object.entries(sequence?.composition ?? {}).map(([letter, count]) => <div key={letter} title={`${letter}: ${count}`}><i style={{ height: `${Math.max(8, (count / Math.max(sequence?.length ?? 1, 1)) * 100)}%` }} /><span>{letter}</span></div>)}</div></div><div className="reproducibility-card"><span className="analysis-label">REPRODUCIBILITY</span><strong>{analysis.result.algorithm}</strong><small>version {analysis.result.version} · cutoff {analysis.result.parameters.contact_cutoff} A</small></div></div>
      </> : <div className="analysis-empty">Run the toolbox to calculate metrics for the selected chain.</div>}
    </div>
  );
}