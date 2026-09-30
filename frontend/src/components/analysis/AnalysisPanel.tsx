import type { ChainAnalysis } from '../../types/analysis';
import './AnalysisPanel.css';

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
    <div className="premium-analysis-container">
      <div className="premium-analysis-header">
        <div>
          <h2 className="premium-title">Chain analysis</h2>
          <p className="premium-subtitle">Reproducible sequence, physicochemical, geometry, and contact metrics.</p>
        </div>
        <button className="premium-run-btn" onClick={onRun} disabled={loading}>
          {loading ? 'Running...' : 'Run analysis'}
        </button>
      </div>

      {analysis ? (
        <>
          <div className="premium-metrics-grid">
            <div className="premium-metric-card">
              <span className="metric-label">Length</span>
              <strong className="metric-value">{sequence?.length}</strong>
              <small className="metric-unit">residues</small>
            </div>
            <div className="premium-metric-card">
              <span className="metric-label">Molecular weight</span>
              <strong className="metric-value">{profile?.molecular_weight?.toLocaleString() ?? '-'}</strong>
              <small className="metric-unit">Da</small>
            </div>
            <div className="premium-metric-card">
              <span className="metric-label">pI</span>
              <strong className="metric-value">{profile?.isoelectric_point ?? '-'}</strong>
              <small className="metric-unit">isoelectric point</small>
            </div>
            <div className="premium-metric-card">
              <span className="metric-label">Hydropathy</span>
              <strong className="metric-value">{profile?.mean_hydropathy ?? '-'}</strong>
              <small className="metric-unit">Kyte-Doolittle mean</small>
            </div>
            <div className="premium-metric-card">
              <span className="metric-label">Contacts</span>
              <strong className="metric-value">{geometry?.contacts}</strong>
              <small className="metric-unit">CA pairs under {analysis.result.parameters.contact_cutoff} Å</small>
            </div>
          </div>

          <div className="premium-lower-section">
            <div className="premium-chart-card">
              <span className="chart-header">AMINO ACID COMPOSITION</span>
              <div className="composition-visualizer">
                {Object.entries(sequence?.composition ?? {}).map(([letter, count]) => {
                  const percentage = (count / Math.max(sequence?.length ?? 1, 1)) * 100;
                  return (
                    <div key={letter} className="bar-wrapper">
                      <div className="bar-tooltip">{letter}: {count} ({percentage.toFixed(1)}%)</div>
                      <i className="bar-fill" style={{ height: `${Math.max(8, percentage)}%` }} />
                      <span className="bar-label">{letter}</span>
                    </div>
                  );
                })}
              </div>
            </div>
            
            <div className="reproducibility-card">
              <span className="chart-header">REPRODUCIBILITY</span>
              <strong className="algo-name">{analysis.result.algorithm}</strong>
              <small className="algo-details">
                version {analysis.result.version} · cutoff {analysis.result.parameters.contact_cutoff} Å
              </small>
            </div>
          </div>
        </>
      ) : (
        <div className="analysis-empty" style={{ 
          display: 'grid', 
          minHeight: '116px', 
          placeItems: 'center', 
          border: '1px dashed rgba(51, 65, 85, 0.8)', 
          borderRadius: '12px', 
          color: '#94a3b8',
          background: 'rgba(15, 23, 42, 0.2)' 
        }}>
          Run the toolbox to calculate metrics for the selected chain.
        </div>
      )}
    </div>
  );
}