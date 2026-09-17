import { useState } from 'react';
import { ProteinViewer } from '../components/viewer/ProteinViewer';
import { SequenceView } from '../components/sequence/SequenceView';
import { ExpandedProteinMap } from '../components/topology/ExpandedProteinMap';
import { AnalysisPanel } from '../components/analysis/AnalysisPanel';
import { analyzeChain, getHealth, uploadStructure } from '../services/api';
import type { ChainAnalysis } from '../types/analysis';
import type { ProteinUpload } from '../types/protein';

export function StructureViewer() {
  const [protein, setProtein] = useState<ProteinUpload | null>(null);
  const [chainId, setChainId] = useState<string | null>(null);
  const [selectedResidue, setSelectedResidue] = useState<number | null>(null);
  const [sequenceOpen, setSequenceOpen] = useState(false);
  const [status, setStatus] = useState('Ready for a structure file');
  const [health, setHealth] = useState('API status unknown');
  const [analysis, setAnalysis] = useState<ChainAnalysis | null>(null);
  const [analysisLoading, setAnalysisLoading] = useState(false);

  const chain = protein?.models[0]?.chains.find((item) => item.id === chainId) ?? protein?.models[0]?.chains[0];
  const chains = protein?.models[0]?.chains ?? [];
  const selectedResidueData = chain?.residues.find((item) => item.id === selectedResidue);

  const handleUpload = async (file: File) => {
    setStatus(`Parsing ${file.name}...`);
    try {
      const result = await uploadStructure(file);
      setProtein(result);
      setChainId(result.models[0]?.chains[0]?.id ?? null);
      setSelectedResidue(null);
      setAnalysis(null);
      setSequenceOpen(false);
      setStatus(`${file.name} loaded`);
    } catch (error) {
      setStatus(error instanceof Error ? error.message : 'Upload failed');
    }
  };

  const runChainAnalysis = async () => {
    if (!protein || !chain) {
      setStatus('Upload a structure and select a chain before analysis');
      return;
    }
    setAnalysisLoading(true);
    try {
      setAnalysis(await analyzeChain(protein.filename, chain.id));
      setStatus(`Analysis complete for chain ${chain.id}`);
    } catch (error) {
      setStatus(error instanceof Error ? error.message : 'Analysis failed');
    } finally {
      setAnalysisLoading(false);
    }
  };

  const checkHealth = async () => {
    try {
      const result = await getHealth();
      setHealth(`${result.service} · ${result.version}`);
    } catch {
      setHealth('Backend unavailable');
    }
  };

  return (
    <main className="app-shell">
      <header className="topbar">
        <div>
          <div className="eyebrow">STRUCTURE LAB / MVP 0.1</div>
          <h1>Visualize Workspace</h1>
          <p>Inspect coordinates, sequence, and overall architecture.</p>
        </div>
        <div className="header-actions">
          <label className="upload-button"><span>Upload PDB / mmCIF</span><input type="file" accept=".pdb,.ent,.cif,.mmcif" onChange={(event) => event.target.files?.[0] && handleUpload(event.target.files[0])} /></label>
          <button className="quiet-button" onClick={checkHealth}>Check API</button>
        </div>
      </header>

      <section className="overview-section panel">
        <div className="panel-heading"><div><span className="section-kicker">FULL ASSEMBLY</span><h2>Biological assembly overview</h2><p className="panel-subtitle">The complete colored protein stays together here, like the reference image.</p></div><span className="tag">{chains.length ? `${chains.length} CHAINS` : 'WAITING'}</span></div>
        <ProteinViewer chain={chain} chains={chains} filename={protein?.filename} variant="overview" selectedResidue={selectedResidue} onSelectResidue={setSelectedResidue} />
      </section>

      <section className="spread-section panel">
        <div className="panel-heading"><div><span className="section-kicker">EXPANDED ARCHITECTURE</span><h2>Spread protein map</h2><p className="panel-subtitle">Compounds are separated into readable lanes while the external ribbons preserve their assembly relationships.</p></div><span className="tag">CLICK A LANE</span></div>
        <ExpandedProteinMap chains={chains} selectedChain={chain?.id} onSelectChain={(nextChain) => { setChainId(nextChain); setSelectedResidue(null); setAnalysis(null); }} />
      </section>

      <section className="focused-workspace panel">
        <div className="panel-heading"><div><span className="section-kicker">INTERACTION</span><h2>Focused chain workspace</h2><p className="panel-subtitle">Choose one chain, read its sequence, then inspect only that chain in 3D.</p></div><span className="tag">{chain ? `CHAIN ${chain.id}` : 'WAITING'}</span></div>
        <div className="focused-sequence-layout">
          <div className="focused-sequence-column">
            <button className={sequenceOpen ? 'sequence-disclosure open' : 'sequence-disclosure'} onClick={() => setSequenceOpen((open) => !open)} aria-expanded={sequenceOpen}>
              <span><span className="section-kicker">SELECTED SEQUENCE</span><strong>Chain {chain?.id ?? '-'}</strong><small>{chain?.residue_count ?? 0} residues · {sequenceOpen ? 'Hide sequence' : 'Show sequence'}</small></span>
              <span className="disclosure-chevron">{sequenceOpen ? '−' : '+'}</span>
            </button>
            {sequenceOpen && <div className="sequence-dropdown-content">
              <div className={selectedResidueData ? 'selection-inspector active' : 'selection-inspector'}>
                <span className="selection-label">CURRENT RESIDUE</span>
                <strong>{selectedResidueData ? `${chain?.id}:${selectedResidueData.id} · ${selectedResidueData.name}` : 'Select a residue below'}</strong>
                <small>{selectedResidueData ? 'Highlighted in the isolated 3D chain and topology map.' : 'Click an amino acid to highlight its exact position in 3D.'}</small>
              </div>
              <SequenceView chain={chain} selectedResidue={selectedResidue} onSelectResidue={setSelectedResidue} />
            </div>}
          </div>
          <div className="focused-viewer-column">
            <div className="panel-heading compact-heading"><div><span className="section-kicker">3D ISOLATE</span><h3>Chain {chain?.id ?? '-'}</h3></div><span className="tag">INTERACTIVE</span></div>
            <ProteinViewer chain={chain} chains={chains} filename={protein?.filename} variant="interactive" focusChainId={chain?.id} selectedResidue={selectedResidue} onSelectResidue={setSelectedResidue} />
          </div>
        </div>
      </section>

      <section className="analysis-section panel"><AnalysisPanel analysis={analysis} loading={analysisLoading} onRun={runChainAnalysis} /></section>

      <footer className="status-bar"><span><b className="status-dot" />{status}</span><span>{health}</span><span>{selectedResidue ? `Selected residue ${chain?.id}:${selectedResidue}` : 'No residue selected'}</span></footer>
    </main>
  );
}
