import { useEffect, useMemo, useState } from 'react';
import { ProteinViewer } from './components/viewer/ProteinViewer';
import { SequenceView } from './components/sequence/SequenceView';
import { TopologyView } from './components/topology/TopologyView';
import { ExpandedProteinMap } from './components/topology/ExpandedProteinMap';
import { AnalysisPanel } from './components/analysis/AnalysisPanel';
import { SecondaryStructureTrack } from './components/analysis/SecondaryStructureTrack';
import { analyzeChain, getHealth, getSecondaryCapabilities, runSecondaryStructure, uploadStructure } from './services/api';
import type { ChainAnalysis } from './types/analysis';
import type { ProteinUpload } from './types/protein';
import type { SecondaryStructureResult } from './types/secondaryStructure';

type Method = 'DSSP' | 'STRIDE' | 'COMPARE' | 'MANUAL';

function App() {
  const [protein, setProtein] = useState<ProteinUpload | null>(null);
  const [chainId, setChainId] = useState<string | null>(null);
  const [selectedResidue, setSelectedResidue] = useState<number | null>(null);
  const [method, setMethod] = useState<Method>('DSSP');
  const [sequenceOpen, setSequenceOpen] = useState(false);
  const [status, setStatus] = useState('Ready for a structure file');
  const [health, setHealth] = useState('API status unknown');
  const [analysis, setAnalysis] = useState<ChainAnalysis | null>(null);
  const [analysisLoading, setAnalysisLoading] = useState(false);
  const [secondaryCapabilities, setSecondaryCapabilities] = useState<Record<string, { available: boolean; executable: string }>>({});
  const [secondaryResult, setSecondaryResult] = useState<SecondaryStructureResult | null>(null);
  const [secondaryError, setSecondaryError] = useState<string | null>(null);
  const chain = useMemo(() => protein?.models[0]?.chains.find((item) => item.id === chainId) ?? protein?.models[0]?.chains[0], [protein, chainId]);
  const chains = protein?.models[0]?.chains ?? [];
  const selectedResidueData = chain?.residues.find((item) => item.id === selectedResidue);

  useEffect(() => {
    void getSecondaryCapabilities().then(setSecondaryCapabilities).catch(() => setSecondaryCapabilities({}));
  }, []);

  const handleUpload = async (file: File) => {
    setStatus(`Parsing ${file.name}...`);
    try {
      const result = await uploadStructure(file);
      setProtein(result);
      setChainId(result.models[0]?.chains[0]?.id ?? null);
      setSelectedResidue(null);
      setAnalysis(null);
      setSecondaryResult(null);
      setSecondaryError(null);
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
      setSecondaryCapabilities(await getSecondaryCapabilities());
    } catch {
      setHealth('Backend unavailable');
    }
  };

  const runAnalysis = async () => {
    if (!protein || (method !== 'DSSP' && method !== 'STRIDE')) {
      setStatus(method === 'COMPARE' ? 'Compare mode will run after both adapters return results' : 'Upload a structure before running analysis');
      return;
    }
    try {
      const result = await runSecondaryStructure(protein.filename, method);
      setSecondaryResult(result);
      setSecondaryError(null);
      setStatus(`${result.method} assigned ${result.residues?.length ?? 0} residues`);
    } catch (error) {
      setSecondaryError(error instanceof Error ? error.message : 'Secondary-structure analysis failed');
      setStatus(error instanceof Error ? error.message : 'Secondary-structure analysis failed');
    }
  };

  return (
    <main className="app-shell">
      <header className="topbar">
        <div>
          <div className="eyebrow">STRUCTURE LAB / MVP 0.1</div>
          <h1>ProteinLab</h1>
          <p>Inspect coordinates, sequence, and structural annotations in one workspace.</p>
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
        <ExpandedProteinMap chains={chains} selectedChain={chain?.id} secondaryResult={secondaryResult} onSelectChain={(nextChain) => { setChainId(nextChain); setSelectedResidue(null); setAnalysis(null); }} />
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

      <section className="lower-grid">
        <div className="panel topology-panel"><div className="panel-heading"><div><span className="section-kicker">RESIDUE MAP</span><h2>Topology</h2><p className="panel-subtitle">A compact residue rail for locating the selected position across the chain.</p></div><span className="legend"><i /> selected residue</span></div><TopologyView chain={chain} selectedResidue={selectedResidue} onSelectResidue={setSelectedResidue} /></div>
        <div className="panel method-panel"><div className="panel-heading"><div><span className="section-kicker">ANNOTATION METHOD</span><h2>Secondary structure</h2><p className="panel-subtitle">Compare helix, strand, turn, and loop assignments from standard coordinate-based methods.</p></div></div><div className="method-options">{(['DSSP', 'STRIDE', 'COMPARE', 'MANUAL'] as Method[]).map((item) => <button key={item} className={method === item ? 'method-option active' : 'method-option'} onClick={() => { setMethod(item); setSecondaryError(null); }}><span className="radio" />{item === 'COMPARE' ? 'DSSP + STRIDE' : item}<small className={secondaryCapabilities[item]?.available ? 'tool-status ready' : 'tool-status'}>{item === 'COMPARE' || item === 'MANUAL' ? 'mode' : secondaryCapabilities[item]?.available ? 'ready' : 'not installed'}</small></button>)}</div><div className="method-note">{method === 'MANUAL' ? 'Manual annotations will map residue ranges to helix, sheet, or turn.' : secondaryError ? secondaryError : secondaryResult ? `${secondaryResult.method} returned ${secondaryResult.residues.length} residue assignments.` : `${method === 'COMPARE' ? 'Comparison' : method} requires its native executable.`}</div>{secondaryError && <div className="tool-error">Install {method === 'DSSP' ? 'mkdssp' : 'stride'} and place it on PATH, then restart the backend.</div>}{secondaryResult && <SecondaryStructureTrack result={secondaryResult} />}<button className="run-button" onClick={runAnalysis}>Run annotation</button></div>
      </section>

      <section className="analysis-section panel"><AnalysisPanel analysis={analysis} loading={analysisLoading} onRun={runChainAnalysis} /></section>

      <footer className="status-bar"><span><b className="status-dot" />{status}</span><span>{health}</span><span>{selectedResidue ? `Selected residue ${chain?.id}:${selectedResidue}` : 'No residue selected'}</span></footer>
    </main>
  );
}

export default App;
