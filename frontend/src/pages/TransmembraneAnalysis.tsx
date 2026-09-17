import { useEffect, useMemo, useState } from 'react';
import { SecondaryStructureTrack } from '../components/analysis/SecondaryStructureTrack';
import { TransmembraneTopologyDiagram } from '../components/topology/TransmembraneTopologyDiagram';
import { getHealth, getSecondaryCapabilities, runSecondaryStructure, uploadStructure } from '../services/api';
import type { ProteinUpload } from '../types/protein';
import type { SecondaryStructureResult } from '../types/secondaryStructure';

type Method = 'DSSP' | 'STRIDE' | 'COMPARE' | 'MANUAL';

export function TransmembraneAnalysis() {
  const [protein, setProtein] = useState<ProteinUpload | null>(null);
  const [chainId, setChainId] = useState<string | null>(null);
  const [selectedResidue, setSelectedResidue] = useState<number | null>(null);
  const [method, setMethod] = useState<Method>('DSSP');
  const [status, setStatus] = useState('Ready for a structure file');
  const [health, setHealth] = useState('API status unknown');
  const [secondaryCapabilities, setSecondaryCapabilities] = useState<Record<string, { available: boolean; executable: string }>>({});
  const [secondaryResult, setSecondaryResult] = useState<SecondaryStructureResult | null>(null);
  const [secondaryError, setSecondaryError] = useState<string | null>(null);
  
  const chain = useMemo(() => protein?.models[0]?.chains.find((item) => item.id === chainId) ?? protein?.models[0]?.chains[0], [protein, chainId]);
  const chains = protein?.models[0]?.chains ?? [];

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
      setSecondaryResult(null);
      setSecondaryError(null);
      setStatus(`${file.name} loaded`);
    } catch (error) {
      setStatus(error instanceof Error ? error.message : 'Upload failed');
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
          <h1>Topology Workspace</h1>
          <p>Inspect secondary structure and transmembrane topology annotations.</p>
        </div>
        <div className="header-actions">
          <label className="upload-button"><span>Upload PDB / mmCIF</span><input type="file" accept=".pdb,.ent,.cif,.mmcif" onChange={(event) => event.target.files?.[0] && handleUpload(event.target.files[0])} /></label>
          <button className="quiet-button" onClick={checkHealth}>Check API</button>
        </div>
      </header>

      {chains.length > 0 && (
        <section className="panel" style={{ padding: '16px', marginBottom: '16px' }}>
          <div style={{ display: 'flex', alignItems: 'center', gap: '16px' }}>
            <strong>Select Chain:</strong>
            <select 
              value={chainId || ''} 
              onChange={(e) => { setChainId(e.target.value); setSelectedResidue(null); }}
              style={{ padding: '8px', borderRadius: '4px', border: '1px solid #ccc', backgroundColor: '#fff', color: '#333' }}
            >
              {chains.map(c => (
                <option key={c.id} value={c.id}>Chain {c.id} ({c.residue_count} residues)</option>
              ))}
            </select>
          </div>
        </section>
      )}

      <section className="lower-grid">
        <div className="panel method-selector-panel">
          <div className="panel-heading">
            <div>
              <span className="section-kicker">ANNOTATION METHOD</span>
              <h2>Secondary structure</h2>
              <p className="panel-subtitle">Choose assignment source.</p>
            </div>
          </div>

          <div className="method-pill-group">
            {(['DSSP', 'STRIDE', 'COMPARE', 'MANUAL'] as Method[]).map((item) => {
              const isTool = item === 'DSSP' || item === 'STRIDE';
              const isReady = secondaryCapabilities[item]?.available;

              return (
                <button
                  key={item}
                  className={`method-pill ${method === item ? 'active' : ''} ${isTool && !isReady ? 'missing-tool' : ''}`}
                  onClick={() => { setMethod(item); setSecondaryError(null); }}
                >
                  <span className="method-pill-name">{item === 'COMPARE' ? 'DSSP + STRIDE' : item}</span>
                  <span className={`method-pill-status ${isReady ? 'ready' : ''}`}>
                    {!isTool ? 'MODE' : isReady ? 'READY' : 'NOT INSTALLED'}
                  </span>
                </button>
              );
            })}
          </div>
        </div>

        <div className="panel method-panel" style={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
          <button className="run-button" style={{ padding: '16px', fontSize: '16px' }} onClick={runAnalysis}>
            Run annotation
          </button>

          <div className="method-note">
            {method === 'MANUAL'
              ? 'Manual annotations will map residue ranges to helix, sheet, or turn.'
              : secondaryError
                ? secondaryError
                : secondaryResult
                  ? `${secondaryResult.method} returned ${secondaryResult.residues.length} residue assignments.`
                  : `${method === 'COMPARE' ? 'Comparison' : method} requires its native executable.`}
          </div>

          {secondaryError && (
            <div className="tool-error">
              Install {method === 'DSSP' ? 'mkdssp' : 'stride'} and place it on PATH, then restart the backend.
            </div>
          )}
        </div>
      </section>

      {secondaryResult && (
        <section className="secondary-result-section" style={{ marginTop: '14px' }}>
          <SecondaryStructureTrack
            result={secondaryResult}
            selectedResidue={selectedResidue}
            onSelectResidue={setSelectedResidue}
          />
        </section>
      )}

      <section className="tm-topology-section">
        <TransmembraneTopologyDiagram
          chain={chain}
          secondaryResult={secondaryResult}
          selectedResidue={selectedResidue}
          onSelectResidue={setSelectedResidue}
          uniprotId={protein?.uniprot_id}
          filename={protein?.filename}
        />
      </section>

      <footer className="status-bar"><span><b className="status-dot" />{status}</span><span>{health}</span><span>{selectedResidue ? `Selected residue ${chain?.id}:${selectedResidue}` : 'No residue selected'}</span></footer>
    </main>
  );
}