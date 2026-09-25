import { useEffect, useMemo, useState } from 'react';
import { SecondaryStructureTrack } from '../components/analysis/SecondaryStructureTrack';
import { TransmembraneStructureTrack } from '../components/analysis/TransmembraneStructureTrack';
import { TransmembraneTopologyDiagram } from '../components/topology/TransmembraneTopologyDiagram';
import { getSecondaryCapabilities, runSecondaryStructure } from '../services/api';
import type { SecondaryStructureResult, UniProtTopologyData } from '../types/secondaryStructure';
import { useProtein } from '../contexts/ProteinContext';
import { Header } from '../components/layout/Header';

type Method = 'DSSP' | 'STRIDE';

export function TransmembraneAnalysis() {
  const { protein, chainId, setChainId, selectedResidue, setSelectedResidue, status, setStatus, health, secondaryResult, setSecondaryResult, activeTopologyData, setActiveTopologyData } = useProtein();
  const [method, setMethod] = useState<Method>('DSSP');
  const [secondaryCapabilities, setSecondaryCapabilities] = useState<Record<string, { available: boolean; executable: string }>>({});
  const [secondaryError, setSecondaryError] = useState<string | null>(null);
  const [tmAlgorithm, setTmAlgorithm] = useState<string>('3d_slab_geom');
  const [topologySource, setTopologySource] = useState<'uniprot' | 'calculated'>('uniprot');
  const [triggerTmRecalc, setTriggerTmRecalc] = useState(0);
  const [distinguishTurns, setDistinguishTurns] = useState(false);
  const [isSSRunning, setIsSSRunning] = useState(false);
  const [isTMRunning, setIsTMRunning] = useState(false);
  
  const chain = useMemo(() => protein?.models[0]?.chains.find((item) => item.id === chainId) ?? protein?.models[0]?.chains[0], [protein, chainId]);
  const chains = protein?.models[0]?.chains ?? [];

  useEffect(() => {
    void getSecondaryCapabilities().then(setSecondaryCapabilities).catch(() => setSecondaryCapabilities({}));
  }, []);

  const runAnalysis = async () => {
    if (!protein || (method !== 'DSSP' && method !== 'STRIDE')) {
      setStatus('Upload a structure before running analysis');
      return;
    }
    setIsSSRunning(true);
    try {
      const result = await runSecondaryStructure(protein.filename, method);
      setSecondaryResult(result);
      setSecondaryError(null);
      setStatus(`${result.method} assigned ${result.residues?.length ?? 0} residues`);
    } catch (error) {
      setSecondaryError(error instanceof Error ? error.message : 'Secondary-structure analysis failed');
      setStatus(error instanceof Error ? error.message : 'Secondary-structure analysis failed');
    } finally {
      setIsSSRunning(false);
    }
  };

  const runTmAnalysis = () => {
    if (!protein) {
      setStatus('Upload a structure before running analysis');
      return;
    }
    setTopologySource('calculated');
    setTriggerTmRecalc(t => t + 1);
  };

  return (
    <main className="app-shell">
      <Header 
        title="Transmembrane Analysis" 
        subtitle="Inspect secondary structure and transmembrane topology annotations." 
      />

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
            {(['DSSP', 'STRIDE'] as Method[]).map((item) => {
              const isTool = item === 'DSSP' || item === 'STRIDE';
              const isReady = secondaryCapabilities[item]?.available;

              return (
                <button
                  key={item}
                  className={`method-pill ${method === item ? 'active' : ''} ${isTool && !isReady ? 'missing-tool' : ''}`}
                  onClick={() => { setMethod(item); setSecondaryError(null); }}
                >
                  <span className="method-pill-name">{item}</span>
                  <span className={`method-pill-status ${isReady ? 'ready' : ''}`}>
                    {!isTool ? 'MODE' : isReady ? 'READY' : 'NOT INSTALLED'}
                  </span>
                </button>
              );
            })}
          </div>
        </div>

        <div className="panel method-panel" style={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>
          <button className="run-button" style={{ padding: '16px', fontSize: '16px' }} onClick={runAnalysis} disabled={isSSRunning}>
            <span>{isSSRunning ? 'Running...' : 'Run SS Analysis'}</span>
          </button>
          {secondaryError && (
            <div className="tool-error">
              {secondaryError}
              {!secondaryCapabilities[method]?.available && (
                <> — install {method === 'DSSP' ? 'mkdssp' : 'stride'} and place it on PATH, then restart the backend.</>
              )}
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

      <section className="lower-grid" style={{ marginTop: '14px' }}>
        <div className="panel method-selector-panel">
          <div className="panel-heading">
            <div>
              <span className="section-kicker">ANNOTATION METHOD</span>
              <h2>Transmembrane topology</h2>
              <p className="panel-subtitle">Choose prediction algorithm.</p>
            </div>
          </div>

          <div className="method-pill-group">
            {[
              { id: '3d_slab_geom', label: '3D Slab Geometry' },
              { id: 'kyte_doolittle_seq', label: 'Kyte-Doolittle' },
              { id: 'uniprot_api', label: 'UniProt API' }
            ].map((item) => (
              <button
                key={item.id}
                className={`method-pill ${tmAlgorithm === item.id ? 'active' : ''}`}
                onClick={() => setTmAlgorithm(item.id)}
              >
                <span className="method-pill-name">{item.label}</span>
                <span className={`method-pill-status ${tmAlgorithm === item.id ? 'ready' : ''}`}>
                  MODE
                </span>
              </button>
            ))}
          </div>
        </div>

        <div className="panel method-panel" style={{ display: 'flex', flexDirection: 'column', justifyContent: 'center' }}>

          <button className="run-button" style={{ padding: '16px', fontSize: '16px' }} onClick={runTmAnalysis} disabled={isTMRunning}>
            <span>{isTMRunning ? 'Running...' : 'Run TM Analysis'}</span>
          </button>
        </div>
      </section>

      {activeTopologyData && chain && (
        <section className="secondary-result-section" style={{ marginTop: '14px' }}>
          <TransmembraneStructureTrack
            chain={chain}
            data={activeTopologyData}
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
          tmAlgorithm={tmAlgorithm}
          onTmAlgorithmChange={setTmAlgorithm}
          topologySource={topologySource}
          onTopologySourceChange={setTopologySource}
          onTopologyDataChange={setActiveTopologyData}
          triggerTmRecalc={triggerTmRecalc}
          distinguishTurns={distinguishTurns}
          onLoadingChange={setIsTMRunning}
        />
      </section>

      <footer className="status-bar"><span><b className="status-dot" />{status}</span><span>{health}</span><span>{selectedResidue ? `Selected residue ${chain?.id}:${selectedResidue}` : 'No residue selected'}</span></footer>
    </main>
  );
}