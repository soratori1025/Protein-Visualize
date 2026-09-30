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
  const [tmAlgorithm, setTmAlgorithm] = useState<string>('uniprot_api');
  const [topologySource, setTopologySource] = useState<'uniprot' | 'calculated'>('calculated');
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

  useEffect(() => {
    if (protein?.filename) {
      runTmAnalysis();
    }
  }, [protein?.filename]);

  return (
    <main className="app-shell">
      <Header 
        title="Transmembrane Analysis" 
        subtitle="Inspect secondary structure and transmembrane topology annotations." 
      />

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