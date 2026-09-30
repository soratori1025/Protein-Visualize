import { useState } from 'react';
import { ProteinViewer } from '../components/viewer/ProteinViewer';
import { FocusedChainWorkspace } from '../components/structure/FocusedChainWorkspace';
import { ExpandedProteinMap } from '../components/structure/ExpandedProteinMap';
import { AnalysisPanel } from '../components/analysis/AnalysisPanel';
import { analyzeChain } from '../services/api';

import { useProtein } from '../contexts/ProteinContext';
import { Header } from '../components/layout/Header';

export function StructureViewer() {
  const { protein, chainId, setChainId, selectedResidue, setSelectedResidue, status, setStatus, health, chainAnalysis, setChainAnalysis } = useProtein();
  const [analysisLoading, setAnalysisLoading] = useState(false);
  const chain = protein?.models[0]?.chains.find((item) => item.id === chainId) ?? protein?.models[0]?.chains[0];
  const chains = protein?.models[0]?.chains ?? [];
  const runChainAnalysis = async () => {
    if (!protein || !chain) {
      setStatus('Upload a structure and select a chain before analysis');
      return;
    }
    setAnalysisLoading(true);
    try {
      setChainAnalysis(await analyzeChain(protein.filename, chain.id));
      setStatus(`Analysis complete for chain ${chain.id}`);
    } catch (error) {
      setStatus(error instanceof Error ? error.message : 'Analysis failed');
    } finally {
      setAnalysisLoading(false);
    }
  };

  return (
    <main className="app-shell">
      <Header 
        title="Visualize Workspace" 
        subtitle="Inspect coordinates, sequence, and overall architecture." 
        showHealth={true} 
      />

      <section className="overview-section panel">
        <div className="panel-heading"><div><span className="section-kicker">FULL ASSEMBLY</span><h2>Biological assembly overview</h2><p className="panel-subtitle">The complete colored protein stays together here, like the reference image.</p></div><span className="tag">{chains.length ? `${chains.length} CHAINS` : 'WAITING'}</span></div>
        <ProteinViewer chain={chain} chains={chains} filename={protein?.filename} variant="overview" selectedResidue={selectedResidue} onSelectResidue={setSelectedResidue} />
      </section>

      <section className="spread-section panel">
        <div className="panel-heading"><div><span className="section-kicker">EXPANDED ARCHITECTURE</span><h2>Spread protein map</h2><p className="panel-subtitle">Compounds are separated into readable lanes while the external ribbons preserve their assembly relationships.</p></div></div>
        <ExpandedProteinMap chains={chains} selectedChain={chain?.id} onSelectChain={(nextChain) => { setChainId(nextChain); setSelectedResidue(null); setChainAnalysis(null); }} />
      </section>

      <FocusedChainWorkspace
        chain={chain}
        chains={chains}
        filename={protein?.filename}
        selectedResidue={selectedResidue}
        onSelectResidue={setSelectedResidue}
      />

      <section className="analysis-section panel"><AnalysisPanel analysis={chainAnalysis} loading={analysisLoading} onRun={runChainAnalysis} /></section>

      <footer className="status-bar"><span><b className="status-dot" />{status}</span><span>{health}</span><span>{selectedResidue ? `Selected residue ${chain?.id}:${selectedResidue}` : 'No residue selected'}</span></footer>
    </main>
  );
}
