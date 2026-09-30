import React, { useState } from 'react';
import { SequenceView } from '../sequence/SequenceView';
import { ProteinViewer } from '../viewer/ProteinViewer';
import type { Chain } from '../../types/protein';
import './FocusedChainWorkspace.css';

interface FocusedChainWorkspaceProps {
  chain?: Chain;
  chains: Chain[];
  filename?: string;
  selectedResidue: number | null;
  onSelectResidue: (resId: number | null) => void;
}

export function FocusedChainWorkspace({
  chain,
  chains,
  filename,
  selectedResidue,
  onSelectResidue,
}: FocusedChainWorkspaceProps) {
  const [sequenceOpen, setSequenceOpen] = useState(false);
  const selectedResidueData = chain?.residues.find((item) => item.id === selectedResidue);

  return (
    <section className="premium-workspace-section">
      <div className="premium-workspace-header">
        <div>
          <h2 className="premium-workspace-title">Focused chain workspace</h2>
          <p className="premium-workspace-subtitle">
            Choose one chain, read its sequence, then inspect only that chain in 3D.
          </p>
        </div>
        <span className="premium-tag">{chain ? `CHAIN ${chain.id}` : 'WAITING'}</span>
      </div>
      
      <div className="premium-workspace-layout">
        <div className="premium-sequence-column">
          <button
            className={`premium-sequence-btn ${sequenceOpen ? 'open' : ''}`}
            onClick={() => setSequenceOpen((open) => !open)}
            aria-expanded={sequenceOpen}
          >
            <span className="premium-btn-content">
              <span className="premium-btn-kicker">SELECTED SEQUENCE</span>
              <strong className="premium-btn-title">Chain {chain?.id ?? '-'}</strong>
              <small className="premium-btn-subtitle">
                {chain?.residue_count ?? 0} residues · {sequenceOpen ? 'Hide sequence' : 'Show sequence'}
              </small>
            </span>
            <span className="premium-chevron">{sequenceOpen ? '▲' : '▼'}</span>
          </button>
          <ProteinViewer
            chain={chain}
            chains={chains}
            filename={filename}
            variant="interactive"
            focusChainId={chain?.id}
            selectedResidue={selectedResidue}
            onSelectResidue={onSelectResidue}
          />
        </div>
        {sequenceOpen && (
            <div className="premium-dropdown-content">
              <div className={`premium-inspector ${selectedResidueData ? 'active' : ''}`}>
                <span className="inspector-label">CURRENT RESIDUE</span>
                <strong className="inspector-value">
                  {selectedResidueData
                    ? `${chain?.id}:${selectedResidueData.id} · ${selectedResidueData.name}`
                    : 'Select a residue below'}
                </strong>
                <small className="inspector-help">
                  {selectedResidueData
                    ? 'Highlighted in the isolated 3D chain and topology map.'
                    : 'Click an amino acid to highlight its exact position in 3D.'}
                </small>
              </div>
              <SequenceView
                chain={chain}
                selectedResidue={selectedResidue}
                onSelectResidue={onSelectResidue}
              />
            </div>
          )}
      </div>
    </section>
  );
}
