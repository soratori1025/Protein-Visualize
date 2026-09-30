import React, { useState } from 'react';
import { SequenceView } from '../sequence/SequenceView';
import { ProteinViewer } from '../viewer/ProteinViewer';
import type { Chain } from '../../types/protein';

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
    <section className="focused-workspace panel">
      <div className="panel-heading">
        <div>
          <span className="section-kicker">INTERACTION</span>
          <h2>Focused chain workspace</h2>
          <p className="panel-subtitle">
            Choose one chain, read its sequence, then inspect only that chain in 3D.
          </p>
        </div>
        <span className="tag">{chain ? `CHAIN ${chain.id}` : 'WAITING'}</span>
      </div>
      <div className="focused-sequence-layout">
        <div className="focused-sequence-column">
          <button
            className={`sequence-disclosure ${sequenceOpen ? 'open' : ''}`}
            onClick={() => setSequenceOpen((open) => !open)}
            aria-expanded={sequenceOpen}
          >
            <span>
              <span className="section-kicker">SELECTED SEQUENCE</span>
              <strong>Chain {chain?.id ?? '-'}</strong>
              <small>
                {chain?.residue_count ?? 0} residues · {sequenceOpen ? 'Hide sequence' : 'Show sequence'}
              </small>
            </span>
            <span className="disclosure-chevron">{sequenceOpen ? '▲' : '▼'}</span>
          </button>
          {sequenceOpen && (
            <div className="sequence-dropdown-content">
              <div className={`selection-inspector ${selectedResidueData ? 'active' : ''}`}>
                <span className="selection-label">CURRENT RESIDUE</span>
                <strong>
                  {selectedResidueData
                    ? `${chain?.id}:${selectedResidueData.id} · ${selectedResidueData.name}`
                    : 'Select a residue below'}
                </strong>
                <small>
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
        <div className="focused-viewer-column">
          <div className="panel-heading compact-heading">
            <div>
              <span className="section-kicker">3D ISOLATE</span>
              <h3>Chain {chain?.id ?? '-'}</h3>
            </div>
            <span className="tag">INTERACTIVE</span>
          </div>
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
      </div>
    </section>
  );
}
