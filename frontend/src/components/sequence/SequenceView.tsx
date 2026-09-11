import { useEffect, useState } from 'react';
import type { Chain } from '../../types/protein';

interface Props {
  chain: Chain | undefined;
  selectedResidue: number | null;
  onSelectResidue: (id: number) => void;
}

export function SequenceView({ chain, selectedResidue, onSelectResidue }: Props) {
  const pageSize = 50;
  const [page, setPage] = useState(0);
  useEffect(() => {
    setPage(0);
  }, [chain?.id]);

  const selectedIndex = chain?.residues.findIndex((residue) => residue.id === selectedResidue) ?? -1;
  useEffect(() => {
    if (selectedIndex >= 0) setPage(Math.floor(selectedIndex / pageSize));
  }, [selectedIndex]);

  if (!chain) return <div className="empty-state">No sequence loaded.</div>;
  const pageCount = Math.max(1, Math.ceil(chain.residues.length / pageSize));
  const start = page * pageSize;
  const residues = chain.residues.slice(start, start + pageSize);

  return (
    <div className="sequence-dropdown-view">
      <div className="sequence-page-toolbar">
        <button type="button" disabled={page === 0} onClick={() => setPage((current) => current - 1)}>Previous</button>
        <span>Residues {residues[0]?.id ?? 0}-{residues[residues.length - 1]?.id ?? 0} <b>·</b> page {page + 1} / {pageCount}</span>
        <button type="button" disabled={page >= pageCount - 1} onClick={() => setPage((current) => current + 1)}>Next</button>
      </div>
      <div className="sequence-view">
      <div className="sequence-ruler" aria-hidden="true">
        {residues.map((residue, index) => index % 10 === 0 ? <span key={residue.id}>{residue.id}</span> : <span key={residue.id} />)}
      </div>
      <div className="sequence-grid">
      {residues.map((residue, index) => (
        <button
          className={selectedResidue === residue.id ? 'residue selected' : 'residue'}
          key={residue.id}
          onClick={() => onSelectResidue(residue.id)}
          title={`${chain.id}:${residue.id} ${residue.name}`}
          aria-label={`${chain.id}:${residue.id} ${residue.name}`}
        >
          <span>{chain.sequence[start + index] ?? 'X'}</span>
          <small>{residue.id}</small>
        </button>
      ))}
      </div>
      </div>
    </div>
  );
}