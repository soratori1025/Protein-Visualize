import { useState } from 'react';
import { useProtein } from '../../contexts/ProteinContext';

interface HeaderProps {
  title: string;
  subtitle: string;
  eyebrow?: string;
  showHealth?: boolean;
}

export function Header({ title, subtitle, eyebrow = 'STRUCTURE LAB / MVP 0.1'}: HeaderProps) {
  const { handleUpload, handleClear, isUploading, protein, setStatus } = useProtein();
  const [pdbInput, setPdbInput] = useState('');
  const [isFetching, setIsFetching] = useState(false);

  const handleFetch = async () => {
    const acc = pdbInput.trim().toUpperCase();
    if (!acc) return;
    setIsFetching(true);
    setStatus(`Fetching ${acc} from RCSB PDB...`);
    try {
      const cifUrl = `https://files.rcsb.org/download/${acc}.cif`;
      const response = await fetch(cifUrl);
      if (!response.ok) {
        throw new Error(`PDB structure for ${acc} not found (Status: ${response.status})`);
      }
      const blob = await response.blob();
      const file = new File([blob], `${acc}.cif`);
      await handleUpload(file);
      setPdbInput('');
    } catch (err) {
      setStatus(err instanceof Error ? err.message : 'Fetch failed');
    } finally {
      setIsFetching(false);
    }
  };

  return (
    <header className="topbar">
      <div>
        <div className="eyebrow">{eyebrow}</div>
        <h1 style={{ fontSize: title.includes('Storyboard') ? '32px' : undefined, margin: title.includes('Storyboard') ? '0 0 16px' : undefined, color: title.includes('Storyboard') ? '#fff' : undefined }}>{title}</h1>
        <p style={{ color: title.includes('Storyboard') ? '#a5b9cb' : undefined, fontSize: title.includes('Storyboard') ? '14px' : undefined }}>{subtitle}</p>
      </div>
      <div className="header-actions" style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
        <div style={{ display: 'flex', alignItems: 'center', border: '1px solid #475569', borderRadius: '7px', overflow: 'hidden', background: '#081218' }}>
          <input
            type="text"
            placeholder="PDB ID"
            value={pdbInput}
            onChange={(e) => setPdbInput(e.target.value)}
            onKeyDown={(e) => e.key === 'Enter' && handleFetch()}
            style={{ background: 'transparent', border: 'none', color: '#fff', padding: '10px 14px', width: '180px', outline: 'none' }}
          />
          <button 
            className="upload-button" 
            style={{ borderRadius: '0', border: 'none', borderLeft: '1px solid #475569' }}
            onClick={handleFetch}
            disabled={isFetching || isUploading}
          >
            {isFetching ? 'Fetching...' : 'Fetch PDB'}
          </button>
        </div>
        {protein && (
          <button 
            className="upload-button" 
            style={{ background: 'transparent', border: '1px solid #475569', color: '#94a3b8' }}
            onClick={handleClear}
          >
            Clear Data
          </button>
        )}
        <label className="upload-button">
          <span>{isUploading && !isFetching ? 'Uploading...' : 'Upload File'}</span>
          <input 
            type="file" 
            accept=".pdb,.ent,.cif,.mmcif" 
            onChange={(event) => event.target.files?.[0] && handleUpload(event.target.files[0])} 
          />
        </label>
      </div>
    </header>
  );
}
