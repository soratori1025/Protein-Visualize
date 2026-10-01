import { useState } from 'react';
import { useProtein } from '../../contexts/ProteinContext';

interface HeaderProps {
  title: string;
  subtitle: string;
  eyebrow?: string;
  showHealth?: boolean;
}

export function Header({ title, subtitle, eyebrow = 'STRUCTURE LAB / MVP 0.1'}: HeaderProps) {
  const { 
    handleUpload, 
    handleFetchRemote, 
    handleClear, 
    isUploading, 
    protein, 
    setStatus,
    fetchInputId: inputId,
    setFetchInputId: setInputId,
    fetchInputType: inputType,
    setFetchInputType: setInputType
  } = useProtein();

  const handleFetch = async () => {
    const acc = inputId.trim().toUpperCase();
    if (!acc) return;
    // Calls the backend endpoint to handle UniProt resolution, PDB fetching, and uploading internally
    await handleFetchRemote(acc, inputType);
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
          <select 
            value={inputType} 
            onChange={e => setInputType(e.target.value as 'PDB' | 'UniProt')}
            style={{ background: '#1e293b', border: 'none', color: '#fff', padding: '10px', outline: 'none', borderRight: '1px solid #475569', cursor: 'pointer', fontSize: '14px' }}
          >
            <option value="PDB">PDB ID</option>
            <option value="UniProt">UniProt ID</option>
          </select>
          <input
            type="text"
            value={inputId}
            onChange={(e) => setInputId(e.target.value)}
            onKeyDown={(e) => e.key === 'Enter' && handleFetch()}
            style={{ background: 'transparent', border: 'none', color: '#fff', padding: '10px 14px', width: '130px', outline: 'none', fontSize: '14px' }}
          />
          <button 
            className="upload-button" 
            style={{ borderRadius: '0', border: 'none', borderLeft: '1px solid #475569' }}
            onClick={handleFetch}
            disabled={isUploading}
          >
            {isUploading ? 'Fetching...' : 'Fetch'}
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
          <span>{isUploading ? 'Uploading...' : 'Upload File'}</span>
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
