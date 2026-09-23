import { useProtein } from '../../contexts/ProteinContext';

interface HeaderProps {
  title: string;
  subtitle: string;
  eyebrow?: string;
  showHealth?: boolean;
}

export function Header({ title, subtitle, eyebrow = 'STRUCTURE LAB / MVP 0.1'}: HeaderProps) {
  const { handleUpload, isUploading, protein } = useProtein();
  return (
    <header className="topbar">
      <div>
        <div className="eyebrow">{eyebrow}</div>
        <h1 style={{ fontSize: title.includes('Storyboard') ? '32px' : undefined, margin: title.includes('Storyboard') ? '0 0 16px' : undefined, color: title.includes('Storyboard') ? '#fff' : undefined }}>{title}</h1>
        <p style={{ color: title.includes('Storyboard') ? '#a5b9cb' : undefined, fontSize: title.includes('Storyboard') ? '14px' : undefined }}>{subtitle}</p>
      </div>
      <div className="header-actions" style={{ display: 'flex', alignItems: 'center', gap: '12px' }}>
        <label className="upload-button">
          <span>{isUploading ? 'Uploading...' : 'Upload Structure'}</span>
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
