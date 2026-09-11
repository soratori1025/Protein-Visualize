import { useEffect, useRef } from 'react';
import * as $3Dmol from '3dmol';
import type { Chain } from '../../types/protein';

interface Props {
  chain: Chain | undefined;
  chains: Chain[];
  filename: string | undefined;
  variant?: 'overview' | 'interactive';
  focusChainId?: string;
  selectedResidue: number | null;
  onSelectResidue: (id: number) => void;
}

const chainColors = [
  '#e63946', '#f4a261', '#2a9d8f', '#457b9d', '#9b5de5', '#f15bb5',
  '#8ab17d', '#e9c46a', '#ef8354', '#00b4d8', '#c77dff', '#90be6d',
];

export function ProteinViewer({ chain, chains, filename, variant = 'interactive', focusChainId, selectedResidue, onSelectResidue }: Props) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<$3Dmol.Viewer>();

  useEffect(() => {
    if (!containerRef.current || !filename) return;
    let viewer: $3Dmol.Viewer | undefined;
    let cancelled = false;

    const loadViewer = async () => {
      const response = await fetch(`http://127.0.0.1:8000/api/structure/file/${encodeURIComponent(filename)}`);
      if (!response.ok || cancelled || !containerRef.current) return;
      const structure = await response.text();
      viewer = $3Dmol.createViewer(containerRef.current, { antialias: true });
      viewerRef.current = viewer;
      viewer.setBackgroundColor('#0b151e');
      viewer.addModel(structure, filename.toLowerCase().endsWith('.pdb') || filename.toLowerCase().endsWith('.ent') ? 'pdb' : 'cif');
      applyChainStyles(viewer, chains, focusChainId);
      if (focusChainId) {
        viewer.setStyle({}, { cartoon: { hidden: true } });
        viewer.setStyle({ chain: focusChainId }, { cartoon: { color: chainColors[chains.findIndex((item) => item.id === focusChainId) % chainColors.length], opacity: 1 } });
        viewer.zoomTo({ chain: focusChainId });
      } else {
        viewer.setStyle({ hetflag: true }, { stick: { colorscheme: 'Jmol', radius: 0.18 } });
      }
      if (variant === 'interactive') {
        viewer.setClickable({}, true, (atom) => {
          if (atom.resi !== undefined) onSelectResidue(Number(atom.resi));
        });
      }
      if (!focusChainId) viewer.zoomTo();
      viewer.render();
    };

    void loadViewer();
    return () => {
      cancelled = true;
      viewer?.clear();
      viewerRef.current = undefined;
      if (containerRef.current) containerRef.current.replaceChildren();
    };
  }, [chains, filename, focusChainId, onSelectResidue, variant]);

  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    applyChainStyles(viewer, chains, focusChainId);
    if (focusChainId) {
      viewer.setStyle({}, { cartoon: { hidden: true } });
      viewer.setStyle({ chain: focusChainId }, { cartoon: { color: chainColors[chains.findIndex((item) => item.id === focusChainId) % chainColors.length], opacity: 1 } });
    } else {
      viewer.setStyle({ hetflag: true }, { stick: { colorscheme: 'Jmol', radius: 0.18 } });
    }
    if (variant === 'interactive' && selectedResidue !== null) {
      viewer.setStyle({ chain: chain?.id, resi: selectedResidue }, {
        cartoon: { color: '#ff6f61', opacity: 1 },
        stick: { colorscheme: 'Jmol', radius: 0.28 },
      });
    }
    viewer.render();
  }, [chain?.id, chains, focusChainId, selectedResidue, variant]);

  if (!chain) return <div className="empty-state">Upload a PDB/mmCIF structure to begin.</div>;

  if (filename) {
    return (
      <div className="viewer-stage molecular-stage">
        <div ref={containerRef} className="molecular-viewer" aria-label={`Interactive 3D structure of ${filename}`} />
        <span className="viewer-caption">{focusChainId ? `Isolated chain ${focusChainId} · drag to rotate · scroll to zoom` : variant === 'overview' ? 'Full biological assembly · colored by chain' : 'Drag to rotate · scroll to zoom · click a residue'}</span>
      </div>
    );
  }

  const points = chain.residues.map((residue) => {
    const atom = residue.atoms.find((item) => item.name === 'CA') ?? residue.atoms[0];
    return atom ? { id: residue.id, x: atom.x, y: atom.y, z: atom.z } : null;
  }).filter((point): point is { id: number; x: number; y: number; z: number } => point !== null);
  const max = Math.max(...points.flatMap((point) => [Math.abs(point.x), Math.abs(point.y), Math.abs(point.z)]), 1);
  const project = (value: number, depth: number) => 50 + ((value + depth * 0.35) / (max * 2.7)) * 100;

  return (
    <div className="viewer-stage">
      <svg viewBox="0 0 200 150" role="img" aria-label={`3D projection of chain ${chain.id}`}>
        <polyline
          className="backbone-line"
          points={points.map((point) => `${project(point.x, point.z)},${project(point.y, -point.z)}`).join(' ')}
        />
        {points.map((point) => (
          <circle
            key={point.id}
            className={selectedResidue === point.id ? 'atom-point selected' : 'atom-point'}
            cx={project(point.x, point.z)}
            cy={project(point.y, -point.z)}
            r={selectedResidue === point.id ? 2.3 : 1.2}
            onClick={() => onSelectResidue(point.id)}
          />
        ))}
      </svg>
      <span className="viewer-caption">CA backbone projection · chain {chain.id}</span>
    </div>
  );
}

function applyChainStyles(viewer: $3Dmol.Viewer, chains: Chain[], focusChainId?: string) {
  viewer.setStyle({}, { cartoon: { color: '#d8dee9', opacity: 1 } });
  chains.forEach((item, index) => {
    if (focusChainId && item.id !== focusChainId) return;
    viewer.setStyle(
      { chain: item.id },
      { cartoon: { color: chainColors[index % chainColors.length], opacity: 1 } },
    );
  });
}