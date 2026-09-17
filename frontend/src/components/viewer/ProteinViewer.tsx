import { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';
import type { Chain } from '../../types/protein';
import type { ColorScheme, RepresentationStyle } from '../../types/viewer';
import { generateRibbonSpline, getHydropathyColor } from '../../utils/ribbonSpline';

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
  const [style, setStyle] = useState<RepresentationStyle>('ribbon');
  const [colorScheme, setColorScheme] = useState<ColorScheme>('chain');

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

      applyStyles(viewer, chains, focusChainId, style, colorScheme, chain);

      if (variant === 'interactive') {
        viewer.setClickable({}, true, (atom) => {
          if (atom.resi !== undefined) onSelectResidue(Number(atom.resi));
        });
      }
      if (focusChainId) {
        viewer.zoomTo({ chain: focusChainId });
      } else {
        viewer.zoomTo();
      }
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

  // Re-apply style when style, colorScheme, or selectedResidue changes
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    applyStyles(viewer, chains, focusChainId, style, colorScheme, chain);

    if (variant === 'interactive' && selectedResidue !== null) {
      viewer.setStyle({ chain: chain?.id, resi: selectedResidue }, {
        cartoon: { color: '#ff6f61', opacity: 1 },
        stick: { colorscheme: 'Jmol', radius: 0.28 },
        sphere: { color: '#ff6f61', scale: 0.9 },
      });
    }
    viewer.render();
  }, [chain, chains, focusChainId, selectedResidue, style, colorScheme, variant]);

  if (!chain) return <div className="empty-state">Upload a PDB/mmCIF structure to begin.</div>;

  if (filename) {
    return (
      <div className="viewer-stage molecular-stage">
        <div className="viewer-toolbar">
          <div className="toolbar-group">
            <span className="toolbar-label">STYLE:</span>
            {(['ribbon', 'stick', 'sphere', 'line', 'pipesAndPlanks'] as RepresentationStyle[]).map((st) => (
              <button
                key={st}
                className={style === st ? 'toolbar-btn active' : 'toolbar-btn'}
                onClick={() => setStyle(st)}
              >
                {st === 'ribbon' ? 'Ribbon' : st === 'stick' ? 'Stick' : st === 'sphere' ? 'Sphere' : st === 'line' ? 'Line' : 'Pipes & Planks'}
              </button>
            ))}
          </div>

          <div className="toolbar-group">
            <span className="toolbar-label">COLOR:</span>
            <button
              className={colorScheme === 'chain' ? 'toolbar-btn active' : 'toolbar-btn'}
              onClick={() => setColorScheme('chain')}
            >
              Chain
            </button>
            <button
              className={colorScheme === 'ss' ? 'toolbar-btn active' : 'toolbar-btn'}
              onClick={() => setColorScheme('ss')}
            >
              Secondary Structure
            </button>
            <button
              className={colorScheme === 'hydropathy' ? 'toolbar-btn active' : 'toolbar-btn'}
              onClick={() => setColorScheme('hydropathy')}
            >
              Hydropathy
            </button>
          </div>
        </div>

        <div ref={containerRef} className="molecular-viewer" aria-label={`Interactive 3D structure of ${filename}`} />
        <span className="viewer-caption">
          {focusChainId ? `Isolated chain ${focusChainId} · style: ${style} · color: ${colorScheme}` : variant === 'overview' ? `Full biological assembly · style: ${style}` : `Drag to rotate · scroll to zoom · click a residue`}
        </span>
      </div>
    );
  }

  // Fallback 2D/3D SVG Backbone Projection using Catmull-Rom Ribbon Spline
  const caPoints = chain.residues.map((residue) => {
    const atom = residue.atoms.find((item) => item.name === 'CA') ?? residue.atoms[0];
    return atom ? { id: residue.id, x: atom.x, y: atom.y, z: atom.z, name: residue.name } : null;
  }).filter((point): point is { id: number; x: number; y: number; z: number; name: string } => point !== null);

  const max = Math.max(...caPoints.flatMap((point) => [Math.abs(point.x), Math.abs(point.y), Math.abs(point.z)]), 1);
  const project = (value: number, depth: number) => 50 + ((value + depth * 0.35) / (max * 2.7)) * 100;

  // Generate smooth Catmull-Rom spline curve points
  const ribbonSpline = generateRibbonSpline(caPoints, 4, 1.4);

  return (
    <div className="viewer-stage">
      <svg viewBox="0 0 200 150" role="img" aria-label={`Catmull-Rom ribbon projection of chain ${chain.id}`}>
        {/* Draw smooth Catmull-Rom Ribbon band */}
        {ribbonSpline.length > 1 && (
          <path
            className="ribbon-spline-path"
            d={
              `M ${ribbonSpline.map((p) => `${project(p.left.x, p.left.z)},${project(p.left.y, -p.left.z)}`).join(' L ')} ` +
              `L ${ribbonSpline.slice().reverse().map((p) => `${project(p.right.x, p.right.z)},${project(p.right.y, -p.right.z)}`).join(' L ')} Z`
            }
            fill="#78d8c1"
            fillOpacity="0.4"
            stroke="#78d8c1"
            strokeWidth="0.8"
          />
        )}

        {/* Backbone center line */}
        <polyline
          className="backbone-line"
          points={caPoints.map((point) => `${project(point.x, point.z)},${project(point.y, -point.z)}`).join(' ')}
          stroke="#457b9d"
          strokeWidth="0.5"
          strokeDasharray="1 1"
        />

        {caPoints.map((point) => (
          <circle
            key={point.id}
            className={selectedResidue === point.id ? 'atom-point selected' : 'atom-point'}
            cx={project(point.x, point.z)}
            cy={project(point.y, -point.z)}
            r={selectedResidue === point.id ? 2.3 : 1.2}
            fill={colorScheme === 'hydropathy' ? getHydropathyColor(point.name) : selectedResidue === point.id ? '#ff6f61' : '#78d8c1'}
            onClick={() => onSelectResidue(point.id)}
          />
        ))}
      </svg>
      <span className="viewer-caption">Catmull-Rom Ribbon Spline Projection · chain {chain.id}</span>
    </div>
  );
}

function applyStyles(
  viewer: $3Dmol.Viewer,
  chains: Chain[],
  focusChainId: string | undefined,
  style: RepresentationStyle,
  colorScheme: ColorScheme,
  selectedChain?: Chain
) {
  (viewer as any).removeAllShapes();
  viewer.setStyle({}, { cartoon: { hidden: true }, stick: { hidden: true }, sphere: { hidden: true }, line: { hidden: true } });

  chains.forEach((item, index) => {
    if (focusChainId && item.id !== focusChainId) return;

    const baseColor = chainColors[index % chainColors.length];
    const selection = focusChainId ? { chain: focusChainId } : { chain: item.id };

    if (style === 'ribbon') {
      if (colorScheme === 'ss') {
        viewer.setStyle(selection, { cartoon: { colorscheme: 'ssPyMOL', opacity: 1 } });
      } else {
        viewer.setStyle(selection, { cartoon: { color: baseColor, opacity: 1 } });
      }
    } else if (style === 'pipesAndPlanks') {
      import('../../utils/PipePlanks').then(({ drawCustomPipesAndPlanks }) => {
        drawCustomPipesAndPlanks(viewer, item.id, baseColor, colorScheme);
        viewer.render();
      });
    } else if (style === 'stick') {
      viewer.setStyle(selection, { stick: { colorscheme: 'Jmol', radius: 0.22 } });
    } else if (style === 'sphere') {
      viewer.setStyle(selection, { sphere: { colorscheme: 'Jmol', scale: 0.75 } });
    } else if (style === 'line') {
      viewer.setStyle(selection, { line: { colorscheme: 'Jmol', linewidth: 1.5 } });
    }
  });

  // Apply hydropathy colors if selected
  if (colorScheme === 'hydropathy' && selectedChain) {
    selectedChain.residues.forEach((res) => {
      const color = getHydropathyColor(res.name);
      const sel = { chain: selectedChain.id, resi: res.id };
      if (style === 'ribbon') {
        viewer.setStyle(sel, { cartoon: { color, opacity: 1 } });
      } else if (style === 'sphere') {
        viewer.setStyle(sel, { sphere: { color, scale: 0.75 } });
      } else if (style === 'stick') {
        viewer.setStyle(sel, { stick: { color, radius: 0.22 } });
      }
    });
  }

  // Always keep heteratoms (ligands, water) readable as sticks
  viewer.setStyle({ hetflag: true }, { stick: { colorscheme: 'Jmol', radius: 0.18 } });
}