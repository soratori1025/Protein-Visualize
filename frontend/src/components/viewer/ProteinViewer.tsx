import { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';
import type { Chain } from '../../types/protein';
import type { ColorScheme, RepresentationStyle } from '../../types/viewer';
import { generateRibbonSpline, getHydropathyColor } from '../../utils/ribbonSpline';
import { API_URL } from '../../services/api';
import type { ConsensusResidue } from '../../types/secondaryStructure';
import { spreadStructure } from '../structure/spreadStructure';

interface Props {
  chain: Chain | undefined;
  chains: Chain[];
  filename: string | undefined;
  variant?: 'overview' | 'interactive';
  focusChainId?: string;
  selectedResidue: number | null;
  onSelectResidue: (id: number) => void;
  consensusMap?: ConsensusResidue[];
  defaultColorScheme?: ColorScheme;
}

const chainColors = [
  '#e63946', '#f4a261', '#2a9d8f', '#457b9d', '#9b5de5', '#f15bb5',
  '#8ab17d', '#e9c46a', '#ef8354', '#00b4d8', '#c77dff', '#90be6d',
];

export function ProteinViewer({ chain, chains, filename, variant = 'interactive', focusChainId, selectedResidue, onSelectResidue, consensusMap, defaultColorScheme }: Props) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<$3Dmol.Viewer>();
  const [style, setStyle] = useState<RepresentationStyle>('ribbon');
  const [colorScheme, setColorScheme] = useState<ColorScheme>(defaultColorScheme || 'chain');
  const [showCustomColors, setShowCustomColors] = useState(false);
  const [layoutMode, setLayoutMode] = useState<'native' | 'spread'>('native');
  const [orthographic, setOrthographic] = useState(false);
  const [ssColors, setSsColors] = useState({
    Helix: '#ff0080',
    Strand: '#ffc107',
    Coil: '#00d2d3'
  });

  const [consensusColors, setConsensusColors] = useState({
    TM_E: '#ff9f43',
    TM_in: '#00d2d3',
    TM_C: '#5f27cd',
    Turn_E: '#ff6b6b',
    Turn_in: '#f368e0',
    Turn_C: '#c44569',
    Coil: '#64748b'
  });

  useEffect(() => {
    if (!containerRef.current || !filename) return;
    let viewer: $3Dmol.Viewer | undefined;
    let cancelled = false;

    const loadViewer = async () => {
      const response = await fetch(`${API_URL}/api/structure/file/${encodeURIComponent(filename)}`);
      if (!response.ok || cancelled || !containerRef.current) return;
      let structure = await response.text();

      if (layoutMode === 'spread' && consensusMap && consensusMap.length > 0) {
        structure = spreadStructure(structure, consensusMap);
      }

      viewer = $3Dmol.createViewer(containerRef.current, { antialias: true });
      viewerRef.current = viewer;
      viewer.setBackgroundColor('#0b151e');
      viewer.addModel(structure, filename.toLowerCase().endsWith('.pdb') || filename.toLowerCase().endsWith('.ent') ? 'pdb' : 'cif');

      if (typeof (viewer as any).setProjection === 'function') {
        (viewer as any).setProjection(orthographic ? 'orthographic' : 'perspective');
      }

      applyStyles(viewer, chains, focusChainId, style, colorScheme, chain, consensusMap, consensusColors);

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
  }, [chains, filename, focusChainId, onSelectResidue, variant, layoutMode, consensusMap]); // Intentionally omitting orthographic here, handled in the next effect

  // Handle orthographic toggle separately without reloading the model
  useEffect(() => {
    if (viewerRef.current && typeof (viewerRef.current as any).setProjection === 'function') {
      (viewerRef.current as any).setProjection(orthographic ? 'orthographic' : 'perspective');
      viewerRef.current.render();
    }
  }, [orthographic]);

  // Re-apply style when style, colorScheme, or selectedResidue changes
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    applyStyles(viewer, chains, focusChainId, style, colorScheme, chain, consensusMap, consensusColors, ssColors);

    if (variant === 'interactive' && selectedResidue !== null) {
      viewer.setStyle({ chain: chain?.id, resi: selectedResidue }, {
        cartoon: { color: '#ff6f61', opacity: 1 },
        stick: { colorscheme: 'Jmol', radius: 0.28 },
        sphere: { color: '#ff6f61', scale: 0.9 },
      });
    }
    viewer.render();
  }, [chain, chains, focusChainId, selectedResidue, style, colorScheme, variant, consensusMap, consensusColors]);

  const handleExportPNG = () => {
    if (viewerRef.current) {
      const imgURI = (viewerRef.current as any).pngURI();
      const a = document.createElement('a');
      a.href = imgURI;
      a.download = `${filename || 'structure'}_snapshot.png`;
      a.click();
    }
  };

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
            <span className="toolbar-label">ANALYSIS METHOD:</span>
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
            {consensusMap && (
              <button
                className={colorScheme === 'consensus' ? 'toolbar-btn active' : 'toolbar-btn'}
                onClick={() => setColorScheme('consensus')}
              >
                TM & SS Mapping
              </button>
            )}
            <button
              className={colorScheme === 'hydropathy' ? 'toolbar-btn active' : 'toolbar-btn'}
              onClick={() => setColorScheme('hydropathy')}
            >
              Hydropathy
            </button>
          </div>
          {consensusMap && consensusMap.length > 0 && (
            <div className="toolbar-group">
              <span className="toolbar-label">LAYOUT:</span>
              <button
                className={layoutMode === 'native' ? 'toolbar-btn active' : 'toolbar-btn'}
                onClick={() => setLayoutMode('native')}
              >
                Native 3D
              </button>
              <button
                className={layoutMode === 'spread' ? 'toolbar-btn active' : 'toolbar-btn'}
                onClick={() => {
                  setLayoutMode('spread');
                  setOrthographic(true); // Default to orthographic for 2D look
                }}
              >
                Spread Out (2D-like)
              </button>
            </div>
          )}
          <div className="toolbar-group">
            <span className="toolbar-label">VIEW:</span>
            <button
              className={orthographic ? 'toolbar-btn active' : 'toolbar-btn'}
              onClick={() => setOrthographic(!orthographic)}
            >
              2D Projection
            </button>
            <button className="toolbar-btn" onClick={handleExportPNG}>
              Export PNG
            </button>
            <button
              className={showCustomColors ? 'toolbar-btn active' : 'toolbar-btn'}
              onClick={() => setShowCustomColors(!showCustomColors)}
              style={{ marginLeft: '8px', border: '1px solid #78d8c1' }}
            >
              {showCustomColors ? 'Hide Colors' : 'Customize Colors'}
            </button>
          </div>

          {showCustomColors && (
            <div style={{ flexBasis: '100%', display: 'flex', flexWrap: 'wrap', gap: '16px', alignItems: 'center', marginTop: '8px', paddingTop: '8px', borderTop: '1px dashed rgba(120, 216, 193, 0.3)' }}>
              <span className="toolbar-label" style={{ color: '#78d8c1' }}>CUSTOM COLORS:</span>

              {colorScheme === 'ss' && (
                <div style={{ display: 'flex', gap: '16px', flexWrap: 'wrap' }}>
                  {Object.entries(ssColors).map(([key, value]) => (
                    <label key={key} style={{ display: 'flex', alignItems: 'center', gap: '6px', color: '#e7edf4', fontSize: '11px', fontWeight: 600, cursor: 'pointer' }}>
                      <input
                        type="color"
                        value={value}
                        onChange={(e) => setSsColors({ ...ssColors, [key]: e.target.value })}
                        style={{ border: 'none', padding: 0, width: '18px', height: '18px', cursor: 'pointer', background: 'transparent', borderRadius: '4px' }}
                        title={`Color for ${key}`}
                      />
                      {key}
                    </label>
                  ))}
                </div>
              )}

              {colorScheme === 'consensus' && (
                <div style={{ display: 'flex', gap: '16px', flexWrap: 'wrap' }}>
                  {Object.entries(consensusColors).map(([key, value]) => (
                    <label key={key} style={{ display: 'flex', alignItems: 'center', gap: '6px', color: '#e7edf4', fontSize: '11px', fontWeight: 600, cursor: 'pointer' }}>
                      <input
                        type="color"
                        value={value}
                        onChange={(e) => setConsensusColors({ ...consensusColors, [key]: e.target.value })}
                        style={{ border: 'none', padding: 0, width: '18px', height: '18px', cursor: 'pointer', background: 'transparent', borderRadius: '4px' }}
                        title={`Color for ${key}`}
                      />
                      {key.replace('_', ' ')}
                    </label>
                  ))}
                </div>
              )}

              {(colorScheme === 'chain' || colorScheme === 'hydropathy') && (
                <div style={{ color: '#829ba9', fontSize: '11px', fontStyle: 'italic' }}>
                  Select <b>Secondary Structure</b> or <b>TM & SS Mapping</b> algorithm to customize structural colors.
                </div>
              )}
            </div>
          )}
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
  selectedChain?: Chain,
  consensusMap?: ConsensusResidue[],
  consensusColors?: Record<string, string>,
  ssColors?: Record<string, string>
) {
  (viewer as any).removeAllShapes();
  viewer.setStyle({}, { cartoon: { hidden: true }, stick: { hidden: true }, sphere: { hidden: true }, line: { hidden: true } });

  chains.forEach((item, index) => {
    if (focusChainId && item.id !== focusChainId) return;

    const baseColor = chainColors[index % chainColors.length];
    const selection = focusChainId ? { chain: focusChainId } : { chain: item.id };

    if (style === 'ribbon') {
      if (colorScheme === 'ss') {
        const colorfunc = (atom: any) => {
          if (atom.ss === 'h') return ssColors?.Helix || '#ff0080';
          if (atom.ss === 's') return ssColors?.Strand || '#ffc107';
          return ssColors?.Coil || '#00d2d3';
        };
        viewer.setStyle(selection, { cartoon: { colorfunc, opacity: 1 } });
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

  // Apply consensus colors if selected
  if (colorScheme === 'consensus' && selectedChain && consensusMap) {
    // Default color for residues not in the membrane
    const defaultColor = '#64748b';
    const selAll = { chain: selectedChain.id };

    if (style === 'ribbon') {
      viewer.setStyle(selAll, { cartoon: { color: defaultColor, opacity: 1 } });
    } else if (style === 'sphere') {
      viewer.setStyle(selAll, { sphere: { color: defaultColor, scale: 0.75 } });
    } else if (style === 'stick') {
      viewer.setStyle(selAll, { stick: { color: defaultColor, radius: 0.22 } });
    } else if (style === 'line') {
      viewer.setStyle(selAll, { line: { color: defaultColor, linewidth: 1.5 } });
    }

    consensusMap.forEach((res) => {
      let color = consensusColors?.Coil || defaultColor;
      if (res.label === 'TM_E') color = consensusColors?.TM_E || '#ff9f43';
      if (res.label === 'TM_in') color = consensusColors?.TM_in || '#00d2d3';
      if (res.label === 'TM_C') color = consensusColors?.TM_C || '#5f27cd';
      if (res.label === 'Turn_E') color = consensusColors?.Turn_E || '#ff6b6b';
      if (res.label === 'Turn_in') color = consensusColors?.Turn_in || '#f368e0';
      if (res.label === 'Turn_C') color = consensusColors?.Turn_C || '#c44569';

      const sel = { chain: selectedChain.id, resi: res.residue_number };
      if (style === 'ribbon') {
        viewer.setStyle(sel, { cartoon: { color, opacity: 1 } });
      } else if (style === 'sphere') {
        viewer.setStyle(sel, { sphere: { color, scale: 0.75 } });
      } else if (style === 'stick') {
        viewer.setStyle(sel, { stick: { color, radius: 0.22 } });
      } else if (style === 'line') {
        viewer.setStyle(sel, { line: { color, linewidth: 1.5 } });
      }
    });
  }

  // Always keep heteratoms (ligands, water) readable as sticks
  viewer.setStyle({ hetflag: true }, { stick: { colorscheme: 'Jmol', radius: 0.18 } });
}

