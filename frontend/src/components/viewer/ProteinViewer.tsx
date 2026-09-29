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
  customHelixColors?: Record<string, string>;
  customRegionColors?: Record<string, string>;
  figureTheme?: 'lab' | 'publication';
}

const chainColors = [
  '#e63946', '#f4a261', '#2a9d8f', '#457b9d', '#9b5de5', '#f15bb5',
  '#8ab17d', '#e9c46a', '#ef8354', '#00b4d8', '#c77dff', '#90be6d',
];

export function ProteinViewer({ chain, chains, filename, variant = 'interactive', focusChainId, selectedResidue, onSelectResidue, consensusMap, defaultColorScheme, customHelixColors, customRegionColors, figureTheme }: Props) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<$3Dmol.Viewer>();
  const [style, setStyle] = useState<RepresentationStyle>('ribbon');
  const [colorScheme, setColorScheme] = useState<ColorScheme>(defaultColorScheme || 'chain');
  const [showCustomColors, setShowCustomColors] = useState(false);
  const [layoutMode, setLayoutMode] = useState<'native' | 'aligned' | 'spread'>('native');
  const [showMembraneRegions, setShowMembraneRegions] = useState(false);
  const [orthographic, setOrthographic] = useState(false);

  useEffect(() => {
    if (defaultColorScheme) {
      setColorScheme(defaultColorScheme);
    }
  }, [defaultColorScheme]);

  useEffect(() => {
    if (customRegionColors) {
      setRegionColors(prev => ({
        ...prev,
        Membrane: customRegionColors.Membrane || prev.Membrane,
        Extracellular: customRegionColors.ExtracellularText || prev.Extracellular,
        Cytoplasmic: customRegionColors.CytoplasmicText || prev.Cytoplasmic
      }));
    }
  }, [customRegionColors]);
  const [ssColors, setSsColors] = useState({
    Helix: '#ff0080',
    Strand: '#ffc107',
    Coil: '#00d2d3'
  });
  const [helixColors, setHelixColors] = useState<Record<string, string>>({});
  const [regionColors, setRegionColors] = useState<Record<string, string>>({
    Extracellular: '#64748b',
    Membrane: '#ffa600',
    Cytoplasmic: '#64748b'
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

      if ((layoutMode === 'spread' || layoutMode === 'aligned') && consensusMap && consensusMap.length > 0) {
        structure = spreadStructure(structure, consensusMap, { mode: layoutMode === 'spread' ? 'unrolled' : 'aligned' });
      }

      viewer = $3Dmol.createViewer(containerRef.current, { antialias: true });
      viewerRef.current = viewer;
      viewer.setBackgroundColor(figureTheme === 'publication' ? '#ffffff' : '#0b151e');
      viewer.addModel(structure, filename.toLowerCase().endsWith('.pdb') || filename.toLowerCase().endsWith('.ent') ? 'pdb' : 'cif');

      if (typeof (viewer as any).setProjection === 'function') {
        (viewer as any).setProjection(orthographic ? 'orthographic' : 'perspective');
      }

      applyStyles(viewer, chains, focusChainId, style, colorScheme, chain, consensusMap, consensusColors, ssColors, helixColors, showMembraneRegions, layoutMode, regionColors);

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

  useEffect(() => {
    if (customHelixColors !== undefined) {
      setHelixColors(customHelixColors);
    }
  }, [customHelixColors]);

  // Re-apply style when style, colorScheme, or selectedResidue changes
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    applyStyles(viewer, chains, focusChainId, style, colorScheme, chain, consensusMap, consensusColors, ssColors, helixColors, showMembraneRegions, layoutMode, regionColors);

    if (variant === 'interactive' && selectedResidue !== null) {
      viewer.setStyle({ chain: chain?.id, resi: selectedResidue }, {
        cartoon: { color: '#ff6f61', opacity: 1 },
        stick: { colorscheme: 'Jmol', radius: 0.28 },
        sphere: { color: '#ff6f61', scale: 0.9 },
      });
    }
    viewer.render();
  }, [chain, chains, focusChainId, selectedResidue, style, colorScheme, variant, consensusMap, consensusColors, ssColors, helixColors, showMembraneRegions, layoutMode, regionColors]);

  useEffect(() => {
    if (viewerRef.current) {
      viewerRef.current.setBackgroundColor(figureTheme === 'publication' ? '#ffffff' : '#0b151e');
      viewerRef.current.render();
    }
  }, [figureTheme]);

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
      <div className={`viewer-stage molecular-stage ${figureTheme === 'publication' ? 'light-mode' : ''}`}>
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
              <>
                <button
                  className={colorScheme === 'helices' ? 'toolbar-btn active' : 'toolbar-btn'}
                  onClick={() => setColorScheme('helices')}
                >
                  Individual Helices
                </button>
              </>
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
                  setOrthographic(true);
                }}
              >
                Spread Out (2D)
              </button>
            </div>
          )}
          {consensusMap && consensusMap.length > 0 && (layoutMode === 'spread' || layoutMode === 'aligned') && (
            <div className="toolbar-group">
              <span className="toolbar-label">REGIONS:</span>
              <button
                className={showMembraneRegions ? 'toolbar-btn active' : 'toolbar-btn'}
                onClick={() => setShowMembraneRegions(!showMembraneRegions)}
              >
                Show Membrane Regions
              </button>
            </div>
          )}
          <div className="toolbar-group">
            <span className="toolbar-label">EXPORT:</span>
            <button className="toolbar-btn" onClick={handleExportPNG}>
              Export PNG
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

              {colorScheme === 'helices' && consensusMap && (
                <div style={{ display: 'flex', gap: '16px', flexWrap: 'wrap' }}>
                  {Array.from({ length: new Set(consensusMap.map(r => r.crossing).filter(c => c != null)).size }).map((_, i) => {
                    const cNum = i + 1;
                    const defaultHelixColor = chainColors[(cNum - 1) % chainColors.length];
                    const currentColor = helixColors[cNum] || defaultHelixColor;
                    return (
                      <label key={cNum} style={{ display: 'flex', alignItems: 'center', gap: '6px', color: '#e7edf4', fontSize: '11px', fontWeight: 600, cursor: 'pointer' }}>
                        <input
                          type="color"
                          value={currentColor}
                          onChange={(e) => setHelixColors({ ...helixColors, [cNum]: e.target.value })}
                          style={{ border: 'none', padding: 0, width: '18px', height: '18px', cursor: 'pointer', background: 'transparent', borderRadius: '4px' }}
                          title={`Color for Helix ${cNum}`}
                        />
                        Helix {cNum}
                      </label>
                    );
                  })}
                </div>
              )}

              {(colorScheme === 'chain' || colorScheme === 'hydropathy') && (
                <div style={{ color: '#829ba9', fontSize: '11px', fontStyle: 'italic' }}>
                  Select <b>Secondary Structure</b>, <b>TM & SS Mapping</b>, or <b>Individual Helices</b> to customize colors.
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
    <div className={`viewer-stage ${figureTheme === 'publication' ? 'light-mode' : ''}`}>
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
  ssColors?: Record<string, string>,
  helixColors?: Record<string, string>,
  showMembraneRegions?: boolean,
  layoutMode?: 'native' | 'aligned' | 'spread',
  regionColors?: Record<string, string>
) {
  (viewer as any).removeAllShapes();
  if (typeof (viewer as any).removeAllLabels === 'function') {
    (viewer as any).removeAllLabels();
  }
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

  // Apply helices colors if selected
  if (colorScheme === 'helices' && selectedChain && consensusMap) {
    const defaultColor = '#64748b';

    // Find all crossings and their start/end residues
    const crossings = new Map<number, { start: number, end: number }>();
    consensusMap.forEach(r => {
      if (r.crossing) {
        if (!crossings.has(r.crossing)) {
          crossings.set(r.crossing, { start: r.residue_number, end: r.residue_number });
        } else {
          const c = crossings.get(r.crossing)!;
          c.start = Math.min(c.start, r.residue_number);
          c.end = Math.max(c.end, r.residue_number);
        }
      }
    });

    const crossingList = Array.from(crossings.entries()).sort((a, b) => a[0] - b[0]);
    const getHelixColor = (cNum: number) => {
      if (helixColors && helixColors[cNum]) return helixColors[cNum];
      return chainColors[(cNum - 1) % chainColors.length];
    };

    selectedChain.residues.forEach((res) => {
      let color = defaultColor;

      if (crossingList.length > 0) {
        let inCrossing = false;
        for (let i = 0; i < crossingList.length; i++) {
          const [cNum, bounds] = crossingList[i];
          if (res.id >= bounds.start && res.id <= bounds.end) {
            color = getHelixColor(cNum);
            inCrossing = true;
            break;
          }
        }

        if (!inCrossing) {
          let prev = null;
          let next = null;
          for (let i = 0; i < crossingList.length; i++) {
            const [, bounds] = crossingList[i];
            if (bounds.end < res.id) prev = crossingList[i];
            if (bounds.start > res.id && !next) next = crossingList[i];
          }

          if (prev && next) {
            const [prevNum, prevBounds] = prev;
            const [nextNum, nextBounds] = next;
            const prevColor = getHelixColor(prevNum);
            const nextColor = getHelixColor(nextNum);
            const fraction = (res.id - prevBounds.end) / (nextBounds.start - prevBounds.end);
            color = interpolateColor(prevColor, nextColor, fraction);
          } else if (prev) {
            color = getHelixColor(prev[0]);
          } else if (next) {
            color = getHelixColor(next[0]);
          }
        }
      }

      const sel = { chain: selectedChain.id, resi: res.id };
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
  // Draw Membrane Regions if toggled and aligned/spread
  if (showMembraneRegions && (layoutMode === 'spread' || layoutMode === 'aligned')) {
    const m = (viewer as any).getModel(0);
    if (m) {
      const atoms = m.selectedAtoms({});
      let minX = Infinity, maxX = -Infinity;
      let minZ = Infinity, maxZ = -Infinity;
      let maxY = -Infinity, minY = Infinity;
      atoms.forEach((a: any) => {
        if (a.x < minX) minX = a.x;
        if (a.x > maxX) maxX = a.x;
        if (a.y < minY) minY = a.y;
        if (a.y > maxY) maxY = a.y;
        if (a.z < minZ) minZ = a.z;
        if (a.z > maxZ) maxZ = a.z;
      });

      const pad = 30;
      const w = (maxX - minX) + pad;
      const d = (maxZ - minZ) + pad;
      const cx = (minX + maxX) / 2;
      const cz = (minZ + maxZ) / 2;

      const halfThickness = 15;

      // Membrane Core
      if (typeof (viewer as any).addBox === 'function') {
        (viewer as any).addBox({
          center: { x: cx, y: 0, z: cz },
          dimensions: { w, h: halfThickness * 2, d },
          color: regionColors?.Membrane || '#95a5a6',
          alpha: 0.5,
        });

        // Extracellular
        if (maxY > halfThickness) {
          const hTop = Math.max(10, maxY - halfThickness + pad / 2);
          const cyTop = halfThickness + hTop / 2;
          (viewer as any).addBox({
            center: { x: cx, y: cyTop, z: cz },
            dimensions: { w, h: hTop, d },
            color: regionColors?.Extracellular || '#ff6b6b',
            alpha: 0.1,
          });
          if (typeof (viewer as any).addLabel === 'function') {
            (viewer as any).addLabel("Extracellular", {
              position: { x: maxX + pad / 2, y: halfThickness + 10, z: cz },
              backgroundColor: 'transparent',
              fontColor: regionColors?.Extracellular || '#ff6b6b',
              backgroundOpacity: 0.0,
              fontSize: 14,
              alignment: 'center'
            });
          }
        }

        // Cytoplasmic
        if (minY < -halfThickness) {
          const hBot = Math.max(10, -halfThickness - minY + pad / 2);
          const cyBot = -halfThickness - hBot / 2;
          (viewer as any).addBox({
            center: { x: cx, y: cyBot, z: cz },
            dimensions: { w, h: hBot, d },
            color: regionColors?.Cytoplasmic || '#c44569',
            alpha: 0.1,
          });
          if (typeof (viewer as any).addLabel === 'function') {
            (viewer as any).addLabel("Cytoplasmic", {
              position: { x: maxX + pad / 2, y: -halfThickness - 10, z: cz },
              backgroundColor: 'transparent',
              fontColor: regionColors?.Cytoplasmic || '#c44569',
              backgroundOpacity: 0.0,
              fontSize: 14,
              alignment: 'center'
            });
          }
        }
      }
    }
  }

  viewer.render();
}

function interpolateColor(c1: string, c2: string, fraction: number): string {
  if (c1.startsWith('#')) c1 = c1.slice(1);
  if (c2.startsWith('#')) c2 = c2.slice(1);
  if (c1.length === 3) c1 = c1.split('').map(c => c + c).join('');
  if (c2.length === 3) c2 = c2.split('').map(c => c + c).join('');

  const r1 = parseInt(c1.slice(0, 2), 16);
  const g1 = parseInt(c1.slice(2, 4), 16);
  const b1 = parseInt(c1.slice(4, 6), 16);
  const r2 = parseInt(c2.slice(0, 2), 16);
  const g2 = parseInt(c2.slice(2, 4), 16);
  const b2 = parseInt(c2.slice(4, 6), 16);

  const r = Math.round(r1 + (r2 - r1) * fraction);
  const g = Math.round(g1 + (g2 - g1) * fraction);
  const b = Math.round(b1 + (b2 - b1) * fraction);

  return `#${(1 << 24 | r << 16 | g << 8 | b).toString(16).slice(1)}`;
}
