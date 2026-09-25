import { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';
import type { Chain } from '../../types/protein';
import type { ColorScheme, RepresentationStyle } from '../../types/viewer';
import { generateRibbonSpline, getHydropathyColor } from '../../utils/ribbonSpline';
import { API_URL } from '../../services/api';
import type { ConsensusResidue } from '../../types/secondaryStructure';

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
  const [layoutMode, setLayoutMode] = useState<'native' | 'spread'>('native');

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

      applyStyles(viewer, chains, focusChainId, style, colorScheme, chain, consensusMap);

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
  }, [chains, filename, focusChainId, onSelectResidue, variant, layoutMode, consensusMap]);

  // Re-apply style when style, colorScheme, or selectedResidue changes
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer) return;
    applyStyles(viewer, chains, focusChainId, style, colorScheme, chain, consensusMap);

    if (variant === 'interactive' && selectedResidue !== null) {
      viewer.setStyle({ chain: chain?.id, resi: selectedResidue }, {
        cartoon: { color: '#ff6f61', opacity: 1 },
        stick: { colorscheme: 'Jmol', radius: 0.28 },
        sphere: { color: '#ff6f61', scale: 0.9 },
      });
    }
    viewer.render();
  }, [chain, chains, focusChainId, selectedResidue, style, colorScheme, variant, consensusMap]);

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
            {consensusMap && (
              <button
                className={colorScheme === 'consensus' ? 'toolbar-btn active' : 'toolbar-btn'}
                onClick={() => setColorScheme('consensus')}
              >
                Consensus
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
                onClick={() => setLayoutMode('spread')}
              >
                Spread Out (2D-like)
              </button>
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
  consensusMap?: ConsensusResidue[]
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
      let color = defaultColor;
      if (res.label === 'TM_E') color = '#ff9f43';
      if (res.label === 'TM_in') color = '#00d2d3';
      if (res.label === 'TM_C') color = '#5f27cd';
      if (res.label === 'Turn_E') color = '#ff6b6b';
      if (res.label === 'Turn_in') color = '#f368e0';
      if (res.label === 'Turn_C') color = '#c44569';
      
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


function spreadStructure(pdbData: string, consensusMap: ConsensusResidue[]): string {
  // 1. Group residues into TM segments
  const tmSegments: { segment: number, start: number, end: number }[] = [];
  let currentSeg: any = null;
  
  for (const res of consensusMap) {
    if (res.tm_segment !== undefined && res.tm_segment !== null) {
      if (!currentSeg || currentSeg.segment !== res.tm_segment) {
        if (currentSeg) tmSegments.push(currentSeg);
        currentSeg = { segment: res.tm_segment, start: res.residue_number, end: res.residue_number };
      } else {
        currentSeg.end = res.residue_number;
      }
    } else {
      if (currentSeg) {
        tmSegments.push(currentSeg);
        currentSeg = null;
      }
    }
  }
  if (currentSeg) tmSegments.push(currentSeg);

  if (tmSegments.length === 0) return pdbData;

  // 2. Parse PDB atoms
  const lines = pdbData.split('\n');
  const atoms: any[] = [];
  const residues = new Map<number, { ca: any, atoms: any[] }>();

  for (let i = 0; i < lines.length; i++) {
    const line = lines[i];
    if (line.startsWith('ATOM  ') || line.startsWith('HETATM')) {
      const resSeq = parseInt(line.substring(22, 26).trim());
      if (isNaN(resSeq)) continue;
      
      const name = line.substring(12, 16).trim();
      const x = parseFloat(line.substring(30, 38));
      const y = parseFloat(line.substring(38, 46));
      const z = parseFloat(line.substring(46, 54));
      
      const atom = { index: i, line, resSeq, name, x, y, z, newX: x, newY: y, newZ: z };
      atoms.push(atom);
      
      if (!residues.has(resSeq)) residues.set(resSeq, { ca: null, atoms: [] });
      residues.get(resSeq)!.atoms.push(atom);
      if (name === 'CA') residues.get(resSeq)!.ca = atom;
    }
  }

  // Fallback CA if missing
  for (const [resSeq, res] of residues.entries()) {
    if (!res.ca && res.atoms.length > 0) res.ca = res.atoms[0];
  }

  // Math helpers
  const normalize = (v: any) => {
    const len = Math.sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    return len > 0 ? {x: v.x/len, y: v.y/len, z: v.z/len} : {x:0,y:0,z:0};
  };
  const cross = (a: any, b: any) => ({
    x: a.y*b.z - a.z*b.y,
    y: a.z*b.x - a.x*b.z,
    z: a.x*b.y - a.y*b.x
  });
  const dot = (a: any, b: any) => a.x*b.x + a.y*b.y + a.z*b.z;
  const rotate = (p: any, R: any) => ({
    x: p.x*R[0][0] + p.y*R[0][1] + p.z*R[0][2],
    y: p.x*R[1][0] + p.y*R[1][1] + p.z*R[1][2],
    z: p.x*R[2][0] + p.y*R[2][1] + p.z*R[2][2]
  });

  const getRotation = (a: any, b: any) => {
    const v = cross(a, b);
    const s = Math.sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
    const c = dot(a, b);
    if (s < 1e-6) {
      return c > 0 ? [[1,0,0],[0,1,0],[0,0,1]] : [[-1,0,0],[0,-1,0],[0,0,-1]];
    }
    const vx = [[0, -v.z, v.y], [v.z, 0, -v.x], [-v.y, v.x, 0]];
    const factor = (1 - c) / (s * s);
    const R = [[1,0,0],[0,1,0],[0,0,1]];
    for (let i=0; i<3; i++) {
      for (let j=0; j<3; j++) {
        R[i][j] += vx[i][j];
        let vx2 = 0;
        for (let k=0; k<3; k++) vx2 += vx[i][k] * vx[k][j];
        R[i][j] += vx2 * factor;
      }
    }
    return R;
  };

  // 3. Process TM segments
  const tmAnchors = new Map<number, { startCA: any, endCA: any }>();
  let currentX = 0;

  for (let i = 0; i < tmSegments.length; i++) {
    const seg = tmSegments[i];
    const firstRes = residues.get(seg.start);
    const lastRes = residues.get(seg.end);
    if (!firstRes || !lastRes) continue;

    if (i > 0) {
      const prevSeg = tmSegments[i-1];
      const numRes = seg.start - prevSeg.end; // number of bonds in loop
      const maxReach = numRes * 3.2; // tighter spacing to ensure no breaks
      const spacingX = Math.min(30, maxReach); // Cap horizontal spacing
      currentX += spacingX;
    }

    // Center of mass
    let cx=0, cy=0, cz=0, count=0;
    for (let r = seg.start; r <= seg.end; r++) {
      const rr = residues.get(r);
      if (rr) {
        for (const a of rr.atoms) { cx+=a.x; cy+=a.y; cz+=a.z; count++; }
      }
    }
    if(count > 0) { cx/=count; cy/=count; cz/=count; }

    const vOrig = normalize({
      x: lastRes.ca.x - firstRes.ca.x,
      y: lastRes.ca.y - firstRes.ca.y,
      z: lastRes.ca.z - firstRes.ca.z
    });
    
    // Determine orientation based on membrane topology.
    // We alternate. If first is IN->OUT, it points UP (0, 1, 0).
    const vTarget = i % 2 === 0 ? {x: 0, y: 1, z: 0} : {x: 0, y: -1, z: 0};
    const R = getRotation(vOrig, vTarget);
    
    for (let r = seg.start; r <= seg.end; r++) {
      const rr = residues.get(r);
      if (rr) {
        for (const a of rr.atoms) {
          const shifted = { x: a.x - cx, y: a.y - cy, z: a.z - cz };
          const rotated = rotate(shifted, R);
          a.newX = rotated.x + currentX;
          a.newY = rotated.y;
          a.newZ = rotated.z;
        }
      }
    }
    
    tmAnchors.set(i, { 
      startCA: { x: firstRes.ca.newX, y: firstRes.ca.newY, z: firstRes.ca.newZ },
      endCA: { x: lastRes.ca.newX, y: lastRes.ca.newY, z: lastRes.ca.newZ }
    });
  }

  // 4. Process loops
  for (const [resSeq, res] of residues.entries()) {
    // Check if in TM
    let inTM = false;
    for (const seg of tmSegments) {
      if (resSeq >= seg.start && resSeq <= seg.end) {
        inTM = true;
        break;
      }
    }
    if (inTM) continue;

    // Find bounding TMs
    let prevTM = -1, nextTM = -1;
    for (let i = 0; i < tmSegments.length; i++) {
      if (tmSegments[i].end < resSeq) prevTM = i;
      if (tmSegments[i].start > resSeq && nextTM === -1) nextTM = i;
    }

    let targetCA = { x: res.ca.x, y: res.ca.y, z: res.ca.z };

    if (prevTM !== -1 && nextTM !== -1) {
      // Loop between TMs
      const startRes = tmSegments[prevTM].end;
      const endRes = tmSegments[nextTM].start;
      const numRes = endRes - startRes; 
      const f = (resSeq - startRes) / numRes;
      
      const p1 = tmAnchors.get(prevTM)!.endCA;
      const p2 = tmAnchors.get(nextTM)!.startCA;
      
      const dist = Math.sqrt(Math.pow(p2.x-p1.x, 2) + Math.pow(p2.y-p1.y, 2) + Math.pow(p2.z-p1.z, 2));
      const S = numRes * 3.7; // Target path length per bond is ~3.7A
      
      let H = 0;
      if (S > dist) {
        // Calculate bulge height to consume slack
        H = Math.sqrt((4 * dist * (S - dist)) / (Math.PI * Math.PI));
      }
      
      // Bulge in Y direction (Rainbow arch) so it's fully visible in 2D projection
      const bulgeDir = prevTM % 2 === 0 ? 1 : -1;
      const arcY = Math.sin(f * Math.PI) * H * bulgeDir;
      
      // Add a tiny spiral/wiggle in Z to keep the 3D aesthetic
      const arcZ = Math.sin(f * Math.PI * 2) * (H * 0.15); 
      
      targetCA = {
        x: p1.x + f * (p2.x - p1.x),
        y: p1.y + f * (p2.y - p1.y) + arcY,
        z: p1.z + f * (p2.z - p1.z) + arcZ
      };
    } else if (prevTM !== -1) {
      // C-term tail
      const p1 = tmAnchors.get(prevTM)!.endCA;
      const step = resSeq - tmSegments[prevTM].end;
      const bulgeDir = prevTM % 2 === 0 ? 1 : -1;
      // Spiral outward and upwards/downwards
      targetCA = { 
        x: p1.x + step * 2.2, 
        y: p1.y + step * 1.5 * bulgeDir + Math.sin(step * 0.6) * 6, 
        z: p1.z + Math.cos(step * 0.6) * 6 
      };
    } else if (nextTM !== -1) {
      // N-term tail
      const p2 = tmAnchors.get(nextTM)!.startCA;
      const step = tmSegments[nextTM].start - resSeq;
      const startDir = nextTM % 2 === 0 ? -1 : 1; // If TM0 points UP, start is at BOTTOM (-1)
      targetCA = { 
        x: p2.x - step * 2.2, 
        y: p2.y + step * 1.5 * startDir + Math.sin(step * 0.6) * 6, 
        z: p2.z + Math.cos(step * 0.6) * 6 
      };
    }

    // Apply translation to all atoms in this loop residue
    const dx = targetCA.x - res.ca.x;
    const dy = targetCA.y - res.ca.y;
    const dz = targetCA.z - res.ca.z;
    
    for (const a of res.atoms) {
      a.newX = a.x + dx;
      a.newY = a.y + dy;
      a.newZ = a.z + dz;
    }
  }

  // 5. Rewrite PDB string
  for (const a of atoms) {
    const line = lines[a.index];
    const newX = a.newX.toFixed(3).padStart(8, ' ');
    const newY = a.newY.toFixed(3).padStart(8, ' ');
    const newZ = a.newZ.toFixed(3).padStart(8, ' ');
    lines[a.index] = line.substring(0, 30) + newX + newY + newZ + line.substring(54);
  }

  return lines.join('\n');
}
