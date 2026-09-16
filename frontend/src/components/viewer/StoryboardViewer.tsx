import { useEffect, useRef, useState } from 'react';
import * as $3Dmol from '3dmol';

export interface KeyframeConfig {
  id: string;
  title: string;
  description: string;
  focus?: { chain?: string; resi?: number[] };
  style: 'ribbon' | 'surface' | 'stick' | 'cartoon';
  colorScheme: 'chain' | 'ss' | 'hydropathy' | 'bFactor';
  highlights: { chain: string; resi: number[] }[];
}

interface Props {
  filename: string | undefined;
  activeKeyframe: KeyframeConfig;
}

export function StoryboardViewer({ filename, activeKeyframe }: Props) {
  const containerRef = useRef<HTMLDivElement>(null);
  const viewerRef = useRef<$3Dmol.Viewer>();
  const [loaded, setLoaded] = useState(false);

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
      
      setLoaded(true);
    };

    void loadViewer();
    return () => {
      cancelled = true;
      viewer?.clear();
      viewerRef.current = undefined;
      if (containerRef.current) containerRef.current.replaceChildren();
    };
  }, [filename]);

  // Apply Keyframe Logic
  useEffect(() => {
    const viewer = viewerRef.current;
    if (!viewer || !loaded) return;

    // 1. Reset all styles
    (viewer as any).removeAllShapes();
    viewer.setStyle({}, { cartoon: { hidden: true }, stick: { hidden: true }, sphere: { hidden: true }, surface: { hidden: true } });
    (viewer as any).removeAllSurfaces();

    // 2. Apply base style
    if (activeKeyframe.style === 'ribbon') {
      viewer.setStyle({}, { cartoon: { color: 'spectrum', opacity: 1 } });
    } else if (activeKeyframe.style === 'surface') {
      viewer.setStyle({}, { cartoon: { color: 'spectrum', opacity: 0.3 } });
      (viewer as any).addSurface(($3Dmol as any).SurfaceType?.VDW || 1, { opacity: 0.8, color: 'white' }, { hetflag: false }, { hetflag: false });
    } else if (activeKeyframe.style === 'stick') {
      viewer.setStyle({}, { stick: { colorscheme: 'Jmol', radius: 0.15 } });
    } else if (activeKeyframe.style === 'cartoon') {
      viewer.setStyle({}, { cartoon: { colorscheme: 'ssPyMOL' } });
    }

    // 3. Apply Highlights
    if (activeKeyframe.highlights.length > 0) {
      // Make everything else transparent/dim
      if (activeKeyframe.style === 'ribbon') {
        viewer.setStyle({}, { cartoon: { color: '#2a4156', opacity: 0.5 } });
      }

      let validHighlights = false;
      activeKeyframe.highlights.forEach(hl => {
        hl.resi.forEach(res => {
          const sel = { chain: hl.chain, resi: res };
          if ((viewer as any).selectedAtoms(sel).length > 0) {
            validHighlights = true;
          }
          viewer.setStyle(sel, {
            cartoon: { color: '#ff6f61', opacity: 1 },
            stick: { colorscheme: 'Jmol', radius: 0.3 }
          });
        });
      });

      // Fallback: If the uploaded PDB doesn't have the hardcoded residues, pick some dynamically!
      if (!validHighlights) {
        const allAtoms = (viewer as any).getModel().selectedAtoms({});
        if (allAtoms.length > 0) {
          const firstChain = allAtoms[0].chain;
          const chainAtoms = allAtoms.filter((a: any) => a.chain === firstChain);
          const resis = Array.from(new Set(chainAtoms.map((a: any) => a.resi))).filter(r => r !== undefined);
          
          if (resis.length > 2) {
            const middle = Math.floor(resis.length / 2);
            const r1 = resis[middle];
            const r2 = resis[middle + 1];
            
            [{ chain: firstChain, resi: r1 }, { chain: firstChain, resi: r2 }].forEach(sel => {
              viewer.setStyle(sel, {
                cartoon: { color: '#ff6f61', opacity: 1 },
                stick: { colorscheme: 'Jmol', radius: 0.3 }
              });
            });
            
            // Override focus to zoom into these dynamic residues
            activeKeyframe.focus = { chain: firstChain, resi: [Number(r1), Number(r2)] };
          }
        }
      }
    }

    // 4. Animate Camera Transition
    if (activeKeyframe.focus) {
      const atoms = (viewer as any).selectedAtoms(activeKeyframe.focus);
      if (atoms.length > 0) {
        (viewer as any).zoomTo(activeKeyframe.focus, 1000); // 1000ms animation
      } else {
        viewer.zoomTo(); // Fallback if selection doesn't exist
      }
    } else {
      viewer.zoomTo();
    }

    viewer.render();
  }, [activeKeyframe, loaded]);

  if (!filename) return <div className="empty-state">No structure loaded for Storyboard.</div>;

  return (
    <div ref={containerRef} style={{ width: '100%', height: '100%' }} />
  );
}
