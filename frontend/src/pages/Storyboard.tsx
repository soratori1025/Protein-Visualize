import { useState, useEffect, useRef } from 'react';
import { StoryboardViewer, KeyframeConfig } from '../components/viewer/StoryboardViewer';
import { useProtein } from '../contexts/ProteinContext';
import { Header } from '../components/layout/Header';

const DEMO_KEYFRAMES: KeyframeConfig[] = [
  {
    id: 'intro',
    title: '1. Overall Structure',
    description: "This is a comprehensive view of the molecule's 3D spatial structure. The Ribbon representation clearly shows major folds, including alpha helices and beta sheets, allowing you to see the overall shape of the molecule.",
    style: 'ribbon',
    colorScheme: 'chain',
    highlights: []
  },
  {
    id: 'pocket',
    title: '2. Core and Binding Pockets',
    description: 'Switching to the Surface representation reveals the spatial mesh of the molecular surface. This mode is typically used to investigate grooves and binding pockets - where substrates or ligands attach.',
    style: 'surface',
    colorScheme: 'chain',
    highlights: []
  },
  {
    id: 'interaction',
    title: '3. Crucial Residues',
    description: 'Returning to the Ribbon representation, we highlight important residues. In this example, we zoom into specific residues on the first chain.',
    focus: { resi: [10, 45] },
    style: 'ribbon',
    colorScheme: 'chain',
    highlights: [
      { chain: 'A', resi: [10, 45] }
    ]
  }
];

export function Storyboard() {
  const { protein } = useProtein();
  const [activeId, setActiveId] = useState<string>(DEMO_KEYFRAMES[0].id);
  const activeKeyframe = DEMO_KEYFRAMES.find(k => k.id === activeId) || DEMO_KEYFRAMES[0];
  const containerRef = useRef<HTMLDivElement>(null);

  // Reset active ID when protein changes
  useEffect(() => {
    if (protein?.filename) {
      setActiveId(DEMO_KEYFRAMES[0].id);
    }
  }, [protein?.filename]);

  // IntersectionObserver for "Scroll-telling"
  useEffect(() => {
    if (!containerRef.current) return;
    const observer = new IntersectionObserver(
      (entries) => {
        entries.forEach((entry) => {
          if (entry.isIntersecting) {
            const id = entry.target.getAttribute('data-id');
            if (id) setActiveId(id);
          }
        });
      },
      {
        root: containerRef.current, // Use the scrollable container
        rootMargin: '-40% 0px -40% 0px', // Trigger when card reaches the middle of the container
        threshold: 0
      }
    );
    const elements = containerRef.current.querySelectorAll('.keyframe-card');
    elements.forEach((el) => observer.observe(el));
    return () => observer.disconnect();
  }, [protein?.filename]); // Re-run if filename changes just in case DOM changes, but it's mostly static.

  return (
    <main className="app-shell" style={{ overflow: 'hidden' }}>
      <Header 
        title="Interactive Molecular Storyboarding" 
        subtitle="Scroll through the content blocks below. The 3D model alongside will automatically adjust the camera, representation, and highlight molecules corresponding to what you are reading." 
      />
      <div className="storyboard-layout" style={{ height: 'calc(100vh - 120px)' }}>
        <div className="storyboard-content" ref={containerRef}>


        {DEMO_KEYFRAMES.map((kf) => (
          <div 
            key={kf.id}
            data-id={kf.id}
            className={`keyframe-card ${activeId === kf.id ? 'active' : ''}`}
            onClick={() => {
              setActiveId(kf.id);
              // Scroll the clicked card into the center of the view smoothly
              document.querySelector(`[data-id="${kf.id}"]`)?.scrollIntoView({ behavior: 'smooth', block: 'center' });
            }}
            style={{ minHeight: '50vh', display: 'flex', flexDirection: 'column', justifyContent: 'center', marginBottom: '30vh' }}
          >
            <h3>{kf.title}</h3>
            <p>{kf.description}</p>
          </div>
        ))}
        
        {/* Extra padding at the bottom so the last card can reach the center of the screen */}
        <div style={{ height: '30vh' }}></div>
      </div>
      
      <div className="storyboard-viewer-container">
        {!protein?.filename ? (
           <div className="empty-state" style={{ height: '100%', display: 'flex', alignItems: 'center', justifyContent: 'center', color: '#a5b9cb' }}>
             Please upload a structural file (PDB/mmCIF) to start the story.
           </div>
        ) : (
          <StoryboardViewer filename={protein.filename} activeKeyframe={activeKeyframe} />
        )}
      </div>
    </div>
    </main>
  );
}
