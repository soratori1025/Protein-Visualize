import { useState, useEffect, useRef } from 'react';
import { StoryboardViewer } from '../components/viewer/StoryboardViewer';
import { useProtein } from '../contexts/ProteinContext';
import { Header } from '../components/layout/Header';
import { useScrollSpy } from '../hooks/useScrollSpy';
import { DEMO_KEYFRAMES } from '../data/storyboardDemo';
import { StoryboardCard } from '../components/structure/StoryboardCard';

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
  useScrollSpy({
    containerRef,
    selector: '.keyframe-card',
    onIntersect: setActiveId,
    dependencies: [protein?.filename]
  });

  return (
    <main className="app-shell" style={{ overflow: 'hidden' }}>
      <Header 
        title="Interactive Molecular Storyboarding" 
        subtitle="Scroll through the content blocks below. The 3D model alongside will automatically adjust the camera, representation, and highlight molecules corresponding to what you are reading." 
      />
      <div className="storyboard-layout" style={{ height: 'calc(100vh - 120px)' }}>
        <div className="storyboard-content" ref={containerRef}>

        {DEMO_KEYFRAMES.map((kf) => (
          <StoryboardCard
            key={kf.id}
            kf={kf}
            isActive={activeId === kf.id}
            onSelect={setActiveId}
          />
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
