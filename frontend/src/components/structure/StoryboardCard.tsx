import React from 'react';
import type { KeyframeConfig } from '../viewer/StoryboardViewer';

interface StoryboardCardProps {
  kf: KeyframeConfig;
  isActive: boolean;
  onSelect: (id: string) => void;
}

export function StoryboardCard({ kf, isActive, onSelect }: StoryboardCardProps) {
  return (
    <div
      data-id={kf.id}
      className={`keyframe-card ${isActive ? 'active' : ''}`}
      onClick={() => {
        onSelect(kf.id);
        // Scroll the clicked card into the center of the view smoothly
        document
          .querySelector(`[data-id="${kf.id}"]`)
          ?.scrollIntoView({ behavior: 'smooth', block: 'center' });
      }}
      style={{
        minHeight: '50vh',
        display: 'flex',
        flexDirection: 'column',
        justifyContent: 'center',
        marginBottom: '30vh',
      }}
    >
      <h3>{kf.title}</h3>
      <p>{kf.description}</p>
    </div>
  );
}
