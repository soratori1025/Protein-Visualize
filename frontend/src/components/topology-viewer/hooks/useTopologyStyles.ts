import { useState } from 'react';
import type { FigureTheme, CustomResidueColorRule } from '../topologyUtils';

export function useTopologyStyles() {
  const [figureTheme, setFigureTheme] = useState<FigureTheme>('publication');
  const [visualStyle, setVisualStyle] = useState<'cylinder' | 'ribbon' | 'wire' | 'flat' | 'beads'>('cylinder');
  const [selectedPaletteKey, setSelectedPaletteKey] = useState<string>('PAPER_DEFAULT');
  const [customHelixColors, setCustomHelixColors] = useState<Record<string, string>>({});
  const [customRegionColors, setCustomRegionColors] = useState<Record<string, string>>({
    Membrane: '#00d2d3',
    ExtracellularText: '#ff9f43',
    CytoplasmicText: '#5f27cd',
  });
  const [customResidueRules, setCustomResidueRules] = useState<CustomResidueColorRule[]>([]);
  const [hoverBrightness, setHoverBrightness] = useState<number>(1.1);
  const [hoverShadow, setHoverShadow] = useState<number>(0.3);

  const handleResetColors = () => {
    setSelectedPaletteKey('PAPER_DEFAULT');
    setCustomHelixColors({});
    setCustomResidueRules([]);
  };

  const getCustomResidueColorForHelix = (startRes: number, endRes: number): string | null => {
    const matchingRule = customResidueRules.find(
      (rule) => rule.startRes <= endRes && rule.endRes >= startRes
    );
    return matchingRule ? matchingRule.color : null;
  };

  const handleAddResidueRule = (start: number, end: number, color: string) => {
    if (isNaN(start)) return;
    setCustomResidueRules((prev) => [
      ...prev,
      {
        id: `rule-${Date.now()}`,
        label: start === end ? `Residue ${start}` : `Residues ${start}–${end}`,
        startRes: Math.min(start, end),
        endRes: Math.max(start, end),
        color,
      },
    ]);
  };

  const handleRemoveResidueRule = (id: string) => {
    setCustomResidueRules((prev) => prev.filter((r) => r.id !== id));
  };

  return {
    figureTheme,
    setFigureTheme,
    visualStyle,
    setVisualStyle,
    selectedPaletteKey,
    setSelectedPaletteKey,
    customHelixColors,
    setCustomHelixColors,
    customRegionColors,
    setCustomRegionColors,
    customResidueRules,
    setCustomResidueRules,
    hoverBrightness,
    setHoverBrightness,
    hoverShadow,
    setHoverShadow,
    handleResetColors,
    getCustomResidueColorForHelix,
    handleAddResidueRule,
    handleRemoveResidueRule,
  };
}
