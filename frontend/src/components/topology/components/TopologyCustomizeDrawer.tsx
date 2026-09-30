import React, { useState } from 'react';
import { buildPalette, PALETTES, TMHelix } from '../topologyUtils';

interface Props {
  isPub: boolean;
  activeTab: string;
  setActiveTab: (tab: 'preset' | 'helices' | 'residues' | 'effects' | 'regions') => void;
  visualStyle: string;
  setVisualStyle: (v: any) => void;
  selectedPaletteKey: string;
  setSelectedPaletteKey: (v: string) => void;
  helices: TMHelix[];
  customHelixColors: Record<string, string>;
  setCustomHelixColors: (v: (prev: any) => any) => void;
  customRegionColors: Record<string, string>;
  setCustomRegionColors: (v: (prev: any) => any) => void;
  customResidueRules: any[];
  resStartInput: string;
  setResStartInput: (v: string) => void;
  resEndInput: string;
  setResEndInput: (v: string) => void;
  resColorInput: string;
  setResColorInput: (v: string) => void;
  handleAddResidueRule: () => void;
  handleRemoveResidueRule: (id: string) => void;
  handleResetColors: () => void;
  columnCount: number;
}

export function TopologyCustomizeDrawer({
  isPub, activeTab, setActiveTab, visualStyle, setVisualStyle, selectedPaletteKey, setSelectedPaletteKey,
  helices, customHelixColors, setCustomHelixColors, customRegionColors, setCustomRegionColors,
  customResidueRules, resStartInput, setResStartInput, resEndInput, setResEndInput, resColorInput, setResColorInput,
  handleAddResidueRule, handleRemoveResidueRule, handleResetColors, columnCount
}: Props) {
  return (
    <div className={`tm-color-customizer-drawer ${isPub ? 'publication' : 'lab'}`}>
      <div className="tm-drawer-tabs">
        <button className={`tm-tab-btn ${activeTab === 'preset' ? 'active' : ''}`} onClick={() => setActiveTab('preset')}>Palettes</button>
        <button className={`tm-tab-btn ${activeTab === 'helices' ? 'active' : ''}`} onClick={() => setActiveTab('helices')}>Individual helices</button>
        <button className={`tm-tab-btn ${activeTab === 'residues' ? 'active' : ''}`} onClick={() => setActiveTab('residues')}>Residue ranges</button>
        <button className={`tm-tab-btn ${activeTab === 'effects' ? 'active' : ''}`} onClick={() => setActiveTab('effects')}>Styles display</button>
        <button className={`tm-tab-btn ${activeTab === 'regions' ? 'active' : ''}`} onClick={() => setActiveTab('regions')}>Regions</button>
        <button className="tm-reset-btn" onClick={handleResetColors}>Reset setting</button>
      </div>

      {activeTab === 'effects' && (
        <div className="tm-effects-section" style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
            <label style={{ fontSize: '11px', color: isPub ? '#334155' : '#94a3b8', fontWeight: 600 }}>TM Helix Style</label>
            <div style={{ display: 'flex', flexWrap: 'wrap', gap: '12px', marginTop: '4px' }}>
              {[{ id: 'cylinder', label: 'Cylinder' }, { id: 'ribbon', label: 'Ribbon' }, { id: 'flat', label: 'Flat Block' }].map(style => (
                <label key={style.id} style={{ display: 'flex', alignItems: 'center', gap: '6px', fontSize: '12px', color: isPub ? '#334155' : '#cbd5e1', cursor: 'pointer' }}>
                  <input type="radio" name="visualStyle" checked={visualStyle === style.id} onChange={() => setVisualStyle(style.id)} /> {style.label}
                </label>
              ))}
            </div>
          </div>
        </div>
      )}

      {activeTab === 'regions' && (
        <div className="tm-effects-section" style={{ display: 'flex', flexDirection: 'column', gap: '12px' }}>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
            <label style={{ fontSize: '11px', color: isPub ? '#334155' : '#94a3b8', fontWeight: 600 }}>Membrane Background</label>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <input type="color" value={customRegionColors.Membrane || '#00d2d3'} onChange={(e) => setCustomRegionColors(prev => ({ ...prev, Membrane: e.target.value }))} style={{ width: '24px', height: '24px', padding: '0', border: 'none', cursor: 'pointer', background: 'transparent' }} />
              <span style={{ fontSize: '11px', color: isPub ? '#64748b' : '#94a3b8' }}>{customRegionColors.Membrane || '#00d2d3'}</span>
            </div>
          </div>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
            <label style={{ fontSize: '11px', color: isPub ? '#334155' : '#94a3b8', fontWeight: 600 }}>Extracellular Label Color</label>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <input type="color" value={customRegionColors.ExtracellularText || '#ff9f43'} onChange={(e) => setCustomRegionColors(prev => ({ ...prev, ExtracellularText: e.target.value }))} style={{ width: '24px', height: '24px', padding: '0', border: 'none', cursor: 'pointer', background: 'transparent' }} />
              <span style={{ fontSize: '11px', color: isPub ? '#64748b' : '#94a3b8' }}>{customRegionColors.ExtracellularText || '#ff9f43'}</span>
            </div>
          </div>
          <div style={{ display: 'flex', flexDirection: 'column', gap: '4px' }}>
            <label style={{ fontSize: '11px', color: isPub ? '#334155' : '#94a3b8', fontWeight: 600 }}>Cytoplasmic Label Color</label>
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <input type="color" value={customRegionColors.CytoplasmicText || '#5f27cd'} onChange={(e) => setCustomRegionColors(prev => ({ ...prev, CytoplasmicText: e.target.value }))} style={{ width: '24px', height: '24px', padding: '0', border: 'none', cursor: 'pointer', background: 'transparent' }} />
              <span style={{ fontSize: '11px', color: isPub ? '#64748b' : '#94a3b8' }}>{customRegionColors.CytoplasmicText || '#5f27cd'}</span>
            </div>
          </div>
        </div>
      )}

      {activeTab === 'preset' && (
        <div className="tm-palette-grid">
          {Object.entries(PALETTES).map(([key, palette]) => (
            <button key={key} className={`tm-palette-btn ${selectedPaletteKey === key ? 'active' : ''}`} onClick={() => { setSelectedPaletteKey(key); setCustomHelixColors(() => ({})); }}>
              <div className="tm-palette-preview">
                {buildPalette(palette.colors, Math.max(6, Math.min(12, columnCount))).map((c, i) => (
                  <span key={i} style={{ backgroundColor: c }} />
                ))}
              </div>
              <span>{palette.label}</span>
            </button>
          ))}
        </div>
      )}

      {activeTab === 'helices' && (
        <div className="tm-helix-color-grid">
          {helices.map((h) => (
            <div key={`color-${h.id}`} className="tm-helix-color-item">
              <span>TM{h.subLabel}</span>
              <input type="color" className="tm-color-input" value={h.color} onChange={(e) => setCustomHelixColors((prev) => ({ ...prev, [h.subLabel]: e.target.value }))} />
            </div>
          ))}
        </div>
      )}

      {activeTab === 'residues' && (
        <div className="tm-residue-color-section">
          <div className="tm-residue-form">
            <input type="number" className="tm-input-field" placeholder="Start" value={resStartInput} onChange={(e) => setResStartInput(e.target.value)} />
            <span>–</span>
            <input type="number" className="tm-input-field" placeholder="End (optional)" value={resEndInput} onChange={(e) => setResEndInput(e.target.value)} />
            <input type="color" className="tm-color-input" value={resColorInput} onChange={(e) => setResColorInput(e.target.value)} />
            <button className="tm-add-btn" onClick={handleAddResidueRule}>Add highlight</button>
          </div>
          {customResidueRules.length > 0 ? (
            <div className="tm-residue-tag-list">
              {customResidueRules.map((rule) => (
                <span key={rule.id} className="tm-residue-tag" style={{ borderLeftColor: rule.color }}>
                  <span className="tm-rule-color-dot" style={{ backgroundColor: rule.color }} />
                  {rule.label}
                  <button onClick={() => handleRemoveResidueRule(rule.id)}>×</button>
                </span>
              ))}
            </div>
          ) : (
            <small style={{ color: '#64748b', fontSize: '11px' }}>Enter a residue number or range to highlight positions on the map.</small>
          )}
        </div>
      )}
    </div>
  );
}
