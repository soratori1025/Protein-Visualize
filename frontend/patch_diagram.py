import sys
import re

file_path = 'src/components/topology/TransmembraneTopologyDiagram.tsx'

with open(file_path, 'r', encoding='utf-8') as f:
    content = f.read()

# 1. Update calcAlgorithm state type
content = content.replace(
    "const [calcAlgorithm, setCalcAlgorithm] = useState<'dssp_slab' | 'kd_slab'>('dssp_slab');",
    "const [calcAlgorithm, setCalcAlgorithm] = useState<'tmhmm_seq' | 'dssp_slab' | 'kd_slab'>('tmhmm_seq');\n  const [overlayType, setOverlayType] = useState<'none' | 'dssp' | 'stride'>('none');"
)

# 2. Add hatched pattern def
hatched_pattern = '''
            <pattern id="hatch-warning" width="8" height="8" patternTransform="rotate(45 0 0)" patternUnits="userSpaceOnUse">
              <line x1="0" y1="0" x2="0" y2="8" stroke="#ff0000" strokeWidth="3" opacity="0.8" />
            </pattern>
'''
content = content.replace('</defs>', hatched_pattern + '\n          </defs>')

# 3. Add the Dropdowns
calc_toolbar = '''          ) : (
            <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
              <select 
                className="tm-input-field" 
                value={calcAlgorithm}
                onChange={(e) => setCalcAlgorithm(e.target.value as any)}
                style={{ padding: '4px 8px' }}
              >
                <option value="tmhmm_seq">TM Predictor (Sequence HMM)</option>
                <option value="dssp_slab">Geometry (DSSP + Slab Fit)</option>
                <option value="kd_slab">Geometry (KD + Slab Fit)</option>
              </select>
              <button
                onClick={() => filename && fetchCalculatedTopology(filename, calcAlgorithm)}
                disabled={loadingCalculated || !filename}
                className="tm-add-btn"
              >
                {loadingCalculated ? 'Computing…' : 'Recalculate'}
              </button>
            </div>
          )}'''
content = re.sub(r'\) : \(\s*<button[\s\S]*?Recalculate\'\}\s*</button>\s*\)\}', calc_toolbar + '\n          )}', content)

# 4. Add the Overlay Dropdown next to "Customize colors"
overlay_dropdown = '''
          <select 
            className="tm-input-field" 
            value={overlayType}
            onChange={(e) => setOverlayType(e.target.value as any)}
            style={{ padding: '4px 8px', width: 'auto' }}
            title="Secondary Structure Overlay"
          >
            <option value="none">No Overlay</option>
            <option value="dssp">Overlay DSSP</option>
            <option value="stride">Overlay STRIDE</option>
          </select>
          <button'''
content = content.replace('<button\n            className={`tm-color-toggle-btn ${colorDrawerOpen ? \'active\' : \'\'}`}', overlay_dropdown + '\n            className={`tm-color-toggle-btn ${colorDrawerOpen ? \'active\' : \'\'}`}')

# 5. Add rendering logic for the disagreement hatching
# We want to render this inside the helix group, over the cylinder.
# So we need to calculate which residues are NOT 'H'/'G'/'I' in the secondaryResult.
hatch_logic = '''
                  {/* Disagreement Overlay */}
                  {overlayType !== 'none' && secondaryResult && (() => {
                     const ssResidues = secondaryResult.residues.filter(r => r.chain_id === chain?.id && r.residue_number >= h.startRes && r.residue_number <= h.endRes);
                     if (ssResidues.length === 0) return null;
                     const mismatched = ssResidues.filter(r => !['H', 'G', 'I'].includes(r.code.toUpperCase()));
                     if (mismatched.length === 0) return null;
                     
                     // Draw bands for mismatched regions
                     const bands = [];
                     let startIdx = -1;
                     for (let i = 0; i < ssResidues.length; i++) {
                       const isMismatch = !['H', 'G', 'I'].includes(ssResidues[i].code.toUpperCase());
                       if (isMismatch && startIdx === -1) startIdx = i;
                       if (!isMismatch && startIdx !== -1) {
                         bands.push({ start: startIdx, end: i - 1 });
                         startIdx = -1;
                       }
                     }
                     if (startIdx !== -1) bands.push({ start: startIdx, end: ssResidues.length - 1 });
                     
                     return bands.map((band, idx) => {
                       const yStart = 6 + (band.start / ssResidues.length) * (cylHeight - 12);
                       const yHeight = Math.max(2, ((band.end - band.start + 1) / ssResidues.length) * (cylHeight - 12));
                       return (
                         <rect
                           key={`mismatch-${idx}`}
                           x="0"
                           y={yStart}
                           width={helixWidth}
                           height={yHeight}
                           fill="url(#hatch-warning)"
                           style={{ pointerEvents: 'none' }}
                         />
                       );
                     });
                  })()}
                  
                  <text'''
content = content.replace('<text\n                    x={cx}\n                    y={cy - 2}', hatch_logic + '\n                  <text\n                    x={cx}\n                    y={cy - 2}')

with open(file_path, 'w', encoding='utf-8') as f:
    f.write(content)
print("Updated TransmembraneTopologyDiagram.tsx successfully")
