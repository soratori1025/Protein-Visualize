import sys
import re

with open("src/components/topology/TransmembraneTopologyDiagram.tsx", "r", encoding="utf-8") as f:
    content = f.read()

# We need to insert a block before {/* Break connectors between the two halves of a discontinuous helix */}

terminals_code = """
          {/* N-terminus and C-terminus tails */}
          {(() => {
            if (helices.length === 0) return null;
            const first = helices[0];
            const last = helices[helices.length - 1];
            const posFirst = helixPositions[first.id];
            const posLast = helixPositions[last.id];
            if (!posFirst || !posLast) return null;

            const nFeatures = getExtraFeatures(1, first.startRes - 1, usingAnnotation ? activeTopologyData : null, secondaryResult, chain?.id) || [];
            const lastResNum = chain?.residues[chain.residues.length - 1]?.residue_number || 10000;
            const cFeatures = getExtraFeatures(last.endRes + 1, lastResNum, usingAnnotation ? activeTopologyData : null, secondaryResult, chain?.id) || [];

            return (
              <g className="tm-terminals">
                {/* N-terminus */}
                <path
                  d={`M 20 ${first.entrySide === 'out' ? membraneTopY - 60 : membraneBottomY + 60} Q ${posFirst.x / 2} ${first.entrySide === 'out' ? membraneTopY - 30 : membraneBottomY + 30} ${posFirst.x + helixWidth / 2} ${posFirst.nEndY}`}
                  fill="none"
                  stroke={first.color}
                  strokeWidth="2.6"
                />
                <text x="20" y={first.entrySide === 'out' ? membraneTopY - 70 : membraneBottomY + 75} textAnchor="middle" style={{ fill: isPub ? '#1e293b' : '#e2e8f0', fontSize: '11px', fontWeight: 'bold' }}>NH2</text>
                
                {nFeatures.map((f, i) => {
                   const y = first.entrySide === 'out' ? membraneTopY - 45 : membraneBottomY + 45;
                   const x = 30 + i * 50;
                   return (
                     <g key={`n-${i}`} transform={`translate(${x}, ${y})`} onMouseEnter={(e) => {
                       e.stopPropagation();
                       setHoveredElement({ title: `${f.type}: ${f.label}`, range: `Residues ${f.startRes}-${f.endRes}`, length: f.endRes - f.startRes + 1, details: 'N-terminus' });
                     }} onMouseLeave={() => setHoveredElement(null)}>
                       <rect x="0" y="0" width="45" height="18" rx={f.type === 'Helix' ? 9 : 4} fill="url(#short-helix-grad)" stroke={isPub ? '#475569' : '#0f172a'} />
                       <text x="22.5" y="13" textAnchor="middle" style={{ fill: '#0f172a', fontSize: '9px' }}>{f.label.substring(0, 7)}</text>
                     </g>
                   );
                })}

                {/* C-terminus */}
                <path
                  d={`M ${posLast.x + helixWidth / 2} ${posLast.cEndY} Q ${(posLast.x + canvasWidth) / 2} ${last.exitSide === 'out' ? membraneTopY - 30 : membraneBottomY + 30} ${canvasWidth - 30} ${last.exitSide === 'out' ? membraneTopY - 60 : membraneBottomY + 60}`}
                  fill="none"
                  stroke={last.color}
                  strokeWidth="2.6"
                />
                <text x={canvasWidth - 30} y={last.exitSide === 'out' ? membraneTopY - 70 : membraneBottomY + 75} textAnchor="middle" style={{ fill: isPub ? '#1e293b' : '#e2e8f0', fontSize: '11px', fontWeight: 'bold' }}>COOH</text>
                
                {cFeatures.map((f, i) => {
                   const y = last.exitSide === 'out' ? membraneTopY - 45 : membraneBottomY + 45;
                   const x = canvasWidth - 80 - i * 50;
                   return (
                     <g key={`c-${i}`} transform={`translate(${x}, ${y})`} onMouseEnter={(e) => {
                       e.stopPropagation();
                       setHoveredElement({ title: `${f.type}: ${f.label}`, range: `Residues ${f.startRes}-${f.endRes}`, length: f.endRes - f.startRes + 1, details: 'C-terminus' });
                     }} onMouseLeave={() => setHoveredElement(null)}>
                       <rect x="0" y="0" width="45" height="18" rx={f.type === 'Helix' ? 9 : 4} fill="url(#short-helix-grad)" stroke={isPub ? '#475569' : '#0f172a'} />
                       <text x="22.5" y="13" textAnchor="middle" style={{ fill: '#0f172a', fontSize: '9px' }}>{f.label.substring(0, 7)}</text>
                     </g>
                   );
                })}
              </g>
            );
          })()}

"""

content = content.replace("{/* Break connectors between the two halves of a discontinuous helix */}", terminals_code + "          {/* Break connectors between the two halves of a discontinuous helix */}")

with open("src/components/topology/TransmembraneTopologyDiagram.tsx", "w", encoding="utf-8") as f:
    f.write(content)
print("Updated successfully")
