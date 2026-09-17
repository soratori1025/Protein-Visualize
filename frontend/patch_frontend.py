import sys
import re

with open("src/components/topology/TransmembraneTopologyDiagram.tsx", "r", encoding="utf-8") as f:
    content = f.read()

# 1. Add state for calcAlgorithm
state_decl = "const [topologySource, setTopologySource] = useState<TopologySource>('uniprot');\n  const [calcAlgorithm, setCalcAlgorithm] = useState<'dssp_slab' | 'kd_slab'>('dssp_slab');"
content = content.replace("const [topologySource, setTopologySource] = useState<TopologySource>('uniprot');", state_decl)

# 2. Update fetchCalculatedTopology
old_fetch = r"""const fetchCalculatedTopology = useCallback(async (filenameToFetch: string) => {
    if (!filenameToFetch.trim()) return;
    setLoadingCalculated(true);
    setCalculatedError(null);
    try {
      const response = await fetch(
        `http://localhost:8000/api/secondary-structure/predict-topology/${encodeURIComponent(filenameToFetch.trim())}`
      );"""

new_fetch = r"""const fetchCalculatedTopology = useCallback(async (filenameToFetch: string, algorithm: string) => {
    if (!filenameToFetch.trim()) return;
    setLoadingCalculated(true);
    setCalculatedError(null);
    try {
      const response = await fetch(
        `http://localhost:8000/api/secondary-structure/predict-topology/${encodeURIComponent(filenameToFetch.trim())}?algorithm=${algorithm}`
      );"""

content = content.replace(old_fetch, new_fetch)

# 3. Fix fetchCalculatedTopology calls in useEffect
content = content.replace("if (topologySource === 'calculated' && filename) fetchCalculatedTopology(filename);", "if (topologySource === 'calculated' && filename) fetchCalculatedTopology(filename, calcAlgorithm);")
content = content.replace("fetchCalculatedTopology(filename);", "fetchCalculatedTopology(filename, calcAlgorithm);")

# 4. Add dependency calcAlgorithm to useEffects
# We need to find the specific useEffects. They are around lines:
# useEffect(() => { ... fetchCalculatedTopology(filename) ... }, [filename])
# Let's replace the whole blocks if possible. Let's just find and replace.

old_effect_1 = """  useEffect(() => {
    setCalculatedData(null);
    setCalculatedError(null);
    if (topologySource === 'calculated' && filename) fetchCalculatedTopology(filename, calcAlgorithm);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [filename]);"""
new_effect_1 = """  useEffect(() => {
    setCalculatedData(null);
    setCalculatedError(null);
    if (topologySource === 'calculated' && filename) fetchCalculatedTopology(filename, calcAlgorithm);
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [filename, calcAlgorithm]);"""
content = content.replace(old_effect_1, new_effect_1)

# 5. Add UI control for algorithm
ui_inject_point = """            <div className="tm-preset-chip" style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
              <button
                className={`tm-tab-btn ${topologySource === 'uniprot' ? 'active' : ''}`}
                onClick={() => setTopologySource('uniprot')}
              >
                UniProt
              </button>
              <button
                className={`tm-tab-btn ${topologySource === 'calculated' ? 'active' : ''}`}
                onClick={() => setTopologySource('calculated')}
                disabled={!filename}
                title={!filename ? 'Upload a structure file to calculate topology' : undefined}
              >
                Calculated (beta)
              </button>
            </div>"""

new_ui = """            <div className="tm-preset-chip" style={{ display: 'flex', gap: '4px', alignItems: 'center' }}>
              <button
                className={`tm-tab-btn ${topologySource === 'uniprot' ? 'active' : ''}`}
                onClick={() => setTopologySource('uniprot')}
              >
                UniProt
              </button>
              <button
                className={`tm-tab-btn ${topologySource === 'calculated' ? 'active' : ''}`}
                onClick={() => setTopologySource('calculated')}
                disabled={!filename}
                title={!filename ? 'Upload a structure file to calculate topology' : undefined}
              >
                Calculated (beta)
              </button>
              {topologySource === 'calculated' && (
                <div style={{ display: 'flex', alignItems: 'center', marginLeft: '12px', gap: '6px' }}>
                  <span style={{ fontSize: '12px', color: '#64748b' }}>Algorithm:</span>
                  <select 
                    value={calcAlgorithm} 
                    onChange={(e) => setCalcAlgorithm(e.target.value as 'dssp_slab' | 'kd_slab')}
                    style={{ fontSize: '12px', padding: '2px 6px', borderRadius: '4px', border: '1px solid #cbd5e1' }}
                  >
                    <option value="dssp_slab">DSSP Slab (Accurate)</option>
                    <option value="kd_slab">KD Slab (Heuristic)</option>
                  </select>
                </div>
              )}
            </div>"""

content = content.replace(ui_inject_point, new_ui)

with open("src/components/topology/TransmembraneTopologyDiagram.tsx", "w", encoding="utf-8") as f:
    f.write(content)
print("Frontend updated successfully!")
