with open("frontend/src/components/topology/TransmembraneTopologyDiagram.tsx", "r", encoding="utf-8") as f:
    content = f.read()

# Replace calcAlgorithm state
old_state = """  const [calcAlgorithm, setCalcAlgorithm] = useState<
    'dssp_ss' | 'stride_ss' | 'dssp_slab' | 'stride_slab' | 'kd_slab' | 'tmhmm_seq'
  >('dssp_ss');"""
new_state = """  const [tmAlgorithm, setTmAlgorithm] = useState<string>('3d_slab_geom');
  const [ssAlgorithm, setSsAlgorithm] = useState<string>('dssp');
  const [flowType, setFlowType] = useState<string>('ss_then_tm');
  const [customUniprotId, setCustomUniprotId] = useState<string>('');"""
content = content.replace(old_state, new_state)

# Replace fetchCalculatedTopology signature and query params
old_fetch = """  const fetchCalculatedTopology = useCallback(async (filenameToFetch: string, algorithm: string) => {
    if (!filenameToFetch.trim()) return;
    setLoadingCalculated(true);
    setCalculatedError(null);
    try {
      // Build query string with algorithm + any non-default advanced params
      const qp = new URLSearchParams({ algorithm });"""
new_fetch = """  const fetchCalculatedTopology = useCallback(async (filenameToFetch: string, tm: string, ss: string, flow: string, cUni: string) => {
    if (!filenameToFetch.trim()) return;
    setLoadingCalculated(true);
    setCalculatedError(null);
    try {
      const qp = new URLSearchParams({ tm_algo: tm, ss_algo: ss, flow_type: flow });
      if (tm === 'uniprot_api' && cUni) qp.set('uniprot_id', cUni);"""
content = content.replace(old_fetch, new_fetch)

# Replace useEffect dependencies and calls
content = content.replace("fetchCalculatedTopology(filename, calcAlgorithm)", "fetchCalculatedTopology(filename, tmAlgorithm, ssAlgorithm, flowType, customUniprotId)")
content = content.replace("}, [filename, calcAlgorithm]);", "}, [filename, tmAlgorithm, ssAlgorithm, flowType, customUniprotId]);")

# Replace UI in render
old_ui = """              <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                <select 
                  className="tm-input-field" 
                  value={calcAlgorithm}
                  onChange={(e) => setCalcAlgorithm(e.target.value as any)}
                  style={{ padding: '4px 8px' }}
                >
                  <option value="dssp_ss">DSSP elements + membrane (recommended)</option>
                  <option value="stride_ss">STRIDE elements + membrane</option>
                  <option value="dssp_slab">DSSP + slab (legacy)</option>
                  <option value="stride_slab">STRIDE + slab (legacy)</option>
                  <option value="kd_slab">Geometry only (no SS)</option>
                  <option value="tmhmm_seq">Sequence only (no 3D)</option>
                </select>
                <button
                  onClick={() => filename && fetchCalculatedTopology(filename, calcAlgorithm)}
                  disabled={loadingCalculated || !filename}
                  className="tm-add-btn"
                >
                  {loadingCalculated ? 'Computing…' : 'Recalculate'}
                </button>
                <button
                  onClick={() => setShowAdvancedParams(!showAdvancedParams)}
                  className="tm-add-btn"
                >
                  ⚙️ Advanced
                </button>
              </div>"""

new_ui = """              <div style={{ display: 'flex', flexDirection: 'column', gap: '8px', border: '1px solid #333', padding: '8px', borderRadius: '4px' }}>
                <div style={{ display: 'flex', alignItems: 'center', gap: '8px', flexWrap: 'wrap' }}>
                  <select 
                    className="tm-input-field" 
                    value={tmAlgorithm}
                    onChange={(e) => setTmAlgorithm(e.target.value)}
                    style={{ padding: '4px 8px' }}
                  >
                    <option value="3d_slab_geom">3D Slab Geometry (Recommended)</option>
                    <option value="kyte_doolittle_seq">Kyte-Doolittle Sequence</option>
                    <option value="uniprot_api">UniProt API</option>
                  </select>

                  {tmAlgorithm === 'uniprot_api' && (
                    <input 
                      type="text"
                      className="tm-input-field"
                      placeholder="UniProt ID (e.g. P31645)"
                      value={customUniprotId}
                      onChange={(e) => setCustomUniprotId(e.target.value)}
                      style={{ width: '150px', padding: '4px 8px' }}
                    />
                  )}

                  <select 
                    className="tm-input-field" 
                    value={ssAlgorithm}
                    onChange={(e) => setSsAlgorithm(e.target.value)}
                    style={{ padding: '4px 8px' }}
                  >
                    <option value="dssp">DSSP</option>
                    <option value="stride">STRIDE</option>
                    <option value="none">None (Only TM boundaries)</option>
                  </select>

                  <select 
                    className="tm-input-field" 
                    value={flowType}
                    onChange={(e) => setFlowType(e.target.value)}
                    style={{ padding: '4px 8px' }}
                  >
                    <option value="ss_then_tm">Filter TM by SS (Recommended)</option>
                    <option value="tm_then_ss">Filter SS by TM</option>
                    <option value="parallel_merge">Parallel Merge (Strict Intersection)</option>
                  </select>
                </div>
                
                <div style={{ display: 'flex', alignItems: 'center', gap: '8px' }}>
                  <button
                    onClick={() => filename && fetchCalculatedTopology(filename, tmAlgorithm, ssAlgorithm, flowType, customUniprotId)}
                    disabled={loadingCalculated || !filename || (tmAlgorithm === 'uniprot_api' && !customUniprotId)}
                    className="tm-add-btn"
                  >
                    {loadingCalculated ? 'Computing…' : 'Recalculate'}
                  </button>
                  <button
                    onClick={() => setShowAdvancedParams(!showAdvancedParams)}
                    className="tm-add-btn"
                  >
                    ⚙️ Advanced Parameters
                  </button>
                  {tmAlgorithm === 'uniprot_api' && !customUniprotId && (
                    <span style={{color: '#ff4444', fontSize: '12px'}}>* UniProt ID is required</span>
                  )}
                </div>
              </div>"""

content = content.replace(old_ui, new_ui)

with open("frontend/src/components/topology/TransmembraneTopologyDiagram.tsx", "w", encoding="utf-8") as f:
    f.write(content)
