import { useState } from 'react';

function App() {
  const [health, setHealth] = useState<string>('Checking...');

  const checkHealth = async () => {
    try {
      const response = await fetch('http://127.0.0.1:8000/api/health');
      const data = await response.json();
      setHealth(JSON.stringify(data));
    } catch (error) {
      setHealth('Backend not running yet');
    }
  };

  return (
    <main className="app-shell">
      <header className="topbar">
        <div>
          <h1>ProteinLab</h1>
          <p>MVP protein visualization platform</p>
        </div>
        <button onClick={checkHealth}>Check API</button>
      </header>

      <section className="panel-grid">
        <div className="panel viewer-panel">
          <h2>3D Viewer</h2>
          <div className="viewer-placeholder">Mol* viewer placeholder</div>
        </div>

        <div className="panel">
          <h2>Sequence</h2>
          <div className="sequence-box">MKTAYIAK... </div>
        </div>
      </section>

      <section className="panel">
        <h2>Topology</h2>
        <div className="topology-box">Domain / motif / contact graph placeholder</div>
      </section>

      <section className="panel status-panel">
        <h2>System status</h2>
        <pre>{health}</pre>
      </section>
    </main>
  );
}

export default App;
