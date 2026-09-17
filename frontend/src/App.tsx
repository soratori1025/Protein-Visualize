import { BrowserRouter, Routes, Route } from 'react-router-dom';
import { AppLayout } from './components/layout/AppLayout';
import { TransmembraneAnalysis } from './pages/TransmembraneAnalysis';
import { Storyboard } from './pages/Storyboard';
import { StructureViewer } from './pages/StructureViewer';

function App() {
  return (
    <BrowserRouter>
      <Routes>
        <Route path="/" element={<AppLayout />}>
          <Route index element={<TransmembraneAnalysis />} />
          <Route path="visualize" element={<StructureViewer />} />
          <Route path="storyboard" element={<Storyboard />} />
        </Route>
      </Routes>
    </BrowserRouter>
  );
}

export default App;
