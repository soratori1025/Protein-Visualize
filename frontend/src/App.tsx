import { BrowserRouter, Routes, Route } from 'react-router-dom';
import { AppLayout } from './components/layout/AppLayout';
import { LabWorkspace } from './pages/LabWorkspace';
import { Storyboard } from './pages/Storyboard';

function App() {
  return (
    <BrowserRouter>
      <Routes>
        <Route path="/" element={<AppLayout />}>
          <Route index element={<LabWorkspace />} />
          <Route path="storyboard" element={<Storyboard />} />
        </Route>
      </Routes>
    </BrowserRouter>
  );
}

export default App;
