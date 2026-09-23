import { Outlet, NavLink } from 'react-router-dom';
import { ProteinProvider } from '../../contexts/ProteinContext';

export function AppLayout() {
  return (
    <ProteinProvider>
      <div className="app-container">
        <nav className="global-nav">
          <div className="global-nav-brand">🧬 ProteinJournal Framework</div>
          <div className="global-nav-links">
            <NavLink to="/" className={({ isActive }) => (isActive ? 'global-nav-link active' : 'global-nav-link')}>Transmembrane Analysis</NavLink>
            <NavLink to="/visualize" className={({ isActive }) => (isActive ? 'global-nav-link active' : 'global-nav-link')}>Structure Viewer</NavLink>
            <NavLink to="/storyboard" className={({ isActive }) => (isActive ? 'global-nav-link active' : 'global-nav-link')}>Interactive Storyboard</NavLink>
          </div>
        </nav>
        <Outlet />
      </div>
    </ProteinProvider>
  );
}
