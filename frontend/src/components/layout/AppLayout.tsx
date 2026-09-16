import { Outlet, NavLink } from 'react-router-dom';

export function AppLayout() {
  return (
    <div className="app-container">
      <nav className="global-nav">
        <div className="global-nav-brand">🧬 ProteinJournal Framework</div>
        <div className="global-nav-links">
          <NavLink to="/" className={({ isActive }) => (isActive ? 'global-nav-link active' : 'global-nav-link')}>Lab Workspace</NavLink>
          <NavLink to="/storyboard" className={({ isActive }) => (isActive ? 'global-nav-link active' : 'global-nav-link')}>Interactive Storyboard</NavLink>
        </div>
      </nav>
      <Outlet />
    </div>
  );
}
