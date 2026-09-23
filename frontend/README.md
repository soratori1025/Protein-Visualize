# Protein-Visualize Frontend

This is the frontend application for the Protein-Visualize project. It is a React-based web application built with Vite and TypeScript, designed to visualize protein structures, analyze secondary structures, and create interactive molecular storyboards.

## Tech Stack

- **Framework:** React 18
- **Build Tool:** Vite
- **Language:** TypeScript
- **Routing:** React Router v7
- **Data Visualization:** D3.js (for 2D topology diagrams)
- **Molecular Viewer:** 3Dmol.js (for rendering 3D protein structures)
- **API Client:** Axios

## Features and Views

The application provides a unified workspace with several distinct interactive display types (views):

1. **Transmembrane Analysis (`/`)**
   - Upload PDB / mmCIF structural files.
   - Run secondary structure annotations (DSSP / STRIDE) via the backend.
   - Visualize results in an interactive 2D topology diagram linked to 3D structures.

2. **Structure Viewer (`/visualize`)**
   - Inspect loaded coordinates and sequences.
   - View a "Spread protein map" displaying different chains.
   - Includes a focused 3D workspace to isolate and interact with specific molecular chains.

3. **Interactive Storyboard (`/storyboard`)**
   - A scroll-telling interface designed for presentations.
   - Scrolling through text cards automatically drives the 3D molecular viewer (changing representations, camera focus, and highlighting residues).

## Architecture

The frontend uses a unified Application Shell pattern to make extending the app with new "display types" simple and efficient:

- **`src/contexts/ProteinContext.tsx`**: A global state context that holds the currently loaded molecule data. If a user uploads a protein in one view, that protein remains instantly available when navigating to other views.
- **`src/components/layout/AppLayout.tsx`**: The main shell that wraps the router and provides the `ProteinContext`.
- **`src/components/layout/Header.tsx`**: A shared header component handling file uploads and displaying backend API health.

### Adding a New Display Type

Thanks to the centralized architecture, adding a new visualization page is very straightforward:

1. Create a new page component in `src/pages`.
2. Retrieve the loaded protein data seamlessly using the custom hook: `const { protein } = useProtein();`
3. Add a new route in `src/App.tsx`.

## Getting Started

### Prerequisites

Ensure you have [Node.js](https://nodejs.org/) installed on your machine.

### Installation

1. Navigate to the frontend directory:
   ```bash
   cd frontend
   ```
2. Install the dependencies:
   ```bash
   npm install
   ```

### Running for Development

To start the Vite development server (usually available at `http://localhost:5173`):

```bash
npm run dev
```

### Building for Production

To build the static assets for production:

```bash
npm run build
```

This will output the minified files into the `dist/` directory, which can then be served by any static file server (or hosted alongside the Python backend).
