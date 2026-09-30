import { KeyframeConfig } from '../components/viewer/StoryboardViewer';

export const DEMO_KEYFRAMES: KeyframeConfig[] = [
  {
    id: 'intro',
    title: '1. Overall Structure',
    description: "This is a comprehensive view of the molecule's 3D spatial structure. The Ribbon representation clearly shows major folds, including alpha helices and beta sheets, allowing you to see the overall shape of the molecule.",
    style: 'ribbon',
    colorScheme: 'chain',
    highlights: []
  },
  {
    id: 'pocket',
    title: '2. Core and Binding Pockets',
    description: 'Switching to the Surface representation reveals the spatial mesh of the molecular surface. This mode is typically used to investigate grooves and binding pockets - where substrates or ligands attach.',
    style: 'surface',
    colorScheme: 'chain',
    highlights: []
  },
  {
    id: 'interaction',
    title: '3. Crucial Residues',
    description: 'Returning to the Ribbon representation, we highlight important residues. In this example, we zoom into specific residues on the first chain.',
    focus: { resi: [10, 45] },
    style: 'ribbon',
    colorScheme: 'chain',
    highlights: [
      { chain: 'A', resi: [10, 45] }
    ]
  }
];
