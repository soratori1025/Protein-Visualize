export interface SecondaryStructureResidue {
  chain_id: string;
  residue_number: number;
  /** PDB insertion code ("A" in residue 100A); "" when blank. Older backends omit it. */
  insertion_code?: string;
  residue_name: string;
  code: string;
  /** null when undefined (chain ends / breaks — DSSP and STRIDE print 360). */
  phi?: number | null;
  psi?: number | null;
  /** Relative solvent accessibility 0..1 for both DSSP and STRIDE; null when unknown. */
  asa?: number | null;
}

export interface SecondaryStructureResult {
  method: string;
  residues: SecondaryStructureResidue[];
}

/** Region types. The calculated endpoint can also return 'Signal'; other UniProt
 *  feature types may pass through, hence the open string. */
export type TopologyRegionType =
  | 'Transmembrane'
  | 'Topological domain'
  | 'Intramembrane'
  | 'Signal'
  | (string & {});

export type TopologySide = 'membrane' | 'Cytoplasmic' | 'Extracellular';
export type TopologySS = 'Helix' | 'Strand' | 'Loop' | 'Coil';

export interface UniProtTopologyRegion {
  type: TopologyRegionType;
  /** Residue numbers. UniProt source: UniProt sequence positions.
   *  Calculated source: author numbering of the uploaded structure. */
  start: number;
  end: number;
  description: string;
  /** UniProt feature name. The calculated endpoint does not send it. */
  name?: string;
  /** Sent by the calculated endpoint — use these instead of parsing `description`. */
  side?: TopologySide | null;
  ss?: TopologySS | null;
  start_icode?: string | null;
  end_icode?: string | null;
}

export interface UniProtTopologyData {
  uniprot_id: string;
  protein_name: string;
  gene_name: string;
  organism: string;
  regions: UniProtTopologyRegion[];
}

/** Response of /api/secondary-structure/predict-topology/{filename}. */
export interface CalculatedTopologyData extends UniProtTopologyData {
  membrane_score: number;
  labeler: string;
  membrane_normal?: number[] | null;
  parameters_used?: Record<string, number | string> | null;
  /** Chain the backend analysed. */
  chain_id?: string | null;
  /** Things to show the user: SS tool missing, numbering re-mapped, features skipped… */
  warnings?: string[];
}

export type TopologyData = UniProtTopologyData | CalculatedTopologyData;

export function isCalculatedTopology(data: TopologyData | null | undefined): data is CalculatedTopologyData {
  return !!data && typeof (data as CalculatedTopologyData).labeler === 'string';
}

/** "100A" style residue label. */
export function residueLabel(residueNumber: number, insertionCode?: string | null): string {
  return `${residueNumber}${(insertionCode ?? '').trim()}`;
}

/** Sort key that keeps insertion codes in order: 100 < 100A < 100B < 101. */
export function compareResidues(
  a: { residue_number: number; insertion_code?: string },
  b: { residue_number: number; insertion_code?: string }
): number {
  return a.residue_number - b.residue_number || (a.insertion_code ?? '').localeCompare(b.insertion_code ?? '');
}
