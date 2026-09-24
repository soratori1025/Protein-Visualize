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
/** 'Irregular' = in the membrane but not one regular helix/strand (mixed elements,
 *  no SS, or SS and membrane evidence disagree). */
export type TopologySS = 'Helix' | 'Strand' | 'Irregular' | 'Loop' | 'Coil';
export type Confidence = 'high' | 'medium' | 'low';

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
  /** Agreement of TM block, membrane geometry and DSSP/STRIDE (heuristic). */
  confidence?: Confidence | null;
  /** Transmembrane regions of the calculated endpoint: 1-based crossing number.
   *  A discontinuous crossing (TM1a/1b …) comes as several regions sharing it:
   *  part 'a', the unwound stretch (part null, "Transmembrane Unwound"), part 'b'. */
  crossing?: number | null;
  part?: 'a' | 'b' | null;
  /** Membrane role, separate from secondary structure (`ss`). */
  membrane_role?:
    | 'TM_CROSSING'
    | 'BROKEN_TM'
    | 'UNWOUND'
    | 'REENTRANT'
    | 'INTERFACIAL'
    | 'EXTRAMEMBRANE'
    | 'SIGNAL'
    | null;
  /** TM1, TM1a, TM1 unwound, EL2, EL3a, IL1, N-term, C-term, RE1 … */
  topology_label?: string | null;
  parent_tm?: string | null;
  /** Evidence for each fragment junction inside a crossing (why it is / is not a/b). */
  transitions?: TransitionEvidence[] | null;
}

export interface TransitionEvidence {
  classification: 'CONTINUOUS' | 'BROKEN_TM' | 'AMBIGUOUS' | 'TWO_TM' | 'REENTRANT' | 'SEPARATE';
  gap_residues: number;
  gap_in_membrane: boolean;
  gap_in_core: boolean;
  same_orientation: boolean;
  a_tilt_deg?: number | null;
  b_tilt_deg?: number | null;
  single_crossing: boolean;
  a_full_cross: boolean;
  b_full_cross: boolean;
  reason: string;
}

/** One row of the backend's residue-level evidence matrix. */
export interface ResidueAnnotation {
  index: number;
  residue_number: number;
  insertion_code?: string | null;
  aa: string;
  ss_raw?: string | null;
  ss?: 'H' | 'E' | 'C' | null;
  tm_evidence: string;
  /** Signed distance from the bilayer mid-plane (Å). */
  depth?: number | null;
  zone?: 'CORE' | 'EDGE' | 'OUT' | null;
  plddt?: number | null;
  label: string;
  confidence?: Confidence | null;
}

export interface MembranePlacement {
  source: string;
  normal: number[];
  half_thickness: number;
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
  domain_type?: 'alpha_helical' | 'beta_barrel' | 'beta' | 'mixed' | 'irregular' | null;
  membrane?: MembranePlacement | null;
  residues?: ResidueAnnotation[] | null;
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
