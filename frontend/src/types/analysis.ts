export interface ChainAnalysis {
  input: { filename: string; chain_id: string };
  result: {
    chain_id: string;
    algorithm: string;
    version: string;
    parameters: { contact_cutoff: number };
    sequence: { length: number; composition: Record<string, number>; unknown_residues: number };
    physicochemical: {
      molecular_weight: number | null;
      isoelectric_point: number | null;
      aromaticity: number | null;
      instability_index: number | null;
      mean_hydropathy: number | null;
    };
    geometry: { ca_atoms: number; contacts: number; contact_density: number };
    contacts: Array<{ source: number; target: number; distance: number }>;
  };
}