export interface SecondaryStructureResidue {
  chain_id: string;
  residue_number: number;
  residue_name: string;
  code: string;
  phi?: number | null;
  psi?: number | null;
  asa?: number | null;
}

export interface SecondaryStructureResult {
  method: string;
  residues: SecondaryStructureResidue[];
}