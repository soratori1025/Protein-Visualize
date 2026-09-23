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

export interface UniProtTopologyRegion {
  type: 'Transmembrane' | 'Topological domain' | 'Intramembrane';
  start: number;
  end: number;
  description: string;
  name: string;
  start_icode?: string | null;
  end_icode?: string | null;
}

export interface UniProtTopologyData {
  uniprot_id: string;
  protein_name: string;
  gene_name: string;
  organism: string;
  regions: UniProtTopologyRegion[];
  chain_id?: string | null;
  warnings?: string[];
}