export interface Atom {
  name: string;
  element: string;
  x: number;
  y: number;
  z: number;
}

export interface Residue {
  id: number;
  name: string;
  atoms: Atom[];
}

export interface Chain {
  id: string;
  sequence: string;
  residue_count: number;
  residues: Residue[];
}

export interface ProteinModel {
  id: number;
  chains: Chain[];
}

export interface ProteinUpload {
  filename: string;
  saved_path: string;
  uniprot_id?: string | null;
  models: ProteinModel[];
}