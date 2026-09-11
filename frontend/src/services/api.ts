import type { ProteinUpload } from '../types/protein';
import type { ChainAnalysis } from '../types/analysis';
import type { SecondaryStructureResult } from '../types/secondaryStructure';

const API_URL = 'http://127.0.0.1:8000';

export async function uploadStructure(file: File): Promise<ProteinUpload> {
  const body = new FormData();
  body.append('file', file);
  const response = await fetch(`${API_URL}/api/structure/upload`, {
    method: 'POST',
    body,
  });

  if (!response.ok) {
    throw new Error(`Upload failed (${response.status})`);
  }
  return response.json() as Promise<ProteinUpload>;
}

export async function getHealth(): Promise<Record<string, string>> {
  const response = await fetch(`${API_URL}/api/health`);
  if (!response.ok) throw new Error('Backend unavailable');
  return response.json();
}

export async function runSecondaryStructure(filename: string, method: 'DSSP' | 'STRIDE'): Promise<SecondaryStructureResult> {
  const response = await fetch(`${API_URL}/api/secondary-structure/${encodeURIComponent(filename)}?method=${method}`, { method: 'POST' });
  const data = await response.json() as SecondaryStructureResult & { detail?: string };
  if (!response.ok) throw new Error(data.detail ?? 'Secondary-structure analysis failed');
  return data;
}

export async function analyzeChain(filename: string, chainId: string, contactCutoff = 8): Promise<ChainAnalysis> {
  const response = await fetch(`${API_URL}/api/analysis/${encodeURIComponent(filename)}/chain/${encodeURIComponent(chainId)}?contact_cutoff=${contactCutoff}`);
  const data = await response.json() as ChainAnalysis | { detail?: string };
  if (!response.ok) throw new Error('detail' in data ? data.detail : 'Analysis failed');
  return data as ChainAnalysis;
}

export async function getSecondaryCapabilities(): Promise<Record<string, { available: boolean; executable: string }>> {
  const response = await fetch(`${API_URL}/api/secondary-structure/capabilities`);
  if (!response.ok) throw new Error('Could not read secondary-structure capabilities');
  return response.json();
}