import React, { createContext, useContext, useState, useEffect } from 'react';
import { uploadStructure, getHealth, fetchRemoteStructure } from '../services/api';
import type { ProteinUpload } from '../types/protein';
import type { SecondaryStructureResult, UniProtTopologyData } from '../types/secondaryStructure';
import type { ChainAnalysis } from '../types/analysis';

interface ProteinContextType {
  protein: ProteinUpload | null;
  chainId: string | null;
  setChainId: (id: string | null) => void;
  selectedResidue: number | null;
  setSelectedResidue: (id: number | null) => void;
  status: string;
  setStatus: (status: string) => void;
  health: string;
  isUploading: boolean;
  handleUpload: (file: File) => Promise<void>;
  handleFetchRemote: (id: string, type: 'PDB' | 'UniProt') => Promise<void>;
  handleClear: () => void;
  checkHealth: () => Promise<void>;
  
  fetchInputId: string;
  setFetchInputId: (id: string) => void;
  fetchInputType: 'PDB' | 'UniProt';
  setFetchInputType: (type: 'PDB' | 'UniProt') => void;

  secondaryResult: SecondaryStructureResult | null;
  setSecondaryResult: (res: SecondaryStructureResult | null) => void;
  activeTopologyData: UniProtTopologyData | null;
  setActiveTopologyData: (data: UniProtTopologyData | null) => void;
  chainAnalysis: ChainAnalysis | null;
  setChainAnalysis: (data: ChainAnalysis | null) => void;
}

const ProteinContext = createContext<ProteinContextType | undefined>(undefined);

export function ProteinProvider({ children }: { children: React.ReactNode }) {
  const [protein, setProtein] = useState<ProteinUpload | null>(null);
  const [chainId, setChainId] = useState<string | null>(null);
  const [selectedResidue, setSelectedResidue] = useState<number | null>(null);
  const [status, setStatus] = useState('Ready for a structure file');
  const [health, setHealth] = useState('API status unknown');
  const [isUploading, setIsUploading] = useState(false);
  
  const [fetchInputId, setFetchInputId] = useState('');
  const [fetchInputType, setFetchInputType] = useState<'PDB' | 'UniProt'>('PDB');

  const [secondaryResult, setSecondaryResult] = useState<SecondaryStructureResult | null>(null);
  const [activeTopologyData, setActiveTopologyData] = useState<UniProtTopologyData | null>(null);
  const [chainAnalysis, setChainAnalysis] = useState<ChainAnalysis | null>(null);

  const checkHealth = async () => {
    try {
      const result = await getHealth();
      setHealth(`${result.service} · ${result.version}`);
    } catch {
      setHealth('Backend unavailable');
    }
  };

  // Restore protein from cache on initial load
  useEffect(() => {
    void checkHealth();
    
    const cachedProtein = localStorage.getItem('protein-cache');
    if (cachedProtein) {
      try {
        const parsed = JSON.parse(cachedProtein) as ProteinUpload;
        setProtein(parsed);
        setChainId(parsed.models[0]?.chains[0]?.id ?? null);
        setStatus(`Restored ${parsed.filename} from cache`);
      } catch (e) {
        console.error('Failed to parse cached protein:', e);
        localStorage.removeItem('protein-cache');
      }
    }
  }, []);

  const _resetStateForNewProtein = (result: ProteinUpload, filename: string) => {
    setProtein(result);
    setChainId(result.models[0]?.chains[0]?.id ?? null);
    setSelectedResidue(null);
    setSecondaryResult(null);
    setActiveTopologyData(null);
    setChainAnalysis(null);
    setStatus(`${filename} loaded`);
    
    try {
      localStorage.setItem('protein-cache', JSON.stringify(result));
    } catch (e) {
      console.warn('Protein too large to cache in localStorage');
    }
  };

  const handleUpload = async (file: File) => {
    setIsUploading(true);
    setStatus(`Parsing ${file.name}...`);
    try {
      const result = await uploadStructure(file);
      _resetStateForNewProtein(result, file.name);
    } catch (error) {
      setStatus(error instanceof Error ? error.message : 'Upload failed');
    } finally {
      setIsUploading(false);
    }
  };

  const handleFetchRemote = async (id: string, type: 'PDB' | 'UniProt') => {
    setIsUploading(true);
    setStatus(`Fetching ${type} ID ${id}...`);
    try {
      const result = await fetchRemoteStructure(id, type);
      _resetStateForNewProtein(result, result.filename);
    } catch (error) {
      const errMsg = error instanceof Error ? error.message : 'Fetch failed';
      setStatus(errMsg);
      window.alert(`Error! Cannot find protein with this ${type} ID "${id}".\nDetails: ${errMsg}`);
    } finally {
      setIsUploading(false);
    }
  };

  const handleClear = () => {
    setProtein(null);
    setChainId(null);
    setSecondaryResult(null);
    setActiveTopologyData(null);
    setChainAnalysis(null);
    setSelectedResidue(null);
    setFetchInputId('');
    setStatus('Ready for a structure file');
    localStorage.removeItem('protein-cache');
  };

  return (
    <ProteinContext.Provider
      value={{
        protein,
        chainId,
        setChainId,
        selectedResidue,
        setSelectedResidue,
        status,
        setStatus,
        health,
        isUploading,
        handleUpload,
        handleFetchRemote,
        handleClear,
        checkHealth,
        fetchInputId,
        setFetchInputId,
        fetchInputType,
        setFetchInputType,
        secondaryResult,
        setSecondaryResult,
        activeTopologyData,
        setActiveTopologyData,
        chainAnalysis,
        setChainAnalysis,
      }}
    >
      {children}
    </ProteinContext.Provider>
  );
}

export function useProtein() {
  const context = useContext(ProteinContext);
  if (context === undefined) {
    throw new Error('useProtein must be used within a ProteinProvider');
  }
  return context;
}
