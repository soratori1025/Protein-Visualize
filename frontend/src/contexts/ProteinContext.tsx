import React, { createContext, useContext, useState, useEffect } from 'react';
import { uploadStructure, getHealth } from '../services/api';
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
  checkHealth: () => Promise<void>;
  
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

  const handleUpload = async (file: File) => {
    setIsUploading(true);
    setStatus(`Parsing ${file.name}...`);
    try {
      const result = await uploadStructure(file);
      setProtein(result);
      setChainId(result.models[0]?.chains[0]?.id ?? null);
      setSelectedResidue(null);
      setSecondaryResult(null);
      setActiveTopologyData(null);
      setChainAnalysis(null);
      setStatus(`${file.name} loaded`);
      
      try {
        localStorage.setItem('protein-cache', JSON.stringify(result));
      } catch (e) {
        console.warn('Protein too large to cache in localStorage');
      }
    } catch (error) {
      setStatus(error instanceof Error ? error.message : 'Upload failed');
    } finally {
      setIsUploading(false);
    }
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
        checkHealth,
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
