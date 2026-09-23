import React, { createContext, useContext, useState, useEffect } from 'react';
import { uploadStructure, getHealth } from '../services/api';
import type { ProteinUpload } from '../types/protein';

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
}

const ProteinContext = createContext<ProteinContextType | undefined>(undefined);

export function ProteinProvider({ children }: { children: React.ReactNode }) {
  const [protein, setProtein] = useState<ProteinUpload | null>(null);
  const [chainId, setChainId] = useState<string | null>(null);
  const [selectedResidue, setSelectedResidue] = useState<number | null>(null);
  const [status, setStatus] = useState('Ready for a structure file');
  const [health, setHealth] = useState('API status unknown');
  const [isUploading, setIsUploading] = useState(false);

  const checkHealth = async () => {
    try {
      const result = await getHealth();
      setHealth(`${result.service} · ${result.version}`);
    } catch {
      setHealth('Backend unavailable');
    }
  };

  useEffect(() => {
    void checkHealth();
  }, []);

  const handleUpload = async (file: File) => {
    setIsUploading(true);
    setStatus(`Parsing ${file.name}...`);
    try {
      const result = await uploadStructure(file);
      setProtein(result);
      setChainId(result.models[0]?.chains[0]?.id ?? null);
      setSelectedResidue(null);
      setStatus(`${file.name} loaded`);
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
