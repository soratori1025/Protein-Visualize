import { useState, useCallback, useEffect } from 'react';
import { API_URL, runSecondaryStructure } from '../../../services/api';
import type { CalculatedTopologyData, UniProtTopologyData, SecondaryStructureResult } from '../../../types/secondaryStructure';
import type { Chain } from '../../../types/protein';

export function useTopologyData(
  filename?: string | null,
  chain?: Chain,
  uniprotId?: string | null,
  triggerTmRecalc?: number,
  distinguishTurns: boolean = false
) {
  const [uniprotData, setUniprotData] = useState<UniProtTopologyData | null>(null);
  const [loadingUniProt, setLoadingUniProt] = useState<boolean>(false);
  const [uniprotError, setUniprotError] = useState<string | null>(null);

  const [calculatedData, setCalculatedData] = useState<CalculatedTopologyData | null>(null);
  const [loadingCalculated, setLoadingCalculated] = useState<boolean>(false);
  const [calculatedError, setCalculatedError] = useState<string | null>(null);

  const fetchUniProtTopology = async (uniprotIdToFetch: string) => {
    if (!uniprotIdToFetch.trim()) return;
    setLoadingUniProt(true);
    setUniprotError(null);
    try {
      const response = await fetch(`${API_URL}/api/secondary-structure/uniprot/${uniprotIdToFetch.trim()}`);
      if (!response.ok) {
        const errData = await response.json().catch(() => ({}));
        throw new Error(errData.detail || `Could not load UniProt entry ${uniprotIdToFetch}`);
      }
      setUniprotData((await response.json()) as UniProtTopologyData);
    } catch (err: any) {
      setUniprotError(err.message || 'Could not load UniProt data');
      setUniprotData(null);
    } finally {
      setLoadingUniProt(false);
    }
  };

  const fetchCalculatedTopology = useCallback(async (
    filenameToFetch: string,
    tm: string,
    ss: string,
    cUni: string,
    chainIdToFetch?: string,
    tmThickness?: string,
    tmMinMembraneScore?: string,
    tmTreatTurnAsHelix?: boolean,
    tmMinCrossSpan?: string,
    tmFullCrossFrac?: string
  ) => {
    if (!filenameToFetch.trim()) return;
    setLoadingCalculated(true);
    setCalculatedError(null);
    try {
      const qp = new URLSearchParams({ tm_algo: tm, ss_algo: ss, flow_type: 'consensus' });
      if (chainIdToFetch) qp.set('chain_id', chainIdToFetch);
      if (tm === 'uniprot_api' && cUni.trim()) qp.set('uniprot_id', cUni.trim().toUpperCase());
      
      const setNumber = (key: string, raw?: string, integer = false) => {
        if (!raw || !raw.trim()) return;
        const value = integer ? parseInt(raw, 10) : parseFloat(raw);
        if (!isNaN(value)) qp.set(key, String(value));
      };
      
      setNumber('thickness', tmThickness);
      setNumber('min_membrane_score', tmMinMembraneScore);
      if (tmTreatTurnAsHelix || distinguishTurns) {
        qp.set('treat_turn_as_helix', 'true');
      }
      if (ss !== 'none') {
        setNumber('min_cross_span', tmMinCrossSpan);
        setNumber('full_cross_frac', tmFullCrossFrac);
      }

      const response = await fetch(
        `${API_URL}/api/secondary-structure/predict-topology/${encodeURIComponent(filenameToFetch.trim())}?${qp.toString()}`
      );
      if (!response.ok) {
        const errData = await response.json().catch(() => ({}));
        throw new Error(errData.detail || `Could not compute topology for ${filenameToFetch}`);
      }
      setCalculatedData((await response.json()) as CalculatedTopologyData);
    } catch (err: any) {
      setCalculatedError(err.message || 'Could not compute topology');
      setCalculatedData(null);
    } finally {
      setLoadingCalculated(false);
    }
  }, [distinguishTurns]);

  useEffect(() => {
    if (uniprotId) {
      fetchUniProtTopology(uniprotId);
    }
  }, [uniprotId]);

  return {
    uniprotData,
    setUniprotData,
    loadingUniProt,
    uniprotError,
    fetchUniProtTopology,
    calculatedData,
    setCalculatedData,
    loadingCalculated,
    calculatedError,
    fetchCalculatedTopology
  };
}
