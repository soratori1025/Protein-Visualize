export type RepresentationStyle = 'ribbon' | 'stick' | 'sphere' | 'line' | 'pipesAndPlanks';

export type ColorScheme = 'chain' | 'ss' | 'hydropathy' | 'bFactor' | 'consensus';

export interface ViewerOptions {
  style: RepresentationStyle;
  colorScheme: ColorScheme;
  showHetAtoms: boolean;
}
