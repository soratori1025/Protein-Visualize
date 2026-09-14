export type RepresentationStyle = 'ribbon' | 'tube' | 'stick' | 'sphere' | 'line' | 'pipesAndPlanks';

export type ColorScheme = 'chain' | 'ss' | 'hydropathy' | 'bFactor';

export interface ViewerOptions {
  style: RepresentationStyle;
  colorScheme: ColorScheme;
  showHetAtoms: boolean;
}
