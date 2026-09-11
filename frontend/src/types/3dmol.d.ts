declare module '3dmol' {
  export interface Viewer {
    addModel(data: string, format: string): unknown;
    setStyle(selection: Record<string, unknown>, style: Record<string, unknown>): void;
    setClickable(selection: Record<string, unknown>, enabled: boolean, callback: (atom: { resi?: number; chain?: string }) => void): void;
    zoomTo(selection?: Record<string, unknown>): void;
    render(): void;
    clear(): void;
    setBackgroundColor(color: string): void;
    resize(): void;
  }

  export function createViewer(element: HTMLElement, config?: Record<string, unknown>): Viewer;
}