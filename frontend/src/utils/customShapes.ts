import * as $3Dmol from '3dmol';
import type { ColorScheme } from '../types/viewer';
import { getHydropathyColor } from './ribbonSpline';

interface Point { x: number; y: number; z: number }

function add(a: Point, b: Point): Point { return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z }; }
function sub(a: Point, b: Point): Point { return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z }; }
function scale(a: Point, s: number): Point { return { x: a.x * s, y: a.y * s, z: a.z * s }; }
function cross(a: Point, b: Point): Point {
  return {
    x: a.y * b.z - a.z * b.y,
    y: a.z * b.x - a.x * b.z,
    z: a.x * b.y - a.y * b.x
  };
}
function dot(a: Point, b: Point): number { return a.x * b.x + a.y * b.y + a.z * b.z; }
function normalize(a: Point): Point {
  const len = Math.sqrt(dot(a, a));
  if (len === 0) return { x: 0, y: 0, z: 0 };
  return scale(a, 1 / len);
}

const ssColors = {
  h: '#ff0080', // Helix (Pinkish/Red)
  s: '#ffc800', // Sheet (Yellow)
  c: '#00ffcc', // Coil (Cyan/Green)
};

export function drawCustomPipesAndPlanks(
  viewer: $3Dmol.Viewer,
  chainId: string,
  baseColor: string,
  colorScheme: ColorScheme
) {
  const v = viewer as any;
  // Extract all CA atoms for this chain
  const atoms = v.selectedAtoms({ chain: chainId, atom: 'CA' });
  if (atoms.length === 0) return;

  // Group continuous atoms by secondary structure
  const segments: { ss: string; atoms: typeof atoms }[] = [];
  let currentSegment: { ss: string; atoms: typeof atoms } | null = null;

  for (const atom of atoms) {
    const ss = (atom.ss as string) || 'c';
    const cleanSs = ss === 'h' ? 'h' : ss === 's' ? 's' : 'c';

    if (!currentSegment || currentSegment.ss !== cleanSs) {
      // Overlap by 1 atom so the rendering connects nicely
      let prevAtom = null;
      if (currentSegment && currentSegment.atoms.length > 0) {
        prevAtom = currentSegment.atoms[currentSegment.atoms.length - 1];
      }
      currentSegment = { ss: cleanSs, atoms: prevAtom ? [prevAtom] : [] };
      segments.push(currentSegment);
    }
    currentSegment.atoms.push(atom);
  }

  const getColor = (atom: any) => {
    if (colorScheme === 'ss') return atom.ss === 'h' ? ssColors.h : atom.ss === 's' ? ssColors.s : ssColors.c;
    if (colorScheme === 'hydropathy') return getHydropathyColor(atom.resn || atom.resname || 'ALA');
    return baseColor;
  };

  segments.forEach((seg) => {
    if (seg.atoms.length < 2) return;
    
    // Predominant color for the segment
    const midAtom = seg.atoms[Math.floor(seg.atoms.length / 2)];
    const color = getColor(midAtom);

    if (seg.ss === 'h') {
      // Draw Cylinder (Pipe) for Helix
      const start = seg.atoms[0];
      const end = seg.atoms[seg.atoms.length - 1];

      // Approximate helix axis with start and end
      v.addCylinder({
        start: { x: start.x, y: start.y, z: start.z },
        end: { x: end.x, y: end.y, z: end.z },
        radius: 1.2,
        color: color,
        fromCap: true,
        toCap: true,
      });

    } else if (seg.ss === 's') {
      // Draw Plank (Box) for Beta Strand
      const start = { x: seg.atoms[0].x, y: seg.atoms[0].y, z: seg.atoms[0].z };
      const end = { x: seg.atoms[seg.atoms.length - 1].x, y: seg.atoms[seg.atoms.length - 1].y, z: seg.atoms[seg.atoms.length - 1].z };
      
      const dir = sub(end, start); // Strand direction vector
      
      // Calculate normal of the beta sheet using zigzag CA pattern
      let avgNormal = { x: 0, y: 0, z: 0 };
      for (let i = 0; i < seg.atoms.length - 2; i++) {
        const p1 = { x: seg.atoms[i].x, y: seg.atoms[i].y, z: seg.atoms[i].z };
        const p2 = { x: seg.atoms[i+1].x, y: seg.atoms[i+1].y, z: seg.atoms[i+1].z };
        const p3 = { x: seg.atoms[i+2].x, y: seg.atoms[i+2].y, z: seg.atoms[i+2].z };
        const v1 = sub(p2, p1);
        const v2 = sub(p3, p2);
        const n = cross(v1, v2);
        if (i > 0 && dot(avgNormal, n) < 0) {
           avgNormal = add(avgNormal, scale(n, -1));
        } else {
           avgNormal = add(avgNormal, n);
        }
      }
      
      avgNormal = normalize(avgNormal);
      if (avgNormal.x === 0 && avgNormal.y === 0 && avgNormal.z === 0) {
        avgNormal = normalize(cross(dir, { x: 1, y: 0, z: 0 }));
        if (avgNormal.x === 0 && avgNormal.y === 0 && avgNormal.z === 0) {
           avgNormal = normalize(cross(dir, { x: 0, y: 1, z: 0 }));
        }
      }

      // Depth vector orthogonal to direction and normal
      const depthVec = normalize(cross(dir, avgNormal));

      const width = 2.4; 
      const thickness = 0.6;
      
      const w = scale(avgNormal, width);
      const h = scale(depthVec, thickness);
      
      const corner = sub(sub(start, scale(w, 0.5)), scale(h, 0.5));

      v.addBox({
        corner: corner,
        dimensions: {
          w: w,
          h: h,
          d: dir
        },
        color: color
      });

    } else {
      // Draw Coil (Tube)
      const points = seg.atoms.map((a: any) => ({ x: a.x, y: a.y, z: a.z }));
      v.addCurve({
        points: points,
        radius: 0.25,
        color: color,
        smooth: 10,
        fromCap: 2,
        toCap: 2
      });
    }
  });
}
