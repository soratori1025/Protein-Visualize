export interface Point3D {
  x: number;
  y: number;
  z: number;
}

export interface RibbonPoint extends Point3D {
  tangent: Point3D;
  normal: Point3D;
  binormal: Point3D;
  left: Point3D;
  right: Point3D;
  t: number;
  residueId?: number;
}

function sub(a: Point3D, b: Point3D): Point3D {
  return { x: a.x - b.x, y: a.y - b.y, z: a.z - b.z };
}

function add(a: Point3D, b: Point3D): Point3D {
  return { x: a.x + b.x, y: a.y + b.y, z: a.z + b.z };
}

function scale(v: Point3D, s: number): Point3D {
  return { x: v.x * s, y: v.y * s, z: v.z * s };
}

function length(v: Point3D): number {
  return Math.sqrt(v.x * v.x + v.y * v.y + v.z * v.z);
}

function normalize(v: Point3D): Point3D {
  const len = length(v);
  if (len < 1e-6) return { x: 0, y: 1, z: 0 };
  return { x: v.x / len, y: v.y / len, z: v.z / len };
}

function cross(a: Point3D, b: Point3D): Point3D {
  return {
    x: a.y * b.z - a.z * b.y,
    y: a.z * b.x - a.x * b.z,
    z: a.x * b.y - a.y * b.x,
  };
}

/**
 * Catmull-Rom spline interpolation between p1 and p2 using control points p0 and p3.
 */
export function catmullRomPoint(p0: Point3D, p1: Point3D, p2: Point3D, p3: Point3D, t: number): Point3D {
  const t2 = t * t;
  const t3 = t2 * t;

  return {
    x: 0.5 * (2 * p1.x + (-p0.x + p2.x) * t + (2 * p0.x - 5 * p1.x + 4 * p2.x - p3.x) * t2 + (-p0.x + 3 * p1.x - 3 * p2.x + p3.x) * t3),
    y: 0.5 * (2 * p1.y + (-p0.y + p2.y) * t + (2 * p0.y - 5 * p1.y + 4 * p2.y - p3.y) * t2 + (-p0.y + 3 * p1.y - 3 * p2.y + p3.y) * t3),
    z: 0.5 * (2 * p1.z + (-p0.z + p2.z) * t + (2 * p0.z - 5 * p1.z + 4 * p2.z - p3.z) * t2 + (-p0.z + 3 * p1.z - 3 * p2.z + p3.z) * t3),
  };
}

/**
 * Derivative (tangent) of Catmull-Rom spline at parameter t.
 */
export function catmullRomTangent(p0: Point3D, p1: Point3D, p2: Point3D, p3: Point3D, t: number): Point3D {
  const t2 = t * t;

  return normalize({
    x: 0.5 * ((-p0.x + p2.x) + 2 * (2 * p0.x - 5 * p1.x + 4 * p2.x - p3.x) * t + 3 * (-p0.x + 3 * p1.x - 3 * p2.x + p3.x) * t2),
    y: 0.5 * ((-p0.y + p2.y) + 2 * (2 * p0.y - 5 * p1.y + 4 * p2.y - p3.y) * t + 3 * (-p0.y + 3 * p1.y - 3 * p2.y + p3.y) * t2),
    z: 0.5 * ((-p0.z + p2.z) + 2 * (2 * p0.z - 5 * p1.z + 4 * p2.z - p3.z) * t + 3 * (-p0.z + 3 * p1.z - 3 * p2.z + p3.z) * t2),
  });
}

/**
 * Generates a smooth Catmull-Rom ribbon band along a series of 3D C-alpha points.
 */
export function generateRibbonSpline(
  points: { id: number; x: number; y: number; z: number }[],
  stepsPerSegment = 5,
  ribbonWidth = 1.2
): RibbonPoint[] {
  if (points.length < 2) return [];

  // Pad endpoints for Catmull-Rom control points
  const pStart = add(points[0], sub(points[0], points[1]));
  const pEnd = add(points[points.length - 1], sub(points[points.length - 1], points[points.length - 2]));
  const padded = [pStart, ...points, pEnd];

  const ribbon: RibbonPoint[] = [];
  let prevNormal: Point3D = { x: 0, y: 1, z: 0 };

  for (let i = 1; i < padded.length - 2; i++) {
    const p0 = padded[i - 1];
    const p1 = padded[i];
    const p2 = padded[i + 1];
    const p3 = padded[i + 2];
    const origResidueId = points[i - 1]?.id;

    for (let step = 0; step < stepsPerSegment; step++) {
      const t = step / stepsPerSegment;
      const pt = catmullRomPoint(p0, p1, p2, p3, t);
      const tangent = catmullRomTangent(p0, p1, p2, p3, t);

      // Compute normal & binormal (Parallel Transport Frame to avoid twisting)
      let normal = cross(tangent, prevNormal);
      if (length(normal) < 1e-4) {
        // Fallback arbitrary vector if tangent parallel to prevNormal
        const up = Math.abs(tangent.z) < 0.9 ? { x: 0, y: 0, z: 1 } : { x: 1, y: 0, z: 0 };
        normal = normalize(cross(tangent, up));
      } else {
        normal = normalize(normal);
      }
      const binormal = normalize(cross(normal, tangent));
      prevNormal = normal;

      // Calculate left and right offset points along binormal
      const left = add(pt, scale(binormal, ribbonWidth / 2));
      const right = add(pt, scale(binormal, -ribbonWidth / 2));

      ribbon.push({
        ...pt,
        tangent,
        normal,
        binormal,
        left,
        right,
        t: (i - 1) + t,
        residueId: origResidueId,
      });
    }
  }

  return ribbon;
}

/**
 * Hydropathy index mapping for amino acid residue types (Kyte-Doolittle scale).
 */
export const hydropathyScale: Record<string, number> = {
  ILE: 4.5, VAL: 4.2, LEU: 3.8, PHE: 2.8, CYS: 2.5, MET: 1.9, ALA: 1.8,
  GLY: -0.4, THR: -0.7, SER: -0.8, TRP: -0.9, TYR: -1.3, PRO: -1.6,
  HIS: -3.2, GLU: -3.5, GLN: -3.5, ASP: -3.5, ASN: -3.5, LYS: -3.9, ARG: -4.5,
};

/**
 * Returns a hex color for a residue based on hydropathy (hydrophobic = orange/red, hydrophilic = cyan/blue).
 */
export function getHydropathyColor(resName: string): string {
  const score = hydropathyScale[resName.toUpperCase()] ?? 0;
  if (score > 1.5) return '#e76f51'; // Hydrophobic (Orange-Red)
  if (score > 0) return '#f4a261';  // Mildly Hydrophobic
  if (score > -2.0) return '#2a9d8f'; // Neutral/Polar (Teal)
  return '#457b9d';                 // Hydrophilic (Blue)
}
