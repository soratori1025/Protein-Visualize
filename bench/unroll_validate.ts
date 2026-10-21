/**
 * Validation harness for the 2-D unrolling (`spreadStructure.ts`).
 *
 * Node 24 can run TypeScript directly (type stripping), so the production module is
 * exercised exactly as the browser does - no re-implementation, no mock.
 *
 *   node --experimental-strip-types bench/unroll_validate.ts <job.json> <out.json>
 *
 * The job file is
 *   { "structure": "<path to a PDB or mmCIF file>",
 *     "consensusMap": [ { "residue_number": 12, "insertion_code": "", "aa": "M",
 *                         "label": "TM_in", "crossing": 1, "part": null } ],
 *     "membrane": { "normal": [0,1,0], "half_thickness": 15 } }
 *
 * and the output file is the `SpreadReport` of the unrolling, plus the report of the
 * rigid 'aligned' mode for comparison.
 */
import { readFileSync, writeFileSync } from 'node:fs';
import { spreadStructureWithReport } from '../frontend/src/components/structure/spreadStructure.ts';

const [, , jobPath, outPath] = process.argv;
if (!jobPath || !outPath) {
  console.error('usage: node --experimental-strip-types bench/unroll_validate.ts <job.json> <out.json>');
  process.exit(2);
}

const job = JSON.parse(readFileSync(jobPath, 'utf8'));
const structure = readFileSync(job.structure, 'utf8');

const unrolled = spreadStructureWithReport(structure, job.consensusMap, {
  mode: 'unrolled',
  chainId: job.chainId ?? null,
  membrane: job.membrane ?? null,
  residues: job.residues ?? null,
});

const aligned = spreadStructureWithReport(structure, job.consensusMap, {
  mode: 'aligned',
  chainId: job.chainId ?? null,
  membrane: job.membrane ?? null,
  residues: job.residues ?? null,
});

writeFileSync(outPath, JSON.stringify({ unrolled: unrolled.report, aligned: aligned.report }, null, 2));
