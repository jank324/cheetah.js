// bundle.js - Simple build script to concatenate and wrap src files into dist/cheetah.js (UMD)

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const srcDir = path.join(__dirname, '../src');
const distDir = path.join(__dirname, '../dist');

// Define dependency order (independent first, dependent later)
const files = [
    'math.js',
    'species.js',
    'beam.js',
    'parameter_beam.js',
    'particle_beam.js',
    'element.js',
    'track_methods.js',
    'drift.js',
    'quadrupole.js',
    'dipole.js',
    'horizontal_corrector.js',
    'vertical_corrector.js',
    'combined_corrector.js',
    'cavity.js',
    'bpm.js',
    'screen.js',
    'aperture.js',
    'solenoid.js',
    'sextupole.js',
    'marker.js',
    'custom_transfer_map.js',
    'segment.js',
    'latticejson.js'
];

const exportsList = [
    'zeros', 'identity', 'transpose', 'matmul', 'matvec', 'clone2d', 'cholesky6',
    'solveTriangularL', 'solveTriangularLT', 'invertTriangularL', 'randn',
    'weightedMean', 'weightedVariance', 'weightedStd', 'weightedCovariance', 'weightedCovarianceMatrix6', 'matchDistributionMoments',
    'si', 'cx_helper', 'log1pdiv', 'si1mdiv', 'sicos1mdiv', 'sipsicos3mdiv', 'cossqrtmcosdivdiff', 'simsidivdiff', 'si2msi2divdiff', 'sqrta2minusbdiva',
    'electron_mass_eV', 'proton_mass_eV', 'deuteron_mass_eV', 'elementary_charge', 'eV_to_kg', 'Species',
    'Beam', 'ParameterBeam', 'ParticleBeam', 'Element',
    'rotationMatrix', 'misalignmentMatrix', 'combinedRotationMisalignmentMatrix',
    'computeRelativisticFactors', 'base_rmatrix', 'base_ttensor', 'drift_matrix',
    'Drift', 'Quadrupole', 'Dipole', 'applyTensorTransform',
    'HorizontalCorrector', 'VerticalCorrector', 'CombinedCorrector',
    'Cavity', 'BPM', 'Screen', 'Aperture', 'Solenoid', 'Sextupole', 'Marker',
    'CustomTransferMap', 'Segment', 'parseLatticeJSON'
];

function bundle() {
    console.log('Bundling cheetah.js...');
    
    let combinedCode = '';
    
    for (const file of files) {
        const filePath = path.join(srcDir, file);
        let content = fs.readFileSync(filePath, 'utf8');
        
        // Strip imports
        content = content.replace(/^import\s+[\s\S]*?;\r?\n/gm, '');
        
        // Strip export prefix (export class, export function, export const, export default)
        content = content.replace(/^export\s+(class|function|const|let|var)\s+/gm, '$1 ');
        // Strip standalone exports (e.g. export { ... })
        content = content.replace(/^export\s+\{[\s\S]*?\};?\r?\n/gm, '');
        
        combinedCode += `// --- File: ${file} ---\n` + content + '\n';
    }
    
    // Wrap in UMD template
    const umdTemplate = `(function (root, factory) {
    if (typeof define === 'function' && define.amd) {
        define([], factory);
    } else if (typeof module === 'object' && module.exports) {
        module.exports = factory();
    } else {
        root.cheetah = factory();
    }
}(typeof self !== 'undefined' ? self : this, function () {
    'use strict';

    ${combinedCode.split('\n').map(line => '    ' + line).join('\n')}

    return {
        ${exportsList.map(exp => `        ${exp}`).join(',\n')}
    };
}));
`;

    if (!fs.existsSync(distDir)) {
        fs.mkdirSync(distDir, { recursive: true });
    }
    
    const outputPath = path.join(distDir, 'cheetah.js');
    fs.writeFileSync(outputPath, umdTemplate, 'utf8');
    console.log(`Successfully bundled to ${outputPath}`);
}

bundle();
