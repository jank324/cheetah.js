// run_tests.js - Test suite for validating cheetah.js outputs against Python reference

import fs from 'fs';
import path from 'path';
import { fileURLToPath } from 'url';
import * as cheetah from '../src/index.js';

const __filename = fileURLToPath(import.meta.url);
const __dirname = path.dirname(__filename);

const testVectorsPath = path.join(__dirname, 'test_vectors.json');
const testCases = JSON.parse(fs.readFileSync(testVectorsPath, 'utf8'));

let passedTests = 0;
let failedTests = 0;

function assertClose(actual, expected, tolerance = 1e-7, pathStr = '') {
    if (expected === null || expected === undefined) {
        if (actual !== null && actual !== undefined) {
            throw new Error(`At ${pathStr}: expected null/undefined, got ${actual}`);
        }
        return;
    }

    if (typeof expected === 'number') {
        if (isNaN(expected)) {
            if (!isNaN(actual)) {
                throw new Error(`At ${pathStr}: expected NaN, got ${actual}`);
            }
            return;
        }
        const diff = Math.abs(actual - expected);
        if (diff > tolerance) {
            throw new Error(`At ${pathStr}: expected ${expected}, got ${actual} (diff: ${diff} > tol ${tolerance})`);
        }
    } else if (Array.isArray(expected) || expected instanceof Float64Array) {
        if (!actual) {
            throw new Error(`At ${pathStr}: expected array, got ${actual}`);
        }
        if (actual.length !== expected.length) {
            throw new Error(`At ${pathStr}: length mismatch, expected ${expected.length}, got ${actual.length}`);
        }
        for (let i = 0; i < expected.length; i++) {
            assertClose(actual[i], expected[i], tolerance, `${pathStr}[${i}]`);
        }
    } else if (typeof expected === 'object') {
        for (const k in expected) {
            assertClose(actual[k], expected[k], tolerance, `${pathStr}.${k}`);
        }
    }
}

console.log('Starting validation tests...');

for (const tc of testCases) {
    const name = `${tc.class} (Case ${tc.id}, ${JSON.stringify(tc.params)})`;
    try {
        // Instantiate the element in JS
        let element;
        const clsName = tc.class;
        const params = tc.params;
        const length = params.length !== undefined ? params.length : 0.0;
        const opts = { ...params };
        delete opts.length;

        if (clsName === "Marker") {
            element = new cheetah.Marker(tc.id.toString());
        } else if (clsName === "BPM") {
            element = new cheetah.BPM(opts);
        } else if (clsName === "Screen") {
            element = new cheetah.Screen(opts);
        } else if (clsName === "Aperture") {
            element = new cheetah.Aperture(opts);
        } else {
            const Cls = cheetah[clsName];
            element = new Cls(length, opts);
        }

        // Test ParameterBeam tracking if available
        if (tc.outgoing_parameter_beam) {
            const ip = tc.incoming_parameter_beam;
            const species = new cheetah.Species(ip.species);
            const inBeam = new cheetah.ParameterBeam(
                new Float64Array(ip.mu),
                ip.cov.map(row => new Float64Array(row)),
                ip.energy,
                ip.total_charge,
                ip.s,
                species
            );

            const outBeam = element.track(inBeam);

            const op = tc.outgoing_parameter_beam;
            assertClose(outBeam.mu, op.mu, 1e-7, 'ParameterBeam.mu');
            assertClose(outBeam.cov, op.cov, 1e-7, 'ParameterBeam.cov');
            assertClose(outBeam.energy, op.energy, 1e-7, 'ParameterBeam.energy');
            assertClose(outBeam.s, op.s, 1e-7, 'ParameterBeam.s');
        }

        // Test ParticleBeam tracking if available
        if (tc.outgoing_particle_beam) {
            const ip = tc.incoming_particle_beam;
            const species = new cheetah.Species(ip.species);
            const inBeam = new cheetah.ParticleBeam(
                ip.particles.map(row => new Float64Array(row)),
                ip.energy,
                new Float64Array(ip.particle_charges),
                new Float64Array(ip.survival_probabilities),
                ip.s,
                species
            );

            const outBeam = element.track(inBeam);

            const op = tc.outgoing_particle_beam;
            
            // Compare particles
            assertClose(outBeam.particles, op.particles, 1e-7, 'ParticleBeam.particles');
            // Compare survival probabilities
            assertClose(outBeam.survival_probabilities, op.survival_probabilities, 1e-7, 'ParticleBeam.survival_probabilities');
            // Compare energy
            assertClose(outBeam.energy, op.energy, 1e-7, 'ParticleBeam.energy');
            // Compare position
            assertClose(outBeam.s, op.s, 1e-7, 'ParticleBeam.s');
        }

        // Special Element Checks
        if (clsName === "BPM" && tc.extra) {
            assertClose(element.reading, tc.extra.reading, 1e-7, 'BPM.reading');
        }

        console.log(`\x1b[32m✔ PASS\x1b[0m ${name}`);
        passedTests++;
    } catch (e) {
        console.error(`\x1b[31m✘ FAIL\x1b[0m ${name}`);
        console.error(`       Error: ${e.message}`);
        failedTests++;
    }
}

console.log(`\nTest results: ${passedTests} passed, ${failedTests} failed.`);
if (failedTests > 0) {
    process.exit(1);
} else {
    console.log('\x1b[32mAll tests passed successfully!\x1b[0m');
}
