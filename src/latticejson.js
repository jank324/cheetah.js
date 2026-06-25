// latticejson.js - LatticeJSON parser

import { Drift } from "./drift.js";
import { Quadrupole } from "./quadrupole.js";
import { Dipole } from "./dipole.js";
import { HorizontalCorrector, VerticalCorrector, CombinedCorrector } from "./corrector.js";
import { Cavity } from "./cavity.js";
import { BPM } from "./bpm.js";
import { Screen } from "./screen.js";
import { Aperture } from "./aperture.js";
import { Solenoid } from "./solenoid.js";
import { Sextupole } from "./sextupole.js";
import { Marker } from "./marker.js";
import { Segment } from "./segment.js";

const classMap = {
    Drift, Quadrupole, Dipole, 
    HorizontalCorrector, VerticalCorrector, CombinedCorrector,
    Cavity, BPM, Screen, Aperture, Solenoid, Sextupole, Marker
};

export function parseLatticeJSON(jsonStr) {
    const data = typeof jsonStr === "string" ? JSON.parse(jsonStr) : jsonStr;
    const elementsData = data.elements;
    const latticesData = data.lattices;
    const rootName = data.root;

    function parseElement(name) {
        if (!elementsData[name]) {
            throw new Error(`Element '${name}' not defined in elements.`);
        }
        const [clsName, params] = elementsData[name];
        const Cls = classMap[clsName];
        if (!Cls) {
            console.warn(`Unsupported element class '${clsName}', falling back to Marker.`);
            return new Marker(name);
        }

        const length = params.length !== undefined ? params.length : 0.0;
        const options = { ...params, name };
        delete options.length;

        // Standardize parameter names mapping from Python to constructor args:
        // in Python it uses dipole_e1, dipole_e2
        // in Screen it might have Resolution, pixel_size
        // in correctors, it has angle
        if (clsName === "Marker") {
            return new Marker(name);
        } else if (clsName === "BPM" || clsName === "Screen" || clsName === "Aperture") {
            return new Cls(options);
        } else {
            return new Cls(length, options);
        }
    }

    function parseLattice(name) {
        if (!latticesData[name]) {
            return parseElement(name);
        }
        const childrenNames = latticesData[name];
        const children = childrenNames.map(childName => {
            if (latticesData[childName]) {
                return parseLattice(childName);
            } else {
                return parseElement(childName);
            }
        });
        return new Segment(children, name);
    }

    return parseLattice(rootName);
}
