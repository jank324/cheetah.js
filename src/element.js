// element.js - Base class for lattice elements

import { ParameterBeam } from "./parameter_beam.js";
import { ParticleBeam } from "./particle_beam.js";
import { matvec, matmul, transpose, zeros, identity } from "./math.js";

export class Element {
    constructor(name = null, length = 0.0) {
        this.name = name || `element_${Math.random().toString(36).substr(2, 9)}`;
        this._length = length;
        this._is_active = true;
    }

    get length() {
        return this._length;
    }

    set length(val) {
        this._length = val;
    }

    first_order_transfer_map(energy, species) {
        throw new Error("first_order_transfer_map not implemented");
    }

    second_order_transfer_map(energy, species) {
        throw new Error("second_order_transfer_map not implemented");
    }

    track(incoming) {
        return this._trackFirstOrder(incoming);
    }

    _trackFirstOrder(incoming) {
        const tm = this.first_order_transfer_map(incoming.energy, incoming.species);
        const newS = incoming.s + this.length;
        if (incoming instanceof ParameterBeam) {
            const newMu = matvec(tm, incoming.mu);
            const newCov = matmul(matmul(tm, incoming.cov), transpose(tm));
            return new ParameterBeam(
                newMu,
                newCov,
                incoming.energy,
                incoming.total_charge,
                newS,
                incoming.species
            );
        } else {
            const N = incoming.particles.length;
            const newParticles = new Array(N);
            for (let i = 0; i < N; i++) {
                newParticles[i] = matvec(tm, incoming.particles[i]);
            }
            return new ParticleBeam(
                newParticles,
                incoming.energy,
                incoming.particle_charges,
                incoming.survival_probabilities,
                newS,
                incoming.species
            );
        }
    }

    _trackSecondOrder(incoming) {
        if (!(incoming instanceof ParticleBeam)) {
            throw new Error("Second-order tracking is currently only supported for ParticleBeam.");
        }
        const N = incoming.particles.length;
        const T = this.second_order_transfer_map(incoming.energy, incoming.species);
        const newParticles = new Array(N);
        for (let i = 0; i < N; i++) {
            const p = incoming.particles[i];
            const p_out = new Float64Array(7);
            for (let coord = 0; coord < 7; coord++) {
                let sum = 0.0;
                for (let j = 0; j < 7; j++) {
                    const pj_val = p[j];
                    if (pj_val === 0.0) continue;
                    for (let k = 0; k < 7; k++) {
                        const pk_val = p[k];
                        if (pk_val === 0.0) continue;
                        sum += T[coord][j][k] * pj_val * pk_val;
                    }
                }
                p_out[coord] = sum;
            }
            newParticles[i] = p_out;
        }
        return new ParticleBeam(
            newParticles,
            incoming.energy,
            incoming.particle_charges,
            incoming.survival_probabilities,
            incoming.s + this.length,
            incoming.species
        );
    }

    get is_skippable() {
        return false;
    }

    get is_active() {
        return this._is_active;
    }

    set is_active(val) {
        this._is_active = val;
    }

    split(resolution) {
        return [this];
    }
}

// ----------------------------------------------------
// Rotation / Misalignment Matrices helpers
// ----------------------------------------------------

export function rotationMatrix(angle) {
    const tm = identity(7);
    if (angle === 0.0) return tm;
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    tm[0][0] = c;
    tm[0][2] = s;
    tm[1][1] = c;
    tm[1][3] = s;
    tm[2][0] = -s;
    tm[2][2] = c;
    tm[3][1] = -s;
    tm[3][3] = c;
    return tm;
}

export function misalignmentMatrix(misalignment) {
    const dx = misalignment[0];
    const dy = misalignment[1];
    const R_entry = identity(7);
    const R_exit = identity(7);
    R_entry[0][6] = -dx;
    R_entry[2][6] = -dy;
    R_exit[0][6] = dx;
    R_exit[2][6] = dy;
    return [R_entry, R_exit];
}

export function combinedRotationMisalignmentMatrix(angle, misalignment) {
    const c = Math.cos(angle);
    const s = Math.sin(angle);
    const dx = misalignment[0];
    const dy = misalignment[1];

    const R_entry = identity(7);
    R_entry[0][0] = c;
    R_entry[0][2] = s;
    R_entry[1][1] = c;
    R_entry[1][3] = s;
    R_entry[2][0] = -s;
    R_entry[2][2] = c;
    R_entry[3][1] = -s;
    R_entry[3][3] = c;

    R_entry[0][6] = -dx * c - dy * s;
    R_entry[2][6] = dx * s - dy * c;

    const R_exit = identity(7);
    R_exit[0][0] = c;
    R_exit[0][2] = -s;
    R_exit[1][1] = c;
    R_exit[1][3] = -s;
    R_exit[2][0] = s;
    R_exit[2][2] = c;
    R_exit[3][1] = s;
    R_exit[3][3] = c;

    R_exit[0][6] = dx;
    R_exit[2][6] = dy;

    return [R_entry, R_exit];
}
