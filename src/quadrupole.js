// quadrupole.js - Quadrupole magnet element

import { Element, combinedRotationMisalignmentMatrix } from "./element.js";
import { base_rmatrix, base_ttensor } from "./track_methods.js";
import { zeros, matmul } from "./math.js";

export class Quadrupole extends Element {
    constructor(length, {
        k1 = 0.0,
        misalignment = [0.0, 0.0],
        tilt = 0.0,
        tracking_method = "linear",
        name = null
    } = {}) {
        super(name, length);
        this.k1 = k1;
        this.misalignment = misalignment;
        this.tilt = tilt;
        this.tracking_method = tracking_method; // "linear" or "second_order"
    }

    first_order_transfer_map(energy, species) {
        const R = base_rmatrix(this.length, this.k1, 0.0, species, energy);
        const [R_entry, R_exit] = combinedRotationMisalignmentMatrix(this.tilt, this.misalignment);
        return matmul(matmul(R_exit, R), R_entry);
    }

    second_order_transfer_map(energy, species) {
        const T = base_ttensor(this.length, this.k1, 0.0, 0.0, species, energy);
        const R = base_rmatrix(this.length, this.k1, 0.0, species, energy);
        for (let i = 0; i < 7; i++) {
            for (let j = 0; j < 7; j++) {
                T[i][6][j] = R[i][j];
            }
        }

        const [R_entry, R_exit] = combinedRotationMisalignmentMatrix(this.tilt, this.misalignment);
        
        // Tensor rotation: T_new = R_exit * T * R_entry * R_entry
        const T_new = new Array(7);
        for (let i = 0; i < 7; i++) T_new[i] = zeros(7, 7);

        for (let i = 0; i < 7; i++) {
            for (let n = 0; n < 7; n++) {
                for (let m = 0; m < 7; m++) {
                    let sum = 0.0;
                    for (let j = 0; j < 7; j++) {
                        const r_exit = R_exit[i][j];
                        if (r_exit === 0.0) continue;
                        for (let k = 0; k < 7; k++) {
                            const r_entry_kn = R_entry[k][n];
                            if (r_entry_kn === 0.0) continue;
                            for (let l = 0; l < 7; l++) {
                                const r_entry_lm = R_entry[l][m];
                                if (r_entry_lm === 0.0) continue;
                                sum += r_exit * T[j][k][l] * r_entry_kn * r_entry_lm;
                            }
                        }
                    }
                    T_new[i][n][m] = sum;
                }
            }
        }
        return T_new;
    }

    track(incoming) {
        if (this.tracking_method === "linear") {
            return this._trackFirstOrder(incoming);
        } else if (this.tracking_method === "second_order") {
            return this._trackSecondOrder(incoming);
        } else {
            throw new Error(`Unsupported tracking method '${this.tracking_method}' in Quadrupole`);
        }
    }

    get is_skippable() {
        return this.tracking_method === "linear";
    }

    get is_active() {
        return this.k1 !== 0.0;
    }

    split(resolution) {
        if (this.length <= 0) return [this];
        const num_splits = Math.ceil(this.length / resolution);
        if (num_splits <= 1) return [this];
        const split_length = this.length / num_splits;
        const splits = [];
        for (let i = 0; i < num_splits; i++) {
            splits.push(new Quadrupole(split_length, {
                k1: this.k1,
                misalignment: this.misalignment,
                tilt: this.tilt,
                tracking_method: this.tracking_method,
                name: `${this.name}_split_${i}`
            }));
        }
        return splits;
    }
}
