// solenoid.js - Solenoid element

import { Element, misalignmentMatrix } from "./element.js";
import { identity, matmul } from "./math.js";

export class Solenoid extends Element {
    constructor(length, { k = 0.0, misalignment = [0.0, 0.0], name = null } = {}) {
        super(name, length);
        this.k = k;
        this.misalignment = misalignment;
    }

    first_order_transfer_map(energy, species) {
        const gamma = energy / species.mass_eV;
        const lk = this.length * this.k;
        const c = Math.cos(lk);
        const s = Math.sin(lk);
        
        const s_k = (this.k === 0.0) ? this.length : Math.sin(lk) / this.k;
        const r56 = (gamma === 1.0) ? 0.0 : this.length / (1.0 - gamma * gamma);

        const R = identity(7);
        R[0][0] = c * c;
        R[0][1] = c * s_k;
        R[0][2] = s * c;
        R[0][3] = s * s_k;
        
        R[1][0] = -this.k * s * c;
        R[1][1] = c * c;
        R[1][2] = -this.k * s * s;
        R[1][3] = s * c;
        
        R[2][0] = -s * c;
        R[2][1] = -s * s_k;
        R[2][2] = c * c;
        R[2][3] = c * s_k;
        
        R[3][0] = this.k * s * s;
        R[3][1] = -s * c;
        R[3][2] = -this.k * s * c;
        R[3][3] = c * c;

        R[4][5] = r56;

        const [R_entry, R_exit] = misalignmentMatrix(this.misalignment);
        return matmul(matmul(R_exit, R), R_entry);
    }

    get is_skippable() {
        return true;
    }

    get is_active() {
        return this.k !== 0.0;
    }

    split(resolution) {
        if (this.length <= 0) return [this];
        const num_splits = Math.ceil(this.length / resolution);
        if (num_splits <= 1) return [this];
        const split_length = this.length / num_splits;
        const splits = [];
        for (let i = 0; i < num_splits; i++) {
            splits.push(new Solenoid(split_length, {
                k: this.k,
                misalignment: this.misalignment,
                name: `${this.name}_split_${i}`
            }));
        }
        return splits;
    }
}
