// dipole.js - Sector bending dipole magnet element

import { Element, rotationMatrix } from "./element.js";
import { base_rmatrix, base_ttensor } from "./track_methods.js";
import { zeros, identity, matmul, transpose } from "./math.js";

// Helper to transform T_ijk using M_exit (exit matrix) and M_entry (entry matrix)
export function applyTensorTransform(T, M_exit, M_entry) {
    const T_new = new Array(7);
    for (let i = 0; i < 7; i++) T_new[i] = zeros(7, 7);

    for (let i = 0; i < 7; i++) {
        for (let n = 0; n < 7; n++) {
            for (let m = 0; m < 7; m++) {
                let sum = 0.0;
                for (let j = 0; j < 7; j++) {
                    const m_exit = M_exit[i][j];
                    if (m_exit === 0.0) continue;
                    for (let k = 0; k < 7; k++) {
                        const m_entry_kn = M_entry[k][n];
                        if (m_entry_kn === 0.0) continue;
                        for (let l = 0; l < 7; l++) {
                            const m_entry_lm = M_entry[l][m];
                            if (m_entry_lm === 0.0) continue;
                            sum += m_exit * T[j][k][l] * m_entry_kn * m_entry_lm;
                        }
                    }
                }
                T_new[i][n][m] = sum;
            }
        }
    }
    return T_new;
}

export class Dipole extends Element {
    constructor(length, {
        angle = 0.0,
        k1 = 0.0,
        dipole_e1 = 0.0,
        dipole_e2 = 0.0,
        tilt = 0.0,
        gap = 0.0,
        gap_exit = null,
        fringe_integral = 0.0,
        fringe_integral_exit = null,
        tracking_method = "linear",
        name = null
    } = {}) {
        super(name, length);
        this.angle = angle;
        this.k1 = k1;
        this.dipole_e1 = dipole_e1;
        this.dipole_e2 = dipole_e2;
        this.tilt = tilt;
        this.gap = gap;
        this.gap_exit = gap_exit !== null ? gap_exit : gap;
        this.fringe_integral = fringe_integral;
        this.fringe_integral_exit = fringe_integral_exit !== null ? fringe_integral_exit : fringe_integral;
        this.tracking_method = tracking_method; // "linear" or "second_order"
    }

    get hx() {
        return this.length > 0 ? this.angle / this.length : 0.0;
    }

    _transfer_map_enter() {
        const tm = identity(7);
        const e1 = this.dipole_e1;
        const cos_e1 = Math.cos(e1);
        if (cos_e1 === 0.0) return tm;
        const sec_e = 1.0 / cos_e1;
        const phi = this.fringe_integral * this.hx * this.gap * sec_e * (1.0 + Math.sin(e1) * Math.sin(e1));
        tm[1][0] = this.hx * Math.tan(e1);
        tm[3][2] = -this.hx * Math.tan(e1 - phi);
        return tm;
    }

    _transfer_map_exit() {
        const tm = identity(7);
        const e2 = this.dipole_e2;
        const cos_e2 = Math.cos(e2);
        if (cos_e2 === 0.0) return tm;
        const sec_e = 1.0 / cos_e2;
        const phi = this.fringe_integral_exit * this.hx * this.gap_exit * sec_e * (1.0 + Math.sin(e2) * Math.sin(e2));
        tm[1][0] = this.hx * Math.tan(e2);
        tm[3][2] = -this.hx * Math.tan(e2 - phi);
        return tm;
    }

    first_order_transfer_map(energy, species) {
        const R_enter = this._transfer_map_enter();
        const R_exit = this._transfer_map_exit();
        const R = base_rmatrix(this.length, this.k1, this.hx, species, energy);
        
        // Fringe fields
        let R_total = matmul(matmul(R_exit, R), R_enter);

        // Rotation
        if (this.tilt !== 0.0) {
            const rot = rotationMatrix(this.tilt);
            R_total = matmul(matmul(transpose(rot), R_total), rot);
        }
        return R_total;
    }

    second_order_transfer_map(energy, species) {
        const R_enter = this._transfer_map_enter();
        const R_exit = this._transfer_map_exit();
        const T = base_ttensor(this.length, this.k1, 0.0, this.hx, species, energy);
        const R = base_rmatrix(this.length, this.k1, this.hx, species, energy);

        for (let i = 0; i < 7; i++) {
            for (let j = 0; j < 7; j++) {
                T[i][6][j] = R[i][j];
            }
        }

        // Apply fringe fields to T
        let T_total = applyTensorTransform(T, R_exit, R_enter);

        // Apply rotation
        if (this.tilt !== 0.0) {
            const rot = rotationMatrix(this.tilt);
            T_total = applyTensorTransform(T_total, transpose(rot), rot);
        }

        return T_total;
    }

    track(incoming) {
        if (this.tracking_method === "linear") {
            return this._trackFirstOrder(incoming);
        } else if (this.tracking_method === "second_order") {
            return this._trackSecondOrder(incoming);
        } else {
            throw new Error(`Unsupported tracking method '${this.tracking_method}' in Dipole`);
        }
    }

    get is_skippable() {
        return this.tracking_method === "linear";
    }

    get is_active() {
        return this.angle !== 0.0;
    }
}
