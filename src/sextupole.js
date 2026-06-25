// sextupole.js - Sextupole magnet element

import { Element, combinedRotationMisalignmentMatrix } from "./element.js";
import { drift_matrix, base_ttensor } from "./track_methods.js";
import { applyTensorTransform } from "./dipole.js";

export class Sextupole extends Element {
    constructor(length, {
        k2 = 0.0,
        misalignment = [0.0, 0.0],
        tilt = 0.0,
        tracking_method = "second_order",
        name = null
    } = {}) {
        super(name, length);
        this.k2 = k2;
        this.misalignment = misalignment;
        this.tilt = tilt;
        this.tracking_method = tracking_method; // "linear" or "second_order"
    }

    first_order_transfer_map(energy, species) {
        return drift_matrix(this.length, energy, species);
    }

    second_order_transfer_map(energy, species) {
        const T = base_ttensor(this.length, 0.0, this.k2, 0.0, species, energy);
        const R = this.first_order_transfer_map(energy, species);
        for (let i = 0; i < 7; i++) {
            for (let j = 0; j < 7; j++) {
                T[i][6][j] = R[i][j];
            }
        }

        const [R_entry, R_exit] = combinedRotationMisalignmentMatrix(this.tilt, this.misalignment);
        return applyTensorTransform(T, R_exit, R_entry);
    }

    track(incoming) {
        if (this.tracking_method === "linear") {
            return this._trackFirstOrder(incoming);
        } else if (this.tracking_method === "second_order") {
            return this._trackSecondOrder(incoming);
        } else {
            throw new Error(`Unsupported tracking method '${this.tracking_method}' in Sextupole`);
        }
    }

    get is_skippable() {
        return false;
    }

    get is_active() {
        return this.k2 !== 0.0;
    }
}
