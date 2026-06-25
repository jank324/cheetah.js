// drift.js - Drift section element

import { Element } from "./element.js";
import { drift_matrix, base_ttensor } from "./track_methods.js";

export class Drift extends Element {
    constructor(length, { tracking_method = "linear", name = null } = {}) {
        super(name, length);
        this.tracking_method = tracking_method; // "linear" or "second_order"
    }

    first_order_transfer_map(energy, species) {
        return drift_matrix(this.length, energy, species);
    }

    second_order_transfer_map(energy, species) {
        const T = base_ttensor(this.length, 0.0, 0.0, 0.0, species, energy);
        const R = this.first_order_transfer_map(energy, species);
        for (let i = 0; i < 7; i++) {
            for (let j = 0; j < 7; j++) {
                T[i][6][j] = R[i][j];
            }
        }
        return T;
    }

    track(incoming) {
        if (this.tracking_method === "linear") {
            return this._trackFirstOrder(incoming);
        } else if (this.tracking_method === "second_order") {
            return this._trackSecondOrder(incoming);
        } else {
            throw new Error(`Unsupported tracking method '${this.tracking_method}' in Drift`);
        }
    }

    get is_skippable() {
        return this.tracking_method === "linear";
    }

    split(resolution) {
        if (this.length <= 0) return [this];
        const num_splits = Math.ceil(this.length / resolution);
        if (num_splits <= 1) return [this];
        const split_length = this.length / num_splits;
        const splits = [];
        for (let i = 0; i < num_splits; i++) {
            splits.push(new Drift(split_length, { tracking_method: this.tracking_method, name: `${this.name}_split_${i}` }));
        }
        return splits;
    }
}
