// bpm.js - Beam Position Monitor element

import { Element } from "./element.js";
import { identity } from "./math.js";

export class BPM extends Element {
    constructor({ is_active = false, misalignment = [0.0, 0.0], name = null } = {}) {
        super(name, 0.0); // BPM length is 0.0
        this.is_active = is_active;
        this.misalignment = misalignment;
        this.reading = [NaN, NaN];
    }

    first_order_transfer_map(energy, species) {
        return identity(7);
    }

    track(incoming) {
        if (this.is_active) {
            this.reading = [
                incoming.mu_x - this.misalignment[0],
                incoming.mu_y - this.misalignment[1]
            ];
        }
        return incoming.clone();
    }

    get is_skippable() {
        return !this.is_active;
    }
}
