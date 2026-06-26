// horizontal_corrector.js - Horizontal corrector magnet

import { Element } from "./element.js";
import { computeRelativisticFactors } from "./track_methods.js";
import { identity } from "./math.js";

export class HorizontalCorrector extends Element {
    constructor(length, { angle = 0.0, name = null } = {}) {
        super(name, length);
        this.angle = angle;
    }

    first_order_transfer_map(energy, species) {
        const { igamma2, beta } = computeRelativisticFactors(energy, species.mass_eV);
        const tm = identity(7);
        tm[0][1] = this.length;
        tm[1][6] = this.angle;
        tm[2][3] = this.length;
        if (beta > 0.0) {
            tm[4][5] = -this.length / (beta * beta) * igamma2;
        }
        return tm;
    }

    get is_skippable() {
        return true;
    }

    get is_active() {
        return this.angle !== 0.0;
    }
}
