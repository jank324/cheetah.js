// marker.js - Marker element

import { Element } from "./element.js";
import { identity } from "./math.js";

export class Marker extends Element {
    constructor(name = null) {
        super(name, 0.0);
    }

    first_order_transfer_map(energy, species) {
        return identity(7);
    }

    track(incoming) {
        return incoming.clone();
    }

    get is_skippable() {
        return true;
    }
}
