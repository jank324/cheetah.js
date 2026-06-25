// custom_transfer_map.js - Merged or custom predefined transfer map element

import { Element } from "./element.js";
import { identity, matmul } from "./math.js";

export class CustomTransferMap extends Element {
    constructor(predefined_transfer_map, { length = 0.0, name = null } = {}) {
        super(name, length);
        this.predefined_transfer_map = predefined_transfer_map; // 7x7 matrix
    }

    static from_merging_elements(elements, incoming_beam) {
        let tm = identity(7);
        let current_beam = incoming_beam;
        
        for (const element of elements) {
            const elm_tm = element.first_order_transfer_map(current_beam.energy, current_beam.species);
            tm = matmul(elm_tm, tm);
            current_beam = element.track(current_beam);
        }

        let combined_length = 0.0;
        for (const element of elements) {
            combined_length += element.length;
        }

        const combined_name = "combined_" + elements.map(e => e.name).join("_");
        return new CustomTransferMap(tm, { length: combined_length, name: combined_name });
    }

    first_order_transfer_map(energy, species) {
        return this.predefined_transfer_map;
    }

    get is_skippable() {
        return true;
    }
}
