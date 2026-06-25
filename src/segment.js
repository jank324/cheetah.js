// segment.js - Segment class representing a sequence of elements

import { Element } from "./element.js";
import { CustomTransferMap } from "./custom_transfer_map.js";
import { identity, matmul } from "./math.js";

export class Segment extends Element {
    constructor(elements, name = null) {
        super(name, 0.0);
        this.elements = elements;

        // Make elements accessible by their name as properties
        for (const element of elements) {
            if (element.name) {
                if (this[element.name]) {
                    if (Array.isArray(this[element.name])) {
                        this[element.name].push(element);
                    } else {
                        this[element.name] = [this[element.name], element];
                    }
                } else {
                    this[element.name] = element;
                }
            }
        }
    }

    get element_names() {
        return this.elements.map(e => e.name);
    }

    get length() {
        let total = 0.0;
        for (const element of this.elements) {
            total += element.length;
        }
        return total;
    }

    get is_skippable() {
        for (const element of this.elements) {
            if (!element.is_skippable) return false;
        }
        return true;
    }

    element_index(element_name) {
        const idx = this.element_names.indexOf(element_name);
        if (idx === -1) {
            throw new Error(`Element '${element_name}' not found in segment.`);
        }
        return idx;
    }

    subcell(start = null, end = null, include_start = true, include_end = true) {
        let is_in_subcell = start === null;
        const subelements = [];

        for (const element of this.elements) {
            if (element.name === start) {
                is_in_subcell = true;
                if (include_start) {
                    subelements.push(element);
                }
                continue;
            }

            if (element.name === end) {
                if (include_end && is_in_subcell) {
                    subelements.push(element);
                }
                break;
            }

            if (is_in_subcell) {
                subelements.push(element);
            }
        }
        return new Segment(subelements, `${this.name}_subcell`);
    }

    flattened() {
        const flattened_elements = [];
        for (const element of this.elements) {
            if (element instanceof Segment) {
                flattened_elements.push(...element.flattened().elements);
            } else {
                flattened_elements.push(element);
            }
        }
        return new Segment(flattened_elements, this.name);
    }

    reversed() {
        const reversed_elements = [...this.elements].reverse().map(element => {
            if (element instanceof Segment) {
                return element.reversed();
            }
            return element;
        });
        return new Segment(reversed_elements, `${this.name}_reversed`);
    }

    first_order_transfer_map(energy, species) {
        if (!this.is_skippable) {
            return null;
        }
        let tm = identity(7);
        for (const element of this.elements) {
            tm = matmul(element.first_order_transfer_map(energy, species), tm);
        }
        return tm;
    }

    transfer_maps_merged(incoming_beam, except_for = []) {
        const merged_elements = [];
        let skippable_elements = [];
        let tracked_beam = incoming_beam;

        for (const element of this.elements) {
            if (element.is_skippable && !except_for.includes(element.name)) {
                skippable_elements.push(element);
            } else {
                if (skippable_elements.length === 1) {
                    merged_elements.push(skippable_elements[0]);
                    tracked_beam = skippable_elements[0].track(tracked_beam);
                } else if (skippable_elements.length > 1) {
                    merged_elements.push(
                        CustomTransferMap.from_merging_elements(skippable_elements, tracked_beam)
                    );
                    tracked_beam = merged_elements[merged_elements.length - 1].track(tracked_beam);
                }
                skippable_elements = [];
                merged_elements.push(element);
                tracked_beam = element.track(tracked_beam);
            }
        }

        if (skippable_elements.length === 1) {
            merged_elements.push(skippable_elements[0]);
        } else if (skippable_elements.length > 1) {
            merged_elements.push(
                CustomTransferMap.from_merging_elements(skippable_elements, tracked_beam)
            );
        }

        return new Segment(merged_elements, this.name);
    }

    track(incoming) {
        if (this.is_skippable) {
            return this._trackFirstOrder(incoming);
        }

        const todos = [];
        let continuous_skippable_elements = [];
        for (const element of this.elements) {
            if (element.is_skippable) {
                continuous_skippable_elements.push(element);
            } else {
                if (continuous_skippable_elements.length > 0) {
                    todos.push(new Segment(continuous_skippable_elements));
                    continuous_skippable_elements = [];
                }
                todos.push(element);
            }
        }
        if (continuous_skippable_elements.length > 0) {
            todos.push(new Segment(continuous_skippable_elements));
        }

        let current_beam = incoming;
        for (const todo of todos) {
            current_beam = todo.track(current_beam);
        }
        return current_beam;
    }

    split(resolution) {
        const split_elements = [];
        for (const element of this.elements) {
            split_elements.push(...element.split(resolution));
        }
        return new Segment(split_elements, this.name);
    }

    get_beams_along_segment(incoming, resolution = null) {
        if (resolution !== null) {
            const split_segment = this.split(resolution);
            return split_segment.get_beams_along_segment(incoming);
        }

        const beams = [incoming];
        let current_beam = incoming;
        for (const element of this.elements) {
            current_beam = element.track(current_beam);
            beams.push(current_beam);
        }
        return beams;
    }

    get_beam_attrs_along_segment(attr_names, incoming, resolution = null) {
        const beams = this.get_beams_along_segment(incoming, resolution);
        const singleAttr = typeof attr_names === "string";
        const attrs = singleAttr ? [attr_names] : attr_names;
        
        const results = {};
        for (const attr of attrs) {
            results[attr] = new Float64Array(beams.length);
        }

        for (let i = 0; i < beams.length; i++) {
            const beam = beams[i];
            for (const attr of attrs) {
                results[attr][i] = beam[attr];
            }
        }

        return singleAttr ? results[attrs[0]] : results;
    }

    clone() {
        return new Segment(this.elements.map(e => e.clone()), this.name);
    }
}
