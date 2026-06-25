// aperture.js - Physical aperture element

import { Element } from "./element.js";
import { ParticleBeam } from "./particle_beam.js";
import { identity } from "./math.js";

export class Aperture extends Element {
    constructor({
        x_max = Infinity,
        y_max = Infinity,
        shape = "rectangular", // "rectangular" or "elliptical"
        is_active = true,
        name = null
    } = {}) {
        super(name, 0.0); // Aperture length is 0.0
        this.x_max = x_max;
        this.y_max = y_max;
        this.shape = shape;
        this.is_active = is_active;
    }

    first_order_transfer_map(energy, species) {
        return identity(7);
    }

    track(incoming) {
        if (!this.is_active) {
            return incoming;
        }
        if (!(incoming instanceof ParticleBeam)) {
            console.warn("Aperture tracking is currently only supported for ParticleBeam.");
            return incoming;
        }

        const N = incoming.particles.length;
        const newProbs = new Float64Array(N);
        const x_lim = this.x_max;
        const y_lim = this.y_max;

        if (this.shape === "rectangular") {
            for (let i = 0; i < N; i++) {
                const px = incoming.particles[i][0];
                const py = incoming.particles[i][2];
                const inside = px > -x_lim && px < x_lim && py > -y_lim && py < y_lim;
                newProbs[i] = incoming.survival_probabilities[i] * (inside ? 1.0 : 0.0);
            }
        } else if (this.shape === "elliptical") {
            const x_lim_sq = x_lim * x_lim;
            const y_lim_sq = y_lim * y_lim;
            for (let i = 0; i < N; i++) {
                const px = incoming.particles[i][0];
                const py = incoming.particles[i][2];
                const val = (px * px) / x_lim_sq + (py * py) / y_lim_sq;
                const inside = val <= 1.0;
                newProbs[i] = incoming.survival_probabilities[i] * (inside ? 1.0 : 0.0);
            }
        } else {
            throw new Error(`Unsupported aperture shape '${this.shape}'`);
        }

        return new ParticleBeam(
            incoming.particles,
            incoming.energy,
            incoming.particle_charges,
            newProbs,
            incoming.s,
            incoming.species
        );
    }

    get is_skippable() {
        return !this.is_active;
    }
}
