// screen.js - Diagnostic screen element

import { Element } from "./element.js";
import { ParticleBeam } from "./particle_beam.js";
import { ParameterBeam } from "./parameter_beam.js";
import { identity } from "./math.js";

export class Screen extends Element {
    constructor({
        resolution = [256, 256],
        pixel_size = [1e-4, 1e-4],
        binning = 1,
        misalignment = [0.0, 0.0],
        is_blocking = false,
        is_active = false,
        name = null
    } = {}) {
        super(name, 0.0); // Screen length is 0.0
        this.resolution = resolution;
        this.pixel_size = pixel_size;
        this.binning = binning;
        this.misalignment = misalignment;
        this.is_blocking = is_blocking;
        this.is_active = is_active;
        this.read_beam = null;
    }

    first_order_transfer_map(energy, species) {
        return identity(7);
    }

    track(incoming) {
        if (this.is_active) {
            const copy_of_incoming = incoming.clone();
            const dx = this.misalignment[0];
            const dy = this.misalignment[1];
            
            if (copy_of_incoming instanceof ParameterBeam) {
                copy_of_incoming.mu[0] -= dx;
                copy_of_incoming.mu[2] -= dy;
            } else if (copy_of_incoming instanceof ParticleBeam) {
                for (let i = 0; i < copy_of_incoming.particles.length; i++) {
                    copy_of_incoming.particles[i][0] -= dx;
                    copy_of_incoming.particles[i][2] -= dy;
                }
            }
            this.read_beam = copy_of_incoming;
        }

        if (this.is_active && this.is_blocking) {
            const blocked = incoming.clone();
            if (blocked instanceof ParameterBeam) {
                blocked.total_charge = 0.0;
            } else if (blocked instanceof ParticleBeam) {
                blocked.survival_probabilities.fill(0.0);
                blocked.total_charge = 0.0;
            }
            return blocked;
        }

        return incoming.clone();
    }

    get_profile_image() {
        if (!this.is_active || !this.read_beam) return null;
        const beam = this.read_beam;
        const [w, h] = this.resolution;
        const [px, py] = this.pixel_size;
        const effective_w = Math.floor(w / this.binning);
        const effective_h = Math.floor(h / this.binning);
        
        const x_min = -w * px / 2;
        const x_max = w * px / 2;
        const y_min = -h * py / 2;
        const y_max = h * py / 2;

        const grid = new Float64Array(effective_w * effective_h);
        
        if (beam instanceof ParticleBeam) {
            for (let i = 0; i < beam.particles.length; i++) {
                const p = beam.particles[i];
                const prob = beam.survival_probabilities[i];
                if (prob <= 0.0) continue;
                
                const px_val = p[0];
                const py_val = p[2];
                
                if (px_val >= x_min && px_val < x_max && py_val >= y_min && py_val < y_max) {
                    const ix = Math.floor((px_val - x_min) / (x_max - x_min) * effective_w);
                    const iy = Math.floor((py_val - y_min) / (y_max - y_min) * effective_h);
                    const cl_ix = Math.min(effective_w - 1, Math.max(0, ix));
                    const cl_iy = Math.min(effective_h - 1, Math.max(0, iy));
                    grid[cl_iy * effective_w + cl_ix] += Math.abs(beam.particle_charges[i]) * prob;
                }
            }
        }
        return {
            width: effective_w,
            height: effective_h,
            data: grid,
            x_min, x_max, y_min, y_max
        };
    }

    get is_skippable() {
        return !this.is_active;
    }
}
