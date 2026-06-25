// cavity.js - RF Accelerating Cavity element

import { Element } from "./element.js";
import { computeRelativisticFactors } from "./track_methods.js";
import { ParameterBeam } from "./parameter_beam.js";
import { ParticleBeam } from "./particle_beam.js";
import { identity, zeros, log1pdiv, matvec, matmul, transpose } from "./math.js";

const speed_of_light = 299792458.0;

export class Cavity extends Element {
    constructor(length, {
        voltage = 0.0,
        phase = 0.0, // phase in degrees
        frequency = 0.0, // frequency in Hz
        cavity_type = "standing_wave", // "standing_wave" or "traveling_wave"
        name = null
    } = {}) {
        super(name, length);
        this.voltage = voltage;
        this.phase = phase;
        this.frequency = frequency;
        this.cavity_type = cavity_type;
    }

    first_order_transfer_map(energy, species) {
        if (!this.is_active) {
            return identity(7);
        }
        return this._cavity_rmatrix(energy, species);
    }

    _cavity_rmatrix(energy, species) {
        const phi = this.phase * Math.PI / 180.0;
        const effective_voltage = -this.voltage * species.num_elementary_charges;
        const delta_energy = effective_voltage * Math.cos(phi);

        const Ei = energy / species.mass_eV;
        const dE = delta_energy / species.mass_eV;
        const Ef = Ei + dE;
        const Ep = dE / this.length;

        const k = 2.0 * Math.PI * this.frequency / speed_of_light;

        let r11 = 1.0, r12 = this.length, r21 = 0.0, r22 = 1.0;
        let r55 = 1.0, r56 = 0.0, r65 = 0.0, r66 = 1.0;

        if (this.cavity_type === "standing_wave") {
            const alpha = Math.sqrt(0.125) * effective_voltage / energy * log1pdiv(delta_energy / energy);
            const beta0 = Math.sqrt(Math.max(0.0, 1.0 - 1.0 / (Ei * Ei)));
            const beta1 = Math.sqrt(Math.max(0.0, 1.0 - 1.0 / (Ef * Ef)));

            r11 = Math.cos(alpha) - Math.sqrt(2.0) * Math.cos(phi) * Math.sin(alpha);
            
            const sinc_val = (alpha === 0.0) ? 1.0 : Math.sin(alpha) / alpha;
            r12 = sinc_val * log1pdiv(delta_energy / energy) * this.length;
            
            r21 = -effective_voltage / ((energy + delta_energy) * Math.sqrt(2.0) * this.length) * (0.5 + Math.cos(phi) * Math.cos(phi)) * Math.sin(alpha);
            r22 = Ei / Ef * (Math.cos(alpha) + Math.sqrt(2.0) * Math.cos(phi) * Math.sin(alpha));

            r55 = 1.0 + (dE !== 0.0 ? k * this.length * beta0 * Math.tan(phi) * (Ei * Ef * (beta0 * beta1 - 1.0) + 1.0) / (beta1 * Ef * dE) : 0.0);
            r56 = -this.length / (Ef * Ef * Ei * beta1) * (Ef + Ei) / (beta1 + beta0);
            r65 = k * Math.sin(phi) * effective_voltage / (beta1 * (energy + delta_energy));
            r66 = Ei / Ef * beta0 / beta1;
        } else if (this.cavity_type === "traveling_wave") {
            // M_body
            const m_body_00 = 1.0;
            const m_body_01 = this.length * log1pdiv(dE / Ei);
            const m_body_10 = 0.0;
            const m_body_11 = Ei / Ef;

            // M_f_entry
            const m_entry_00 = 1.0;
            const m_entry_01 = 0.0;
            const m_entry_10 = -Ep / (2.0 * Ei);
            const m_entry_11 = 1.0;

            // M_f_exit
            const m_exit_00 = 1.0;
            const m_exit_01 = 0.0;
            const m_exit_10 = Ep / (2.0 * Ef);
            const m_exit_11 = 1.0;

            // M_combined = M_f_exit * M_body * M_f_entry
            const t_00 = m_body_00 * m_entry_00 + m_body_01 * m_entry_10;
            const t_01 = m_body_00 * m_entry_01 + m_body_01 * m_entry_11;
            const t_10 = m_body_10 * m_entry_00 + m_body_11 * m_entry_10;
            const t_11 = m_body_10 * m_entry_01 + m_body_11 * m_entry_11;

            r11 = m_exit_00 * t_00 + m_exit_01 * t_10;
            r12 = m_exit_00 * t_01 + m_exit_01 * t_11;
            r21 = m_exit_10 * t_00 + m_exit_11 * t_10;
            r22 = m_exit_10 * t_01 + m_exit_11 * t_11;

            r55 = 1.0;
            r56 = 0.0;
            r65 = k * Math.sin(phi) * effective_voltage / (energy + delta_energy);
            r66 = r22;
        }

        const R = identity(7);
        R[0][0] = r11;
        R[0][1] = r12;
        R[1][0] = r21;
        R[1][1] = r22;
        R[2][2] = r11;
        R[2][3] = r12;
        R[3][2] = r21;
        R[3][3] = r22;
        R[4][4] = r55;
        R[4][5] = r56;
        R[5][4] = r65;
        R[5][5] = r66;
        return R;
    }

    track(incoming) {
        const { gamma: gamma0, igamma2, beta: beta0 } = computeRelativisticFactors(incoming.energy, incoming.species.mass_eV);
        const phi = this.phase * Math.PI / 180.0;
        const tm = this.first_order_transfer_map(incoming.energy, incoming.species);

        const delta_energy = this.voltage * Math.cos(phi) * incoming.species.num_elementary_charges * -1;
        const outgoing_energy = incoming.energy + delta_energy;
        const { gamma: gamma1, beta: beta1 } = computeRelativisticFactors(outgoing_energy, incoming.species.mass_eV);

        let T566 = 1.5 * this.length * igamma2 / (beta0 * beta0 * beta0);
        let T556 = 0.0;
        let T555 = 0.0;

        const k = 2.0 * Math.PI * this.frequency / speed_of_light;
        const dgamma = this.voltage / incoming.species.mass_eV;

        if (delta_energy > 0.0) {
            const b0_3 = beta0 * beta0 * beta0;
            const b1_3 = beta1 * beta1 * beta1;
            const g0_3 = gamma0 * gamma0 * gamma0;
            const g1_3 = gamma1 * gamma1 * gamma1;
            T566 = this.length * (b0_3 * g0_3 - b1_3 * g1_3) / (2.0 * beta0 * b1_3 * gamma0 * (gamma0 - gamma1) * g1_3);
            T556 = beta0 * k * this.length * dgamma * gamma0 * (b1_3 * g1_3 + beta0 * (gamma0 - g1_3)) * Math.sin(phi) / (b1_3 * g1_3 * (gamma0 - gamma1) * (gamma0 - gamma1));
            T555 = (beta0 * beta0) * (k * k) * this.length * dgamma / 2.0 * (
                dgamma * (2.0 * gamma0 * g1_3 * (beta0 * b1_3 - 1.0) + gamma0 * gamma0 + 3.0 * gamma1 * gamma1 - 2.0) / (b1_3 * g1_3 * Math.pow(gamma0 - gamma1, 3)) * Math.sin(phi) * Math.sin(phi)
                - (gamma1 * gamma0 * (beta1 * beta0 - 1.0) + 1.0) / (beta1 * gamma1 * (gamma0 - gamma1) * (gamma0 - gamma1)) * Math.cos(phi)
            );
        }

        const newS = incoming.s + this.length;

        if (incoming instanceof ParameterBeam) {
            const outgoing_mu = matvec(tm, incoming.mu);
            const outgoing_cov = matmul(matmul(tm, incoming.cov), transpose(tm));

            outgoing_mu[5] = incoming.mu[5] * incoming.energy * beta0 / (outgoing_energy * beta1) 
                + this.voltage * beta0 / (outgoing_energy * beta1) * (Math.cos(-incoming.mu[4] * beta0 * k + phi) - Math.cos(phi));
            outgoing_cov[5][5] = incoming.cov[5][5];

            outgoing_mu[4] = outgoing_mu[4] + T566 * incoming.mu[5] * incoming.mu[5] 
                + T556 * incoming.mu[4] * incoming.mu[5] + T555 * incoming.mu[4] * incoming.mu[4];
            
            const cov44 = T566 * incoming.cov[5][5] * incoming.cov[5][5] 
                + T556 * incoming.cov[4][5] * incoming.cov[5][5] + T555 * incoming.cov[4][4] * incoming.cov[4][4];
            outgoing_cov[4][4] = cov44;
            outgoing_cov[4][5] = cov44;
            outgoing_cov[5][4] = cov44;

            return new ParameterBeam(
                outgoing_mu,
                outgoing_cov,
                outgoing_energy,
                incoming.total_charge,
                newS,
                incoming.species
            );
        } else {
            const N = incoming.particles.length;
            const outgoing_particles = new Array(N);
            for (let i = 0; i < N; i++) {
                const p = incoming.particles[i];
                const p_out = matvec(tm, p);

                p_out[5] = p[5] * incoming.energy * beta0 / (outgoing_energy * beta1) 
                    + this.voltage * beta0 / (outgoing_energy * beta1) * (Math.cos(-p[4] * beta0 * k + phi) - Math.cos(phi));
                
                p_out[4] = p_out[4] + T566 * p[5] * p[5] + T556 * p[4] * p[5] + T555 * p[4] * p[4];
                outgoing_particles[i] = p_out;
            }
            return new ParticleBeam(
                outgoing_particles,
                outgoing_energy,
                incoming.particle_charges,
                incoming.survival_probabilities,
                newS,
                incoming.species
            );
        }
    }

    get is_skippable() {
        return !this.is_active;
    }

    get is_active() {
        return this.voltage !== 0.0;
    }
}
