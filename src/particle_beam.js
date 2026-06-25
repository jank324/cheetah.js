// particle_beam.js - Particle beam class simulating individual macroparticles

import { Beam } from "./beam.js";
import { Species } from "./species.js";
import { 
    zeros, randn, weightedMean, weightedStd, weightedCovariance, 
    matchDistributionMoments 
} from "./math.js";

export class ParticleBeam extends Beam {
    constructor(particles, energy, particle_charges = null, survival_probabilities = null, s = 0.0, species = null) {
        const activeSpecies = species || new Species("electron");
        const N = particles.length;
        const charges = particle_charges || new Float64Array(N).fill(activeSpecies.charge_coulomb);
        const probs = survival_probabilities || new Float64Array(N).fill(1.0);
        
        // Compute total charge: sum of charges * survival probabilities
        let sumQ = 0.0;
        for (let i = 0; i < N; i++) {
            sumQ += charges[i] * probs[i];
        }
        
        super(energy, sumQ, s, activeSpecies);
        this.particles = particles; // Array of Float64Arrays [N][7]
        this.particle_charges = charges; // Float64Array [N]
        this.survival_probabilities = probs; // Float64Array [N]
    }

    // Coordinates helpers
    get x() {
        const out = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) out[i] = this.particles[i][0];
        return out;
    }
    get px() {
        const out = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) out[i] = this.particles[i][1];
        return out;
    }
    get y() {
        const out = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) out[i] = this.particles[i][2];
        return out;
    }
    get py() {
        const out = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) out[i] = this.particles[i][3];
        return out;
    }
    get tau() {
        const out = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) out[i] = this.particles[i][4];
        return out;
    }
    get p() {
        const out = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) out[i] = this.particles[i][5];
        return out;
    }

    // Statistical getters
    get mu_x() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][0];
        return weightedMean(vals, this.survival_probabilities);
    }
    get sigma_x() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][0];
        return weightedStd(vals, this.survival_probabilities, this.mu_x);
    }

    get mu_px() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][1];
        return weightedMean(vals, this.survival_probabilities);
    }
    get sigma_px() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][1];
        return weightedStd(vals, this.survival_probabilities, this.mu_px);
    }

    get mu_y() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][2];
        return weightedMean(vals, this.survival_probabilities);
    }
    get sigma_y() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][2];
        return weightedStd(vals, this.survival_probabilities, this.mu_y);
    }

    get mu_py() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][3];
        return weightedMean(vals, this.survival_probabilities);
    }
    get sigma_py() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][3];
        return weightedStd(vals, this.survival_probabilities, this.mu_py);
    }

    get mu_tau() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][4];
        return weightedMean(vals, this.survival_probabilities);
    }
    get sigma_tau() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][4];
        return weightedStd(vals, this.survival_probabilities, this.mu_tau);
    }

    get mu_p() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][5];
        return weightedMean(vals, this.survival_probabilities);
    }
    get sigma_p() {
        const vals = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) vals[i] = this.particles[i][5];
        return weightedStd(vals, this.survival_probabilities, this.mu_p);
    }

    // Covariances
    get cov_xpx() {
        const xs = new Float64Array(this.particles.length);
        const pxs = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) {
            xs[i] = this.particles[i][0];
            pxs[i] = this.particles[i][1];
        }
        return weightedCovariance(xs, pxs, this.survival_probabilities, this.mu_x, this.mu_px);
    }
    get cov_ypy() {
        const ys = new Float64Array(this.particles.length);
        const pys = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) {
            ys[i] = this.particles[i][2];
            pys[i] = this.particles[i][3];
        }
        return weightedCovariance(ys, pys, this.survival_probabilities, this.mu_y, this.mu_py);
    }
    get cov_taup() {
        const taus = new Float64Array(this.particles.length);
        const ps = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) {
            taus[i] = this.particles[i][4];
            ps[i] = this.particles[i][5];
        }
        return weightedCovariance(taus, ps, this.survival_probabilities, this.mu_tau, this.mu_p);
    }
    get cov_xp() {
        const xs = new Float64Array(this.particles.length);
        const ps = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) {
            xs[i] = this.particles[i][0];
            ps[i] = this.particles[i][5];
        }
        return weightedCovariance(xs, ps, this.survival_probabilities, this.mu_x, this.mu_p);
    }
    get cov_pxp() {
        const pxs = new Float64Array(this.particles.length);
        const ps = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) {
            pxs[i] = this.particles[i][1];
            ps[i] = this.particles[i][5];
        }
        return weightedCovariance(pxs, ps, this.survival_probabilities, this.mu_px, this.mu_p);
    }
    get cov_yp() {
        const ys = new Float64Array(this.particles.length);
        const ps = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) {
            ys[i] = this.particles[i][2];
            ps[i] = this.particles[i][5];
        }
        return weightedCovariance(ys, ps, this.survival_probabilities, this.mu_y, this.mu_p);
    }
    get cov_pyp() {
        const pys = new Float64Array(this.particles.length);
        const ps = new Float64Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) {
            pys[i] = this.particles[i][3];
            ps[i] = this.particles[i][5];
        }
        return weightedCovariance(pys, ps, this.survival_probabilities, this.mu_py, this.mu_p);
    }

    clone() {
        const newParticles = new Array(this.particles.length);
        for (let i = 0; i < this.particles.length; i++) {
            newParticles[i] = new Float64Array(7);
            newParticles[i].set(this.particles[i]);
        }
        const newCharges = new Float64Array(this.particle_charges.length);
        newCharges.set(this.particle_charges);
        const newProbs = new Float64Array(this.survival_probabilities.length);
        newProbs.set(this.survival_probabilities);
        return new ParticleBeam(
            newParticles,
            this.energy,
            newCharges,
            newProbs,
            this.s,
            this.species.clone()
        );
    }

    static fromDistribution(mu, cov, num_particles = 100000, energy = 1e8, total_charge = null, s = 0.0, species = null) {
        const activeSpecies = species || new Species("electron");
        const actualCharge = total_charge !== null ? total_charge : activeSpecies.charge_coulomb * num_particles;

        // Sample standard normal in 6D
        const samples = new Array(num_particles);
        for (let i = 0; i < num_particles; i++) {
            const p = new Float64Array(7);
            for (let c = 0; c < 6; c++) {
                p[c] = randn();
            }
            p[6] = 1.0;
            samples[i] = p;
        }

        // Target mu 6D
        const targetMu = new Float64Array(6);
        for (let c = 0; c < 6; c++) targetMu[c] = mu[c];

        // Match distribution moments exactly
        const matched = matchDistributionMoments(samples, targetMu, cov);

        const charges = new Float64Array(num_particles).fill(actualCharge / num_particles);
        const probs = new Float64Array(num_particles).fill(1.0);

        return new ParticleBeam(matched, energy, charges, probs, s, activeSpecies);
    }

    static fromParameters({
        num_particles = 100000,
        mu_x = 0.0, mu_px = 0.0, mu_y = 0.0, mu_py = 0.0, mu_tau = 0.0, mu_p = 0.0,
        sigma_x = 175e-6, sigma_px = 4e-6, sigma_y = 175e-6, sigma_py = 4e-6,
        sigma_tau = 8e-6, sigma_p = 2e-3,
        cov_xpx = 0.0, cov_ypy = 0.0, cov_taup = 0.0,
        cov_xp = 0.0, cov_pxp = 0.0, cov_yp = 0.0, cov_pyp = 0.0,
        energy = 1e8, total_charge = null, s = 0.0, species = null
    } = {}) {
        const mu = new Float64Array([mu_x, mu_px, mu_y, mu_py, mu_tau, mu_p, 1.0]);
        const cov = zeros(7, 7);

        cov[0][0] = sigma_x * sigma_x;
        cov[0][1] = cov_xpx;
        cov[1][0] = cov_xpx;
        cov[1][1] = sigma_px * sigma_px;
        
        cov[2][2] = sigma_y * sigma_y;
        cov[2][3] = cov_ypy;
        cov[3][2] = cov_ypy;
        cov[3][3] = sigma_py * sigma_py;
        
        cov[4][4] = sigma_tau * sigma_tau;
        cov[4][5] = cov_taup;
        cov[5][4] = cov_taup;
        cov[5][5] = sigma_p * sigma_p;
        
        cov[0][5] = cov_xp;
        cov[5][0] = cov_xp;
        cov[1][5] = cov_pxp;
        cov[5][1] = cov_pxp;
        cov[2][5] = cov_yp;
        cov[5][2] = cov_yp;
        cov[3][5] = cov_pyp;
        cov[5][3] = cov_pyp;

        return ParticleBeam.fromDistribution(mu, cov, num_particles, energy, total_charge, s, species);
    }

    static fromTwiss({
        num_particles = 100000,
        beta_x = 1.0, alpha_x = 0.0, emittance_x = 7.1971891e-13,
        beta_y = 1.0, alpha_y = 0.0, emittance_y = 7.1971891e-13,
        sigma_tau = 1e-6, sigma_p = 1e-6, cov_taup = 0.0,
        dispersion_x = 0.0, dispersion_px = 0.0, dispersion_y = 0.0, dispersion_py = 0.0,
        energy = 1e8, total_charge = null, s = 0.0, species = null
    } = {}) {
        if (beta_x <= 0) throw new Error("beta_x must be greater than 0");
        if (beta_y <= 0) throw new Error("beta_y must be greater than 0");

        const sigma_x = Math.sqrt(emittance_x * beta_x + dispersion_x * dispersion_x * sigma_p * sigma_p);
        const sigma_px = Math.sqrt(emittance_x * (1 + alpha_x * alpha_x) / beta_x + dispersion_px * dispersion_px * sigma_p * sigma_p);
        const sigma_y = Math.sqrt(emittance_y * beta_y + dispersion_y * dispersion_y * sigma_p * sigma_p);
        const sigma_py = Math.sqrt(emittance_y * (1 + alpha_y * alpha_y) / beta_y + dispersion_py * dispersion_py * sigma_p * sigma_p);
        
        const cov_xpx = -emittance_x * alpha_x + dispersion_x * dispersion_px * sigma_p * sigma_p;
        const cov_ypy = -emittance_y * alpha_y + dispersion_y * dispersion_py * sigma_p * sigma_p;
        
        const cov_xp = dispersion_x * sigma_p * sigma_p;
        const cov_pxp = dispersion_px * sigma_p * sigma_p;
        const cov_yp = dispersion_y * sigma_p * sigma_p;
        const cov_pyp = dispersion_py * sigma_p * sigma_p;

        return ParticleBeam.fromParameters({
            num_particles,
            sigma_x, sigma_px, sigma_y, sigma_py, sigma_tau, sigma_p,
            cov_xpx, cov_ypy, cov_taup, cov_xp, cov_pxp, cov_yp, cov_pyp,
            energy, total_charge, s, species
        });
    }
}
