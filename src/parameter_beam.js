// parameter_beam.js - Beam represented by statistical parameters

import { Beam } from "./beam.js";
import { Species } from "./species.js";
import { zeros } from "./math.js";

export class ParameterBeam extends Beam {
    constructor(mu, cov, energy, total_charge = 0.0, s = 0.0, species = null) {
        const activeSpecies = species || new Species("electron");
        super(energy, total_charge, s, activeSpecies);
        this.mu = mu; // Float64Array of size 7
        this.cov = cov; // Array of Float64Arrays of size 7x7
    }

    // Index properties
    get mu_x() { return this.mu[0]; }
    get sigma_x() { return Math.sqrt(Math.max(0.0, this.cov[0][0])); }
    get mu_px() { return this.mu[1]; }
    get sigma_px() { return Math.sqrt(Math.max(0.0, this.cov[1][1])); }
    get mu_y() { return this.mu[2]; }
    get sigma_y() { return Math.sqrt(Math.max(0.0, this.cov[2][2])); }
    get mu_py() { return this.mu[3]; }
    get sigma_py() { return Math.sqrt(Math.max(0.0, this.cov[3][3])); }
    get mu_tau() { return this.mu[4]; }
    get sigma_tau() { return Math.sqrt(Math.max(0.0, this.cov[4][4])); }
    get mu_p() { return this.mu[5]; }
    get sigma_p() { return Math.sqrt(Math.max(0.0, this.cov[5][5])); }

    get cov_xpx() { return this.cov[0][1]; }
    get cov_ypy() { return this.cov[2][3]; }
    get cov_taup() { return this.cov[4][5]; }
    get cov_xp() { return this.cov[0][5]; }
    get cov_pxp() { return this.cov[1][5]; }
    get cov_yp() { return this.cov[2][5]; }
    get cov_pyp() { return this.cov[3][5]; }

    clone() {
        const newMu = new Float64Array(7);
        newMu.set(this.mu);
        const newCov = zeros(7, 7);
        for (let i = 0; i < 7; i++) {
            newCov[i].set(this.cov[i]);
        }
        return new ParameterBeam(
            newMu,
            newCov,
            this.energy,
            this.total_charge,
            this.s,
            this.species.clone()
        );
    }

    static fromParameters({
        mu_x = 0.0, mu_px = 0.0, mu_y = 0.0, mu_py = 0.0, mu_tau = 0.0, mu_p = 0.0,
        sigma_x = 175e-6, sigma_px = 4e-6, sigma_y = 175e-6, sigma_py = 4e-6,
        sigma_tau = 8e-6, sigma_p = 2e-3,
        cov_xpx = 0.0, cov_ypy = 0.0, cov_taup = 0.0,
        cov_xp = 0.0, cov_pxp = 0.0, cov_yp = 0.0, cov_pyp = 0.0,
        energy = 1e8, total_charge = 0.0, s = 0.0, species = null
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

        return new ParameterBeam(mu, cov, energy, total_charge, s, species);
    }

    static fromTwiss({
        beta_x = 1.0, alpha_x = 0.0, emittance_x = 7.1971891e-13,
        beta_y = 1.0, alpha_y = 0.0, emittance_y = 7.1971891e-13,
        sigma_tau = 1e-6, sigma_p = 1e-6, cov_taup = 0.0,
        dispersion_x = 0.0, dispersion_px = 0.0, dispersion_y = 0.0, dispersion_py = 0.0,
        energy = 1e8, total_charge = 0.0, s = 0.0, species = null
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

        return ParameterBeam.fromParameters({
            sigma_x, sigma_px, sigma_y, sigma_py, sigma_tau, sigma_p,
            cov_xpx, cov_ypy, cov_taup, cov_xp, cov_pxp, cov_yp, cov_pyp,
            energy, total_charge, s, species
        });
    }
}
