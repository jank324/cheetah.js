// beam.js - Base class for beams

export class Beam {
    constructor(energy, total_charge, s, species) {
        this.energy = energy; // reference energy in eV
        this.total_charge = total_charge; // total charge in Coulombs
        this.s = s; // position along beamline in meters
        this.species = species; // Species object
    }

    get relativistic_gamma() {
        return this.energy / this.species.mass_eV;
    }

    get relativistic_beta() {
        const gamma = this.relativistic_gamma;
        if (gamma <= 1.0) return 0.0;
        return Math.sqrt(1.0 - 1.0 / (gamma * gamma));
    }

    get p0c() {
        return this.relativistic_beta * this.relativistic_gamma * this.species.mass_eV;
    }

    // Abstract properties to be implemented by subclasses
    get mu_x() { throw new Error("Not implemented"); }
    get sigma_x() { throw new Error("Not implemented"); }
    get mu_px() { throw new Error("Not implemented"); }
    get sigma_px() { throw new Error("Not implemented"); }
    get mu_y() { throw new Error("Not implemented"); }
    get sigma_y() { throw new Error("Not implemented"); }
    get mu_py() { throw new Error("Not implemented"); }
    get sigma_py() { throw new Error("Not implemented"); }
    get mu_tau() { throw new Error("Not implemented"); }
    get sigma_tau() { throw new Error("Not implemented"); }
    get mu_p() { throw new Error("Not implemented"); }
    get sigma_p() { throw new Error("Not implemented"); }

    get cov_xpx() { throw new Error("Not implemented"); }
    get cov_ypy() { throw new Error("Not implemented"); }
    get cov_taup() { throw new Error("Not implemented"); }
    get cov_xp() { throw new Error("Not implemented"); }
    get cov_pxp() { throw new Error("Not implemented"); }
    get cov_yp() { throw new Error("Not implemented"); }
    get cov_pyp() { throw new Error("Not implemented"); }

    get projected_emittance_x() {
        const sx = this.sigma_x;
        const spx = this.sigma_px;
        const cxpx = this.cov_xpx;
        const val = sx * sx * spx * spx - cxpx * cxpx;
        return Math.sqrt(val > 0 ? val : 0.0);
    }

    get emittance_x() {
        const sx2 = this.sigma_x * this.sigma_x;
        const spx2 = this.sigma_px * this.sigma_px;
        const sp2 = this.sigma_p * this.sigma_p;
        const cxp = this.cov_xp;
        const cpxp = this.cov_pxp;
        const cxpx = this.cov_xpx;

        if (sp2 <= 0.0) {
            const val = sx2 * spx2 - cxpx * cxpx;
            return Math.sqrt(val > 0 ? val : 1e-30);
        }

        const term1 = sx2 - (cxp * cxp) / sp2;
        const term2 = spx2 - (cpxp * cpxp) / sp2;
        const term3 = cxpx - (cxp * cpxp) / sp2; // Wait, let's verify if Python has cxp * cpxp or cxp * cxp.
        // Python: self.cov_xpx - self.cov_xp * self.cov_pxp / self.sigma_p.square()
        // Yes, term3 is: cxpx - (cxp * cpxp) / sp2!
        const term3_act = cxpx - (cxp * cpxp) / sp2;

        const val = term1 * term2 - term3_act * term3_act;
        return Math.sqrt(val > 0 ? val : 1e-30);
    }

    get normalized_emittance_x() {
        return this.emittance_x * this.relativistic_beta * this.relativistic_gamma;
    }

    get beta_x() {
        const sp2 = this.sigma_p * this.sigma_p;
        const term = sp2 > 0 ? (this.cov_xp * this.cov_xp) / sp2 : 0.0;
        return (this.sigma_x * this.sigma_x - term) / this.emittance_x;
    }

    get alpha_x() {
        const sp2 = this.sigma_p * this.sigma_p;
        const term = sp2 > 0 ? (this.cov_xp * this.cov_pxp) / sp2 : 0.0;
        return -(this.cov_xpx - term) / this.emittance_x;
    }

    get projected_emittance_y() {
        const sy = this.sigma_y;
        const spy = this.sigma_py;
        const cypy = this.cov_ypy;
        const val = sy * sy * spy * spy - cypy * cypy;
        return Math.sqrt(val > 0 ? val : 0.0);
    }

    get emittance_y() {
        const sy2 = this.sigma_y * this.sigma_y;
        const spy2 = this.sigma_py * this.sigma_py;
        const sp2 = this.sigma_p * this.sigma_p;
        const cyp = this.cov_yp;
        const cpyp = this.cov_pyp;
        const cypy = this.cov_ypy;

        if (sp2 <= 0.0) {
            const val = sy2 * spy2 - cypy * cypy;
            return Math.sqrt(val > 0 ? val : 1e-30);
        }

        const term1 = sy2 - (cyp * cyp) / sp2;
        const term2 = spy2 - (cpyp * cpyp) / sp2;
        const term3 = cypy - (cyp * cpyp) / sp2;

        const val = term1 * term2 - term3 * term3;
        return Math.sqrt(val > 0 ? val : 1e-30);
    }

    get normalized_emittance_y() {
        return this.emittance_y * this.relativistic_beta * this.relativistic_gamma;
    }

    get beta_y() {
        const sp2 = this.sigma_p * this.sigma_p;
        const term = sp2 > 0 ? (this.cov_yp * this.cov_yp) / sp2 : 0.0;
        return (this.sigma_y * this.sigma_y - term) / this.emittance_y;
    }

    get alpha_y() {
        const sp2 = this.sigma_p * this.sigma_p;
        const term = sp2 > 0 ? (this.cov_yp * this.cov_pyp) / sp2 : 0.0;
        return -(this.cov_ypy - term) / this.emittance_y;
    }

    get dispersion_x() {
        const sp2 = this.sigma_p * this.sigma_p;
        return sp2 > 0 ? this.cov_xp / sp2 : 0.0;
    }

    get dispersion_px() {
        const sp2 = this.sigma_p * this.sigma_p;
        return sp2 > 0 ? this.cov_pxp / sp2 : 0.0;
    }

    get dispersion_y() {
        const sp2 = this.sigma_p * this.sigma_p;
        return sp2 > 0 ? this.cov_yp / sp2 : 0.0;
    }

    get dispersion_py() {
        const sp2 = this.sigma_p * this.sigma_p;
        return sp2 > 0 ? this.cov_pyp / sp2 : 0.0;
    }
}
