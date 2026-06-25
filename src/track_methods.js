// track_methods.js - Universal transfer map functions for first-order and second-order tracking

import { zeros, identity, si1mdiv, sicos1mdiv, sipsicos3mdiv, cossqrtmcosdivdiff, simsidivdiff, si2msi2divdiff } from "./math.js";

export function computeRelativisticFactors(energy, mass_eV) {
    const gamma = energy / mass_eV;
    const igamma2 = 1.0 / (gamma * gamma);
    const beta = Math.sqrt(Math.max(0.0, 1.0 - igamma2));
    return { gamma, igamma2, beta };
}

export function base_rmatrix(length, k1, hx, species, energy = 0.0) {
    const { igamma2, beta } = computeRelativisticFactors(energy, species.mass_eV);

    const kx2 = k1 + hx * hx;
    const ky2 = -k1;

    let cx = 1.0, sx = length;
    if (kx2 > 0.0) {
        const kx = Math.sqrt(kx2);
        cx = Math.cos(kx * length);
        sx = Math.sin(kx * length) / kx;
    } else if (kx2 < 0.0) {
        const kx = Math.sqrt(-kx2);
        cx = Math.cosh(kx * length);
        sx = Math.sinh(kx * length) / kx;
    }

    let cy = 1.0, sy = length;
    if (ky2 > 0.0) {
        const ky = Math.sqrt(ky2);
        cy = Math.cos(ky * length);
        sy = Math.sin(ky * length) / ky;
    } else if (ky2 < 0.0) {
        const ky = Math.sqrt(-ky2);
        cy = Math.cosh(ky * length);
        sy = Math.sinh(ky * length) / ky;
    }

    let r = 1.0;
    if (kx2 > 0.0) {
        const kx = Math.sqrt(kx2);
        r = Math.sin(kx * length / 2.0) / (kx * length / 2.0);
    } else if (kx2 < 0.0) {
        const kx = Math.sqrt(-kx2);
        r = Math.sinh(kx * length / 2.0) / (kx * length / 2.0);
    }
    const dx = hx * 0.5 * length * length * r * r;

    const r56 = (hx * hx * Math.pow(length, 3) * si1mdiv(kx2 * length * length) / (beta * beta)) - length / (beta * beta) * igamma2;

    const R = identity(7);
    R[0][0] = cx;
    R[0][1] = sx;
    R[0][5] = dx / beta;
    R[1][0] = -kx2 * sx;
    R[1][1] = cx;
    R[1][5] = sx * hx / beta;
    R[2][2] = cy;
    R[2][3] = sy;
    R[3][2] = -ky2 * sy;
    R[3][3] = cy;
    R[4][0] = sx * hx / beta;
    R[4][1] = dx / beta;
    R[4][5] = r56;

    return R;
}

export function base_ttensor(length, k1, k2, hx, species, energy = 0.0) {
    const { igamma2, beta } = computeRelativisticFactors(energy, species.mass_eV);

    const kx2 = k1 + hx * hx;
    const ky2 = -k1;

    let cx = 1.0, sx = length;
    if (kx2 > 0.0) {
        const kx = Math.sqrt(kx2);
        cx = Math.cos(kx * length);
        sx = Math.sin(kx * length) / kx;
    } else if (kx2 < 0.0) {
        const kx = Math.sqrt(-kx2);
        cx = Math.cosh(kx * length);
        sx = Math.sinh(kx * length) / kx;
    }

    let cy = 1.0, sy = length;
    if (ky2 > 0.0) {
        const ky = Math.sqrt(ky2);
        cy = Math.cos(ky * length);
        sy = Math.sin(ky * length) / ky;
    } else if (ky2 < 0.0) {
        const ky = Math.sqrt(-ky2);
        cy = Math.cosh(ky * length);
        sy = Math.sinh(ky * length) / ky;
    }

    let r = 1.0;
    if (kx2 > 0.0) {
        const kx = Math.sqrt(kx2);
        r = Math.sin(kx * length / 2.0) / (kx * length / 2.0);
    } else if (kx2 < 0.0) {
        const kx = Math.sqrt(-kx2);
        r = Math.sinh(kx * length / 2.0) / (kx * length / 2.0);
    }
    const dx = 0.5 * length * length * r * r;

    const fx = Math.pow(length, 3) * si1mdiv(kx2 * length * length);
    const f2y = Math.pow(length, 3) * sicos1mdiv(ky2 * length * length);

    const j1 = fx;
    const j2 = Math.pow(length, 3) * sipsicos3mdiv(kx2 * length * length);
    
    let j3 = 0.0;
    if (kx2 !== 0.0) {
        j3 = (15.0 * length - 22.5 * sx + 9.0 * sx * cx - 1.5 * sx * cx * cx + kx2 * sx * sx * sx) / (6.0 * Math.pow(kx2, 3));
    } else {
        j3 = Math.pow(length, 7) / 56.0;
    }

    const j_denominator = kx2 - 4.0 * ky2;
    const jc = length * length * cossqrtmcosdivdiff(kx2 * length * length, ky2 * length * length);
    const js = Math.pow(length, 3) * simsidivdiff(kx2 * length * length, ky2 * length * length);
    const jd = Math.pow(length, 4) * si2msi2divdiff(kx2 * length * length, ky2 * length * length);
    
    let jf = 0.0;
    if (j_denominator !== 0.0) {
        jf = (f2y - fx) / j_denominator;
    } else {
        jf = Math.pow(length, 5) / 120.0;
    }

    const khk = k2 + 2.0 * hx * k1;

    // Create 3D tensor 7x7x7
    const T = new Array(7);
    for (let i = 0; i < 7; i++) {
        T[i] = zeros(7, 7);
    }

    T[0][0][0] = -khk * (sx * sx + dx) / 6.0 - 0.5 * hx * kx2 * sx * sx;
    T[0][0][1] = 2.0 * (-khk * sx * dx / 6.0 + 0.5 * hx * sx * cx);
    T[0][1][1] = -khk * dx * dx / 6.0 + 0.5 * hx * dx * cx;
    T[0][0][5] = 2.0 * (-hx / 12.0 / beta * khk * (3.0 * sx * j1 - dx * dx) + 0.5 * hx * hx / beta * sx * sx + 0.25 / beta * k1 * length * sx);
    T[0][1][5] = 2.0 * (-hx / 12.0 / beta * khk * (sx * dx * dx - 2.0 * cx * j2) + 0.25 * hx * hx / beta * (sx * dx + cx * j1) - 0.25 / beta * (sx + length * cx));
    T[0][5][5] = -hx * hx / 6.0 / (beta * beta) * khk * (dx * dx * dx - 2.0 * sx * j2) + 0.5 * Math.pow(hx, 3) / (beta * beta) * sx * j1 - 0.5 * hx / (beta * beta) * length * sx - 0.5 * hx / (beta * beta) * igamma2 * dx;
    
    T[0][2][2] = k1 * k2 * jd + 0.5 * (k2 + hx * k1) * dx;
    T[0][2][3] = 2.0 * (0.5 * k2 * js);
    T[0][3][3] = k2 * jd - 0.5 * hx * dx;

    T[1][0][0] = -khk * sx * (1.0 + 2.0 * cx) / 6.0;
    T[1][0][1] = -2.0 * khk * dx * (1.0 + 2.0 * cx) / 6.0;
    T[1][1][1] = -khk * sx * dx / 3.0 - 0.5 * hx * sx;
    T[1][0][5] = 2.0 * (-hx / 12.0 / beta * khk * (3.0 * cx * j1 + sx * dx) - 0.25 / beta * k1 * (sx - length * cx));
    T[1][1][5] = 2.0 * (-hx / 12.0 / beta * khk * (3.0 * sx * j1 + dx * dx) + 0.25 / beta * k1 * length * sx);
    T[1][5][5] = -hx * hx / 6.0 / (beta * beta) * khk * (sx * dx * dx - 2.0 * cx * j2) - 0.5 * hx / (beta * beta) * k1 * (cx * j1 - sx * dx) - 0.5 * hx / (beta * beta) * igamma2 * sx;
    
    T[1][2][2] = k1 * k2 * js + 0.5 * (k2 + hx * k1) * sx;
    T[1][2][3] = 2.0 * (0.5 * k2 * jc);
    T[1][3][3] = k2 * js - 0.5 * hx * sx;

    T[2][0][2] = 2.0 * (0.5 * k2 * (cy * jc - 2.0 * k1 * sy * js) + 0.5 * hx * k1 * sx * sy);
    T[2][0][3] = 2.0 * (0.5 * k2 * (sy * jc - 2.0 * cy * js) + 0.5 * hx * sx * cy);
    T[2][1][2] = 2.0 * (0.5 * k2 * (cy * js - 2.0 * k1 * sy * jd) + 0.5 * hx * k1 * dx * sy);
    T[2][1][3] = 2.0 * (0.5 * k2 * (sy * js - 2.0 * cy * jd) + 0.5 * hx * dx * cy);
    T[2][2][5] = 2.0 * (0.5 * hx / beta * k2 * (cy * jd - 2.0 * k1 * sy * jf) + 0.5 * hx * hx / beta * k1 * j1 * sy - 0.25 / beta * k1 * length * sy);
    T[2][3][5] = 2.0 * (0.5 * hx / beta * k2 * (sy * jd - 2.0 * cy * jf) + 0.5 * hx * hx / beta * j1 * cy - 0.25 / beta * (sy + length * cy));

    T[3][0][2] = 2.0 * (0.5 * k1 * k2 * (2.0 * cy * js - sy * jc) + 0.5 * (k2 + hx * k1) * sx * cy);
    T[3][0][3] = 2.0 * (0.5 * k2 * (2.0 * k1 * sy * js - cy * jc) + 0.5 * (k2 + hx * k1) * sx * sy);
    T[3][1][2] = 2.0 * (0.5 * k1 * k2 * (2.0 * cy * jd - sy * js) + 0.5 * (k2 + hx * k1) * dx * cy);
    T[3][1][3] = 2.0 * (0.5 * k2 * (2.0 * k1 * sy * jd - cy * js) + 0.5 * (k2 + hx * k1) * dx * sy);
    T[3][2][5] = 2.0 * (0.5 * hx / beta * k1 * k2 * (2.0 * cy * jf - sy * jd) + 0.5 * hx / beta * (k2 + hx * k1) * j1 * cy + 0.25 / beta * k1 * (sy - length * cy));
    T[3][3][5] = 2.0 * (0.5 * hx / beta * k2 * (2.0 * k1 * sy * jf - cy * jd) + 0.5 * hx / beta * (k2 + hx * k1) * j1 * sy - 0.25 / beta * k1 * length * sy);

    T[4][0][0] = -(hx / 12.0 / beta * khk * (sx * dx + 3.0 * j1) - 0.25 / beta * k1 * (length - sx * cx));
    T[4][0][1] = -2.0 * (hx / 12.0 / beta * khk * dx * dx + 0.25 / beta * k1 * sx * sx);
    T[4][1][1] = -(hx / 6.0 / beta * khk * j2 - 0.5 / beta * sx - 0.25 / beta * k1 * (j1 - sx * dx));
    T[4][0][5] = -2.0 * (hx * hx / 12.0 / (beta * beta) * khk * (3.0 * dx * j1 - 4.0 * j2) + 0.25 * hx / (beta * beta) * k1 * j1 * (1.0 + cx) + 0.5 * hx / (beta * beta) * igamma2 * sx);
    T[4][1][5] = -2.0 * (hx * hx / 12.0 / (beta * beta) * khk * (dx * dx * dx - 2.0 * sx * j2) + 0.25 * hx / (beta * beta) * k1 * sx * j1 + 0.5 * hx / (beta * beta) * igamma2 * dx);
    T[4][5][5] = -(Math.pow(hx, 3) / 6.0 / Math.pow(beta, 3) * khk * (3.0 * j3 - 2.0 * dx * j2) + hx * hx / 6.0 / Math.pow(beta, 3) * k1 * (sx * dx * dx - j2 * (1.0 + 2.0 * cx)) + 1.5 / Math.pow(beta, 3) * igamma2 * (hx * hx * j1 - length));
    
    T[4][2][2] = -(-hx / beta * k1 * k2 * jf - 0.5 * hx / beta * (k2 + hx * k1) * j1 + 0.25 / beta * k1 * (length - cy * sy));
    T[4][2][3] = -2.0 * (-0.5 * hx / beta * k2 * jd - 0.25 / beta * k1 * sy * sy);
    T[4][3][3] = -(-hx / beta * k2 * jf + 0.5 * hx * hx / beta * j1 - 0.25 / beta * (length + cy * sy));

    return T;
}

export function drift_matrix(length, energy, species) {
    const { igamma2, beta } = computeRelativisticFactors(energy, species.mass_eV);
    const tm = identity(7);
    tm[0][1] = length;
    tm[2][3] = length;
    tm[4][5] = -length / (beta * beta) * igamma2;
    return tm;
}
