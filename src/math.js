// math.js - Core math utilities for cheetah.js

export function zeros(r, c) {
    const mat = new Array(r);
    for (let i = 0; i < r; i++) {
        mat[i] = new Float64Array(c);
    }
    return mat;
}

export function identity(n) {
    const mat = zeros(n, n);
    for (let i = 0; i < n; i++) {
        mat[i][i] = 1.0;
    }
    return mat;
}

export function transpose(A) {
    const r = A.length;
    const c = A[0].length;
    const B = zeros(c, r);
    for (let i = 0; i < r; i++) {
        for (let j = 0; j < c; j++) {
            B[j][i] = A[i][j];
        }
    }
    return B;
}

export function matmul(A, B) {
    const rA = A.length;
    const cA = A[0].length;
    const rB = B.length;
    const cB = B[0].length;
    if (cA !== rB) {
        throw new Error(`Dimension mismatch in matmul: ${cA} != ${rB}`);
    }
    const C = zeros(rA, cB);
    for (let i = 0; i < rA; i++) {
        for (let j = 0; j < cB; j++) {
            let sum = 0.0;
            for (let k = 0; k < cA; k++) {
                sum += A[i][k] * B[k][j];
            }
            C[i][j] = sum;
        }
    }
    return C;
}

export function matvec(A, v) {
    const r = A.length;
    const c = A[0].length;
    if (c !== v.length) {
        throw new Error(`Dimension mismatch in matvec: ${c} != ${v.length}`);
    }
    const out = new Float64Array(r);
    for (let i = 0; i < r; i++) {
        let sum = 0.0;
        for (let j = 0; j < c; j++) {
            sum += A[i][j] * v[j];
        }
        out[i] = sum;
    }
    return out;
}

export function clone2d(A) {
    const r = A.length;
    const c = A[0].length;
    const B = zeros(r, c);
    for (let i = 0; i < r; i++) {
        B[i].set(A[i]);
    }
    return B;
}

// 6x6 Cholesky Decomposition
// A = L * L^T. Returns L.
export function cholesky6(A) {
    const n = 6;
    const L = zeros(n, n);
    for (let i = 0; i < n; i++) {
        for (let j = 0; j <= i; j++) {
            let sum = 0.0;
            for (let k = 0; k < j; k++) {
                sum += L[i][k] * L[j][k];
            }
            if (i === j) {
                const val = A[i][i] - sum;
                L[i][j] = Math.sqrt(val > 0 ? val : 1e-30); // clamp for numerical stability
            } else {
                L[i][j] = (A[i][j] - sum) / L[j][j];
            }
        }
    }
    return L;
}

// Solve L * y = b where L is lower triangular of size n x n
export function solveTriangularL(L, b) {
    const n = L.length;
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
        let sum = 0.0;
        for (let k = 0; k < i; k++) {
            sum += L[i][k] * y[k];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }
    return y;
}

// Solve L^T * y = b where L is lower triangular, so L^T is upper triangular
export function solveTriangularLT(L, b) {
    const n = L.length;
    const y = new Float64Array(n);
    for (let i = n - 1; i >= 0; i--) {
        let sum = 0.0;
        for (let k = i + 1; k < n; k++) {
            sum += L[k][i] * y[k];
        }
        y[i] = (b[i] - sum) / L[i][i];
    }
    return y;
}

// Solve L * Y = B where B is an identity matrix to compute L^{-1}
export function invertTriangularL(L) {
    const n = L.length;
    const inv = zeros(n, n);
    for (let j = 0; j < n; j++) {
        const ej = new Float64Array(n);
        ej[j] = 1.0;
        const col = solveTriangularL(L, ej);
        for (let i = 0; i < n; i++) {
            inv[i][j] = col[i];
        }
    }
    return inv;
}

// Box-Muller transform for standard normal samples
let spareNormal = null;
export function randn() {
    if (spareNormal !== null) {
        const val = spareNormal;
        spareNormal = null;
        return val;
    }
    let u = 0, v = 0;
    while(u === 0) u = Math.random(); 
    while(v === 0) v = Math.random();
    const r = Math.sqrt(-2.0 * Math.log(u));
    const theta = 2.0 * Math.PI * v;
    spareNormal = r * Math.sin(theta);
    return r * Math.cos(theta);
}

// Weighted Statistics
export function weightedMean(values, weights) {
    let sumVal = 0.0;
    let sumWeight = 0.0;
    for (let i = 0; i < values.length; i++) {
        sumVal += values[i] * weights[i];
        sumWeight += weights[i];
    }
    return sumWeight > 0 ? sumVal / sumWeight : 0.0;
}

export function weightedVariance(values, weights, mean) {
    let sumWeight = 0.0;
    let sumSquareWeight = 0.0;
    for (let i = 0; i < weights.length; i++) {
        sumWeight += weights[i];
        sumSquareWeight += weights[i] * weights[i];
    }
    if (sumWeight <= 0) return 0.0;
    const correctionFactor = sumWeight - sumSquareWeight / sumWeight;
    if (correctionFactor <= 0) return 0.0;

    let sumDiffSq = 0.0;
    for (let i = 0; i < values.length; i++) {
        const diff = values[i] - mean;
        sumDiffSq += weights[i] * diff * diff;
    }
    return sumDiffSq / correctionFactor;
}

export function weightedStd(values, weights, mean) {
    return Math.sqrt(weightedVariance(values, weights, mean));
}

export function weightedCovariance(x, y, weights, meanX, meanY) {
    let sumWeight = 0.0;
    let sumSquareWeight = 0.0;
    for (let i = 0; i < weights.length; i++) {
        sumWeight += weights[i];
        sumSquareWeight += weights[i] * weights[i];
    }
    if (sumWeight <= 0) return 0.0;
    const correctionFactor = sumWeight - sumSquareWeight / sumWeight;
    if (correctionFactor <= 0) return 0.0;

    let sumDiffProduct = 0.0;
    for (let i = 0; i < x.length; i++) {
        sumDiffProduct += weights[i] * (x[i] - meanX) * (y[i] - meanY);
    }
    return sumDiffProduct / correctionFactor;
}

// Compute 6x6 covariance matrix of particles (N x 7, only first 6 columns used)
export function weightedCovarianceMatrix6(particles, weights) {
    const N = particles.length;
    const means = new Float64Array(6);
    const sumW = weights.reduce((a, b) => a + b, 0.0);
    
    for (let c = 0; c < 6; c++) {
        let sumVal = 0.0;
        for (let i = 0; i < N; i++) {
            sumVal += particles[i][c] * weights[i];
        }
        means[c] = sumW > 0 ? sumVal / sumW : 0.0;
    }

    let sumSquareW = 0.0;
    for (let i = 0; i < N; i++) {
        sumSquareW += weights[i] * weights[i];
    }
    const correctionFactor = sumW - sumSquareW / sumW;

    const cov = zeros(6, 6);
    if (correctionFactor <= 0) return cov;

    for (let r = 0; r < 6; r++) {
        for (let c = r; c < 6; c++) {
            let sumDiffProduct = 0.0;
            for (let i = 0; i < N; i++) {
                sumDiffProduct += weights[i] * (particles[i][r] - means[r]) * (particles[i][c] - means[c]);
            }
            const val = sumDiffProduct / correctionFactor;
            cov[r][c] = val;
            cov[c][r] = val;
        }
    }
    return cov;
}

// Match distribution moments exactly
export function matchDistributionMoments(samples, targetMu, targetCov, weights = null) {
    const N = samples.length;
    if (N === 0) return samples;
    
    if (!weights) {
        weights = new Float64Array(N).fill(1.0);
    }
    const sumW = weights.reduce((a, b) => a + b, 0.0);

    // Compute sample mean and covariance
    const sampleMu = new Float64Array(6);
    for (let c = 0; c < 6; c++) {
        let sumVal = 0.0;
        for (let i = 0; i < N; i++) {
            sumVal += samples[i][c] * weights[i];
        }
        sampleMu[c] = sumW > 0 ? sumVal / sumW : 0.0;
    }

    const sampleCov = weightedCovarianceMatrix6(samples, weights);
    
    // Cholesky decomposition of sample covariance
    const L_sample = cholesky6(sampleCov);
    // Invert L_sample to whiten the samples
    const L_sample_inv = invertTriangularL(L_sample);

    // Cholesky decomposition of target covariance
    // targetCov is 7x7 in parameter beam, but we only need 6x6 here
    const targetCov6 = zeros(6, 6);
    for (let i = 0; i < 6; i++) {
        for (let j = 0; j < 6; j++) {
            targetCov6[i][j] = targetCov[i][j];
        }
    }
    const L_target = cholesky6(targetCov6);

    // Whiten samples and map to target distribution:
    // x_new = L_target * L_sample_inv * (x - sampleMu) + targetMu
    const transform = matmul(L_target, L_sample_inv);

    const matched = new Array(N);
    for (let i = 0; i < N; i++) {
        const x_centered = new Float64Array(6);
        for (let c = 0; c < 6; c++) {
            x_centered[c] = samples[i][c] - sampleMu[c];
        }
        const x_transformed = matvec(transform, x_centered);
        const p_new = new Float64Array(7);
        for (let c = 0; c < 6; c++) {
            p_new[c] = x_transformed[c] + targetMu[c];
        }
        p_new[6] = 1.0; // 7th coordinate is 1
        matched[i] = p_new;
    }
    return matched;
}

// ----------------------------------------------------
// Autograd forward-only helper functions
// ----------------------------------------------------

export function si(x) {
    if (x === 0.0) return 1.0;
    if (x > 0.0) {
        const y = Math.sqrt(x);
        return Math.sin(y) / y;
    } else {
        const y = Math.sqrt(-x);
        return Math.sinh(y) / y;
    }
}

export function cx_helper(x) {
    if (x === 0.0) return 1.0;
    if (x > 0.0) return Math.cos(Math.sqrt(x));
    return Math.cosh(Math.sqrt(-x));
}

export function log1pdiv(x) {
    if (Math.abs(x) < 1e-12) return 1.0;
    return Math.log1p(x) / x;
}

export function si1mdiv(x) {
    if (Math.abs(x) < 1e-12) return 1.0 / 6.0;
    const sx = si(x);
    return (1.0 - sx) / x;
}

export function sicos1mdiv(x) {
    if (Math.abs(x) < 1e-12) return 1.0 / 6.0;
    const sx = si(x);
    const cx = cx_helper(x);
    return (1.0 - sx * cx) / x;
}

export function sipsicos3mdiv(x) {
    if (Math.abs(x) < 1e-12) return 0.0;
    const sx = si(x);
    const cx = cx_helper(x);
    return (3.0 - 4.0 * sx + sx * cx) / (2.0 * x);
}

export function cossqrtmcosdivdiff(a, b) {
    if (Math.abs(a - b) < 1e-12) {
        return 0.5 * si(a);
    }
    const ca = cx_helper(a);
    const cb = cx_helper(b);
    return (cb - ca) / (a - b);
}

export function simsidivdiff(a, b) {
    if (Math.abs(a - b) < 1e-12) {
        if (Math.abs(b) < 1e-12) return 1.0 / 6.0;
        const sb = si(b);
        const cb = cx_helper(b);
        return 0.5 * (sb - cb) / b;
    }
    const sa = si(a);
    const sb = si(b);
    return (sa - sb) / (b - a);
}

export function si2msi2divdiff(a, b) {
    if (Math.abs(a - b) < 1e-12) {
        if (Math.abs(b) < 1e-12) return 1.0 / 3.0;
        const sb = si(b);
        const cb = cx_helper(b);
        return (1.0 - cb * cb - b * sb * cb) / (b * b);
    }
    const sa = si(a);
    const sb = si(b);
    return (sb * sb - sa * sa) / (a - b);
}

export function sqrta2minusbdiva(a, b) {
    if (Math.abs(b) < 1e-12) {
        return 1.0 / (2.0 * a);
    }
    return (Math.sqrt(a * a + b) - a) / b;
}
