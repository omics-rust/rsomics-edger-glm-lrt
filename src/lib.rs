//! edgeR glmFit + glmLRT: negative-binomial GLM fit and likelihood-ratio test
//! of a coefficient or contrast. Method: McCarthy, Chen & Smyth (2012), NAR
//! 40:4288-4297. Per gene the full design is fit by NB IRLS (log link), the LR
//! statistic is the drop in NB deviance from the reduced model, and its p-value
//! is the chi-square upper tail on the tested degrees of freedom. logFC is the
//! tested coefficient/contrast from edgeR's predFC (prior.count 0.125 shrinkage)
//! rescaled to log2; logCPM is aveLogCPM (prior.count 2, dispersion 0.05).

mod special;

use std::fs::File;
use std::io::{BufRead, BufReader, Write};
use std::path::Path;

use rsomics_common::{Result, RsomicsError};

const LN2: f64 = std::f64::consts::LN_2;
const PRIOR_COUNT: f64 = 0.125;
const AVELOGCPM_PRIOR: f64 = 2.0;
const AVELOGCPM_DISP: f64 = 0.05;

pub struct Matrix {
    pub header: String,
    pub genes: Vec<String>,
    pub counts: Vec<f64>,
    pub n_samples: usize,
}

impl Matrix {
    pub fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty count matrix".into()))?
            .map_err(RsomicsError::Io)?;
        let n_samples = header.split('\t').count() - 1;
        if n_samples == 0 {
            return Err(RsomicsError::InvalidInput(
                "count matrix has no sample columns".into(),
            ));
        }
        let mut genes = Vec::new();
        let mut counts = Vec::new();
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.is_empty() {
                continue;
            }
            let mut fields = line.split('\t');
            let gene = fields
                .next()
                .ok_or_else(|| RsomicsError::InvalidInput("row without a gene id".into()))?;
            genes.push(gene.to_string());
            let before = counts.len();
            for f in fields {
                counts.push(f.parse::<f64>().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-numeric count '{f}' for gene {gene}"))
                })?);
            }
            if counts.len() - before != n_samples {
                return Err(RsomicsError::InvalidInput(format!(
                    "gene {gene}: {} values, header has {n_samples} samples",
                    counts.len() - before
                )));
            }
        }
        Ok(Self {
            header,
            genes,
            counts,
            n_samples,
        })
    }

    pub fn n_genes(&self) -> usize {
        self.genes.len()
    }
    fn row(&self, g: usize) -> &[f64] {
        &self.counts[g * self.n_samples..(g + 1) * self.n_samples]
    }
}

/// Design matrix: row-major n_samples × n_coef, plus column names from the header.
pub struct Design {
    pub data: Vec<f64>,
    pub n_samples: usize,
    pub n_coef: usize,
    pub coef_names: Vec<String>,
}

impl Design {
    fn load(path: &Path) -> Result<Self> {
        let file = File::open(path)
            .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
        let mut lines = BufReader::new(file).lines();
        let header = lines
            .next()
            .ok_or_else(|| RsomicsError::InvalidInput("empty design matrix".into()))?
            .map_err(RsomicsError::Io)?;
        let coef_names: Vec<String> = header.split('\t').map(str::to_string).collect();
        if coef_names.is_empty() {
            return Err(RsomicsError::InvalidInput(
                "design has no coefficient columns".into(),
            ));
        }
        let n_coef = coef_names.len();
        let mut data = Vec::new();
        let mut n_samples = 0;
        for line in lines {
            let line = line.map_err(RsomicsError::Io)?;
            if line.is_empty() {
                continue;
            }
            let before = data.len();
            for f in line.split('\t') {
                data.push(f.parse::<f64>().map_err(|_| {
                    RsomicsError::InvalidInput(format!("non-numeric design value '{f}'"))
                })?);
            }
            if data.len() - before != n_coef {
                return Err(RsomicsError::InvalidInput(format!(
                    "design row {n_samples}: {} values, header has {n_coef} columns",
                    data.len() - before
                )));
            }
            n_samples += 1;
        }
        Ok(Self {
            data,
            n_samples,
            n_coef,
            coef_names,
        })
    }

    fn row(&self, s: usize) -> &[f64] {
        &self.data[s * self.n_coef..(s + 1) * self.n_coef]
    }
}

fn load_norm_factors(path: &Path, n_samples: usize) -> Result<Vec<f64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut factors = Vec::with_capacity(n_samples);
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let val = line.rsplit('\t').next().unwrap_or(line);
        factors.push(
            val.parse::<f64>().map_err(|_| {
                RsomicsError::InvalidInput(format!("non-numeric norm factor '{val}'"))
            })?,
        );
    }
    if factors.len() != n_samples {
        return Err(RsomicsError::InvalidInput(format!(
            "{} norm factors for {n_samples} samples",
            factors.len()
        )));
    }
    Ok(factors)
}

fn load_contrast(path: &Path, n_coef: usize) -> Result<Vec<f64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut c = Vec::with_capacity(n_coef);
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let val = line.rsplit('\t').next().unwrap_or(line);
        c.push(
            val.parse::<f64>()
                .map_err(|_| RsomicsError::InvalidInput(format!("non-numeric contrast '{val}'")))?,
        );
    }
    if c.len() != n_coef {
        return Err(RsomicsError::InvalidInput(format!(
            "contrast has {} entries, design has {n_coef} coefficients",
            c.len()
        )));
    }
    Ok(c)
}

const MAXIT: usize = 200;
const TOL: f64 = 1e-10;

/// NB deviance, edgeR nbinomDeviance: 2·Σ [ y·log(y/μ) − (y+1/φ)·log((y+1/φ)/(μ+1/φ)) ].
fn nb_deviance(y: &[f64], mu: &[f64], dispersion: f64) -> f64 {
    let mut dev = 0.0;
    for (&yi, &mui) in y.iter().zip(mu) {
        let r = if dispersion > 0.0 {
            1.0 / dispersion
        } else {
            f64::INFINITY
        };
        let term_y = if yi > 0.0 { yi * (yi / mui).ln() } else { 0.0 };
        let term = if dispersion > 0.0 {
            (yi + r) * ((yi + r) / (mui + r)).ln()
        } else {
            mui - yi
        };
        dev += if dispersion > 0.0 {
            term_y - term
        } else {
            term_y - (yi - mui)
        };
    }
    2.0 * dev
}

/// Per-thread scratch covering all three fits a gene needs (full + reduced + the
/// prior-augmented full), so a worker allocates once, not per gene.
struct GeneScratch {
    full: Scratch,
    reduced: Scratch,
    row_aug: Vec<f64>,
}

/// Reusable per-thread scratch so the hot IRLS loop allocates nothing per gene.
struct Scratch {
    beta: Vec<f64>,
    trial: Vec<f64>,
    mu: Vec<f64>,
    mu_t: Vec<f64>,
    xtwx: Vec<f64>,
    a: Vec<f64>,
    xtr: Vec<f64>,
    rhs: Vec<f64>,
    step: Vec<f64>,
}

impl Scratch {
    fn new(n: usize, p: usize) -> Self {
        Scratch {
            beta: vec![0.0; p],
            trial: vec![0.0; p],
            mu: vec![0.0; n],
            mu_t: vec![0.0; n],
            xtwx: vec![0.0; p * p],
            a: vec![0.0; p * p],
            xtr: vec![0.0; p],
            rhs: vec![0.0; p],
            step: vec![0.0; p],
        }
    }
}

/// IRLS NB GLM with log link and per-sample offset, Fisher scoring with a
/// Levenberg ridge (matches edgeR's mglmLevenberg). `tight` converges to the
/// β-MLE (for logFC); else stops on edgeR's tol=1e-6 deviance criterion (the LR
/// fits, where only the deviance matters). Returns (β, deviance).
fn fit_nb_glm(
    y: &[f64],
    design: &Design,
    offset: &[f64],
    dispersion: f64,
    tight: bool,
    start_dir: &[f64],
    sc: &mut Scratch,
) -> (Vec<f64>, f64) {
    let n = design.n_samples;
    let p = design.n_coef;

    // edgeR start.method="null": intercept-only NB fit sets a common mean; the
    // OLS projection of that constant predictor is b0·(XᵀX)⁻¹Xᵀ1 = b0·start_dir.
    let b0 = mglm_one_group(y, offset, dispersion);
    for (b, &d) in sc.beta[..p].iter_mut().zip(start_dir) {
        *b = b0 * d;
    }

    eta_mu(design, offset, &sc.beta[..p], &mut sc.mu);
    let mut dev = nb_deviance(y, &sc.mu[..n], dispersion);
    let mut lambda = 0.0f64;

    for _ in 0..MAXIT {
        for v in sc.xtwx[..p * p].iter_mut() {
            *v = 0.0;
        }
        for v in sc.xtr[..p].iter_mut() {
            *v = 0.0;
        }
        let rows = design.data[..n * p].chunks_exact(p);
        for ((&yi, &mui), xr) in y.iter().zip(&sc.mu[..n]).zip(rows) {
            let denom = 1.0 + dispersion * mui;
            let w = mui / denom;
            let resid = (yi - mui) / denom;
            for (j, &xj) in xr.iter().enumerate() {
                sc.xtr[j] += xj * resid;
                let xjw = xj * w;
                let rowj = &mut sc.xtwx[j * p..j * p + p];
                for (rk, &xk) in rowj.iter_mut().zip(xr) {
                    *rk += xjw * xk;
                }
            }
        }
        let mut accepted = false;
        for _ in 0..20 {
            sc.a[..p * p].copy_from_slice(&sc.xtwx[..p * p]);
            for d in 0..p {
                sc.a[d * p + d] += lambda * sc.xtwx[d * p + d].max(1e-6);
            }
            sc.rhs[..p].copy_from_slice(&sc.xtr[..p]);
            if !solve(&mut sc.a[..p * p], &mut sc.rhs[..p], &mut sc.step[..p], p) {
                lambda = if lambda == 0.0 { 1.0 } else { lambda * 2.0 };
                continue;
            }
            for (t, (&b, &s)) in sc.trial[..p]
                .iter_mut()
                .zip(sc.beta[..p].iter().zip(&sc.step[..p]))
            {
                *t = b + s;
            }
            eta_mu(design, offset, &sc.trial[..p], &mut sc.mu_t);
            let dev_t = nb_deviance(y, &sc.mu_t[..n], dispersion);
            if dev_t <= dev + 1e-8 * (1.0 + dev.abs()) {
                let max_step = sc.step[..p].iter().fold(0.0f64, |m, s| m.max(s.abs()));
                let improved = dev - dev_t;
                sc.beta[..p].copy_from_slice(&sc.trial[..p]);
                sc.mu[..n].copy_from_slice(&sc.mu_t[..n]);
                dev = dev_t;
                lambda *= 0.5;
                accepted = true;
                let done = if tight {
                    max_step < TOL
                } else {
                    improved < 1e-7 * (dev + 1.0)
                };
                if done {
                    return (sc.beta[..p].to_vec(), dev);
                }
                break;
            }
            lambda = if lambda == 0.0 { 1.0 } else { lambda * 4.0 };
        }
        if !accepted {
            break;
        }
    }
    (sc.beta[..p].to_vec(), dev)
}

fn eta_mu(design: &Design, offset: &[f64], beta: &[f64], mu: &mut [f64]) {
    for s in 0..design.n_samples {
        let xr = design.row(s);
        let mut eta = offset[s];
        for (&x, &b) in xr.iter().zip(beta) {
            eta += x * b;
        }
        mu[s] = eta.exp();
    }
}

/// Per-design start direction d = (XᵀX)⁻¹Xᵀ1, so the null-method start for a gene
/// is b0·d (the OLS projection of a constant linear predictor b0 onto the design).
fn start_direction(design: &Design) -> Vec<f64> {
    let p = design.n_coef;
    let mut xtx = vec![0.0f64; p * p];
    let mut xt1 = vec![0.0f64; p];
    for s in 0..design.n_samples {
        let xr = design.row(s);
        for (j, &xj) in xr.iter().enumerate() {
            xt1[j] += xj;
            let rowj = &mut xtx[j * p..j * p + p];
            for (rk, &xk) in rowj.iter_mut().zip(xr) {
                *rk += xj * xk;
            }
        }
    }
    let mut d = vec![0.0f64; p];
    if !solve(&mut xtx, &mut xt1, &mut d, p) {
        return vec![0.0f64; p];
    }
    d
}

/// Design columns that are not estimable because the design is rank-deficient
/// (limma nonEstimable). Left-to-right Gram-Schmidt keeps a column only if its
/// residual against already-kept columns clears R's qr tol=1e-7; a column
/// dependent on earlier ones is flagged, so the later of two collinear columns
/// is the one reported — matching R's `qr()` pivoting for exact dependence.
fn non_estimable(design: &Design) -> Vec<String> {
    let p = design.n_coef;
    let n = design.n_samples;
    let mut basis: Vec<Vec<f64>> = Vec::new();
    let mut ne = Vec::new();
    for k in 0..p {
        let mut col: Vec<f64> = (0..n).map(|s| design.data[s * p + k]).collect();
        let orig = col.iter().map(|x| x * x).sum::<f64>().sqrt();
        for b in &basis {
            let d: f64 = col.iter().zip(b).map(|(a, c)| a * c).sum();
            for (ci, bi) in col.iter_mut().zip(b) {
                *ci -= d * bi;
            }
        }
        let resid = col.iter().map(|x| x * x).sum::<f64>().sqrt();
        if resid <= 1e-7 * orig {
            ne.push(design.coef_names[k].clone());
        } else {
            for x in &mut col {
                *x /= resid;
            }
            basis.push(col);
        }
    }
    ne
}

/// Solve A x = b (A row-major p×p), Gaussian elimination with partial pivoting.
/// `a` and `rhs` are clobbered. Returns false if singular.
fn solve(a: &mut [f64], rhs: &mut [f64], x: &mut [f64], p: usize) -> bool {
    for col in 0..p {
        let mut piv = col;
        let mut best = a[col * p + col].abs();
        for r in (col + 1)..p {
            let v = a[r * p + col].abs();
            if v > best {
                best = v;
                piv = r;
            }
        }
        if best < 1e-12 {
            return false;
        }
        if piv != col {
            for k in 0..p {
                a.swap(col * p + k, piv * p + k);
            }
            rhs.swap(col, piv);
        }
        let d = a[col * p + col];
        for r in (col + 1)..p {
            let f = a[r * p + col] / d;
            if f == 0.0 {
                continue;
            }
            for k in col..p {
                a[r * p + k] -= f * a[col * p + k];
            }
            rhs[r] -= f * rhs[col];
        }
    }
    for col in (0..p).rev() {
        let mut s = rhs[col];
        for k in (col + 1)..p {
            s -= a[col * p + k] * x[k];
        }
        x[col] = s / a[col * p + col];
    }
    true
}

/// aveLogCPM for one gene (edgeR aveLogCPM, prior.count 2, dispersion 0.05),
/// fit by the one-group NB model against offsets that include twice the prior.
fn ave_log_cpm_gene(row: &[f64], lib: &[f64]) -> f64 {
    let mean_lib = lib.iter().sum::<f64>() / lib.len() as f64;
    let prior: Vec<f64> = lib
        .iter()
        .map(|&l| AVELOGCPM_PRIOR * l / mean_lib)
        .collect();
    let off_aug: Vec<f64> = lib
        .iter()
        .zip(&prior)
        .map(|(&l, &p)| (l + 2.0 * p).ln())
        .collect();
    let aug: Vec<f64> = row.iter().zip(&prior).map(|(&c, &p)| c + p).collect();
    let beta = mglm_one_group(&aug, &off_aug, AVELOGCPM_DISP);
    (beta + 1e6f64.ln()) / LN2
}

/// One-group NB fit (edgeR mglmOneGroup): Fisher scoring for the single
/// coefficient β where μ[j] = exp(β + offset[j]). Natural-log scale.
fn mglm_one_group(row: &[f64], offset: &[f64], dispersion: f64) -> f64 {
    let total: f64 = row.iter().sum();
    if total == 0.0 {
        return f64::NEG_INFINITY;
    }
    let mean_off = offset.iter().sum::<f64>() / offset.len() as f64;
    let mut beta = (total / row.len() as f64).ln() - mean_off;
    for _ in 0..MAXIT {
        let mut dl = 0.0;
        let mut info = 0.0;
        for (&y, &off) in row.iter().zip(offset) {
            let mu = (beta + off).exp();
            let denom = 1.0 + mu * dispersion;
            dl += (y - mu) / denom;
            info += mu / denom;
        }
        let s = dl / info;
        beta += s;
        if s.abs() < TOL {
            break;
        }
    }
    beta
}

fn bh_fdr(pvals: &[f64]) -> Vec<f64> {
    let n = pvals.len();
    let mut order: Vec<usize> = (0..n).collect();
    order.sort_by(|&a, &b| pvals[b].partial_cmp(&pvals[a]).unwrap());
    let mut adj = vec![0.0f64; n];
    let mut cummin = f64::INFINITY;
    for (rank, &i) in order.iter().enumerate() {
        let m = n - rank;
        let v = (pvals[i] * n as f64 / m as f64).min(1.0);
        cummin = cummin.min(v);
        adj[i] = cummin;
    }
    adj
}

enum Test {
    Coef(usize),
    Contrast(Vec<f64>),
}

enum Dispersion {
    Common(f64),
    PerGene(Vec<f64>),
}

pub struct GlmLrtArgs<'a> {
    pub counts: &'a Path,
    pub design: &'a Path,
    pub norm_factors: Option<&'a Path>,
    pub coef: Option<usize>,
    pub contrast: Option<&'a Path>,
    pub dispersion: f64,
    pub dispersion_file: Option<&'a Path>,
    pub fdr: bool,
}

/// Reduced design for testing `coef` (1-based): drop that column. For a contrast,
/// edgeR rotates the design so the contrast becomes one coefficient and tests it;
/// the LR is identical to dropping the rotated column, which for a single
/// non-zero-contrast direction reduces to refitting under the linear constraint.
/// We reparametrize by an orthonormal basis of the contrast's null space.
fn reduced_design_coef(design: &Design, drop: usize) -> Design {
    let p = design.n_coef;
    let mut data = Vec::with_capacity(design.n_samples * (p - 1));
    for s in 0..design.n_samples {
        let xr = design.row(s);
        for (k, &v) in xr.iter().enumerate() {
            if k != drop {
                data.push(v);
            }
        }
    }
    let coef_names: Vec<String> = design
        .coef_names
        .iter()
        .enumerate()
        .filter(|(k, _)| *k != drop)
        .map(|(_, n)| n.clone())
        .collect();
    Design {
        data,
        n_samples: design.n_samples,
        n_coef: p - 1,
        coef_names,
    }
}

/// Reparametrize the design for a contrast test (edgeR glmLRT contrast path):
/// build Q = QR(contrast) so the first new coefficient is along the contrast and
/// the rest span its complement. The full design is X·Q⁻¹·(reordered); the
/// reduced design drops the contrast direction. We return (full_reparam,
/// reduced) where full has the contrast as its first column.
fn contrast_designs(design: &Design, contrast: &[f64]) -> (Design, Design) {
    let p = design.n_coef;
    let n = design.n_samples;
    // Householder QR of the contrast (a p×1 vector) gives an orthonormal basis
    // whose first column is contrast/||contrast||; columns 2..p span its
    // orthogonal complement.
    let q = householder_basis(contrast);
    // X* = X Q. The first column of X* corresponds to the contrast direction.
    let mut full = vec![0.0f64; n * p];
    for s in 0..n {
        let xr = design.row(s);
        for j in 0..p {
            let mut v = 0.0;
            for k in 0..p {
                v += xr[k] * q[k * p + j];
            }
            full[s * p + j] = v;
        }
    }
    let full_d = Design {
        data: full.clone(),
        n_samples: n,
        n_coef: p,
        coef_names: (0..p).map(|j| format!("c{j}")).collect(),
    };
    let mut reduced = Vec::with_capacity(n * (p - 1));
    for s in 0..n {
        for j in 1..p {
            reduced.push(full[s * p + j]);
        }
    }
    let reduced_d = Design {
        data: reduced,
        n_samples: n,
        n_coef: p - 1,
        coef_names: (1..p).map(|j| format!("c{j}")).collect(),
    };
    (full_d, reduced_d)
}

/// Orthonormal p×p basis (row-major, columns are basis vectors) whose first
/// column is `v` normalized; the rest complete it via Gram-Schmidt on the
/// standard basis.
fn householder_basis(v: &[f64]) -> Vec<f64> {
    let p = v.len();
    let mut cols: Vec<Vec<f64>> = Vec::with_capacity(p);
    let norm = v.iter().map(|x| x * x).sum::<f64>().sqrt();
    cols.push(v.iter().map(|x| x / norm).collect());
    for e in 0..p {
        let mut cand = vec![0.0f64; p];
        cand[e] = 1.0;
        for c in &cols {
            let d: f64 = cand.iter().zip(c).map(|(a, b)| a * b).sum();
            for i in 0..p {
                cand[i] -= d * c[i];
            }
        }
        let nrm = cand.iter().map(|x| x * x).sum::<f64>().sqrt();
        if nrm > 1e-9 {
            for x in &mut cand {
                *x /= nrm;
            }
            cols.push(cand);
            if cols.len() == p {
                break;
            }
        }
    }
    let mut q = vec![0.0f64; p * p];
    for (j, c) in cols.iter().enumerate() {
        for (i, &val) in c.iter().enumerate() {
            q[i * p + j] = val;
        }
    }
    q
}

pub fn glm_lrt(args: &GlmLrtArgs, output: &mut dyn Write) -> Result<u64> {
    let &GlmLrtArgs {
        counts: counts_path,
        design: design_path,
        norm_factors: norm_factors_path,
        coef,
        contrast: contrast_path,
        dispersion: dispersion_arg,
        dispersion_file: dispersion_path,
        fdr,
    } = args;
    let m = Matrix::load(counts_path)?;
    let design = Design::load(design_path)?;
    if design.n_samples != m.n_samples {
        return Err(RsomicsError::InvalidInput(format!(
            "design has {} rows but matrix has {} samples",
            design.n_samples, m.n_samples
        )));
    }
    if design.n_coef < 2 {
        return Err(RsomicsError::InvalidInput(
            "Need at least two columns for design, usually the first is the intercept column"
                .into(),
        ));
    }
    let not_estimable = non_estimable(&design);
    if !not_estimable.is_empty() {
        return Err(RsomicsError::InvalidInput(format!(
            "Design matrix not of full rank. The following coefficients not estimable: {}",
            not_estimable.join(" ")
        )));
    }

    let norm_factors = match norm_factors_path {
        Some(p) => load_norm_factors(p, m.n_samples)?,
        None => vec![1.0; m.n_samples],
    };

    let mut lib = vec![0.0f64; m.n_samples];
    for row in m.counts.chunks_exact(m.n_samples) {
        for (s, &c) in lib.iter_mut().zip(row) {
            *s += c;
        }
    }
    let eff_lib: Vec<f64> = lib
        .iter()
        .zip(&norm_factors)
        .map(|(&l, &f)| l * f)
        .collect();
    let offset: Vec<f64> = eff_lib.iter().map(|&l| l.ln()).collect();
    if !offset.iter().all(|o| o.is_finite()) {
        return Err(RsomicsError::InvalidInput(
            "offsets must be finite values (a library size is zero)".into(),
        ));
    }

    // edgeR glmFit reports coefficients (hence logFC) from a prior.count=0.125
    // augmented fit, while the deviance and LR come from the un-augmented fit.
    let mean_eff = eff_lib.iter().sum::<f64>() / eff_lib.len() as f64;
    let prior: Vec<f64> = eff_lib
        .iter()
        .map(|&l| PRIOR_COUNT * l / mean_eff)
        .collect();
    let offset_aug: Vec<f64> = eff_lib
        .iter()
        .zip(&prior)
        .map(|(&l, &p)| (l + 2.0 * p).ln())
        .collect();

    let test = match (coef, contrast_path) {
        (Some(_), Some(_)) => {
            return Err(RsomicsError::InvalidInput(
                "give --coef or --contrast, not both".into(),
            ));
        }
        (Some(c), None) => {
            if c == 0 || c > design.n_coef {
                return Err(RsomicsError::InvalidInput(format!(
                    "--coef {c} out of range 1..={}",
                    design.n_coef
                )));
            }
            Test::Coef(c - 1)
        }
        (None, Some(p)) => Test::Contrast(load_contrast(p, design.n_coef)?),
        (None, None) => Test::Coef(design.n_coef - 1),
    };

    let dispersions = match dispersion_path {
        Some(p) => {
            let v = load_dispersions(p, m.n_genes())?;
            Dispersion::PerGene(v)
        }
        None => Dispersion::Common(dispersion_arg),
    };

    // Reparametrize once: for a coef test we fit the original design and a
    // column-dropped reduced one; for a contrast we rotate so the tested
    // direction is the first coefficient and report its log2FC.
    let (full_design, reduced_design, fc_coef, df) = match &test {
        Test::Coef(c) => (
            DesignRef::Borrowed(&design),
            reduced_design_coef(&design, *c),
            *c,
            1usize,
        ),
        Test::Contrast(c) => {
            let (full, reduced) = contrast_designs(&design, c);
            (DesignRef::Owned(full), reduced, 0usize, 1usize)
        }
    };
    let full = full_design.as_ref();

    // For a contrast, the logFC is the contrast's first-coefficient estimate
    // rescaled, but the natural contrast value is c·β on the ORIGINAL design.
    let contrast_norm = match &test {
        Test::Contrast(c) => c.iter().map(|x| x * x).sum::<f64>().sqrt(),
        Test::Coef(_) => 1.0,
    };

    let n = m.n_samples;
    let start_full = start_direction(full);
    let start_reduced = start_direction(&reduced_design);

    let per_gene = |w: &mut GeneScratch, g: usize| -> (f64, f64, f64, f64) {
        let row = m.row(g);
        let disp = match &dispersions {
            Dispersion::Common(d) => *d,
            Dispersion::PerGene(v) => v[g],
        };
        let (_b, dev_full) = fit_nb_glm(row, full, &offset, disp, false, &start_full, &mut w.full);
        let (_b0, dev_null) = fit_nb_glm(
            row,
            &reduced_design,
            &offset,
            disp,
            false,
            &start_reduced,
            &mut w.reduced,
        );

        for (s, (&c, &p)) in row.iter().zip(&prior).enumerate() {
            w.row_aug[s] = c + p;
        }
        let (beta_aug, _) = fit_nb_glm(
            &w.row_aug,
            full,
            &offset_aug,
            disp,
            true,
            &start_full,
            &mut w.full,
        );
        let logfc = beta_aug[fc_coef] * contrast_norm / LN2;

        let logcpm = ave_log_cpm_gene(row, &eff_lib);
        let stat = (dev_null - dev_full).max(0.0);
        (logfc, logcpm, stat, special::pchisq_upper(stat, df as f64))
    };

    let make = || GeneScratch {
        full: Scratch::new(n, full.n_coef),
        reduced: Scratch::new(n, reduced_design.n_coef.max(1)),
        row_aug: vec![0.0; n],
    };

    let rows: Vec<(f64, f64, f64, f64)> = if rayon::current_num_threads() > 1 {
        use rayon::prelude::*;
        (0..m.n_genes())
            .into_par_iter()
            .map_init(make, |w, g| per_gene(w, g))
            .collect()
    } else {
        let mut w = make();
        (0..m.n_genes()).map(|g| per_gene(&mut w, g)).collect()
    };

    let logfc: Vec<f64> = rows.iter().map(|r| r.0).collect();
    let logcpm: Vec<f64> = rows.iter().map(|r| r.1).collect();
    let lr: Vec<f64> = rows.iter().map(|r| r.2).collect();
    let pvals: Vec<f64> = rows.iter().map(|r| r.3).collect();

    let fdr_vals = if fdr { Some(bh_fdr(&pvals)) } else { None };
    let gene_col = m.header.split('\t').next().unwrap_or("gene");
    let mut header = format!("{gene_col}\tlogFC\tlogCPM\tLR\tPValue");
    if fdr {
        header.push_str("\tFDR");
    }
    writeln!(output, "{header}").map_err(RsomicsError::Io)?;
    for g in 0..m.n_genes() {
        write!(
            output,
            "{}\t{:.7}\t{:.6}\t{:.6}\t{:.6e}",
            m.genes[g], logfc[g], logcpm[g], lr[g], pvals[g]
        )
        .map_err(RsomicsError::Io)?;
        if let Some(f) = &fdr_vals {
            write!(output, "\t{:.6e}", f[g]).map_err(RsomicsError::Io)?;
        }
        writeln!(output).map_err(RsomicsError::Io)?;
    }
    Ok(m.n_genes() as u64)
}

fn load_dispersions(path: &Path, n_genes: usize) -> Result<Vec<f64>> {
    let file = File::open(path)
        .map_err(|e| RsomicsError::InvalidInput(format!("{}: {e}", path.display())))?;
    let mut v = Vec::with_capacity(n_genes);
    for line in BufReader::new(file).lines() {
        let line = line.map_err(RsomicsError::Io)?;
        let line = line.trim();
        if line.is_empty() || line.starts_with('#') {
            continue;
        }
        let val = line.rsplit('\t').next().unwrap_or(line);
        v.push(
            val.parse::<f64>().map_err(|_| {
                RsomicsError::InvalidInput(format!("non-numeric dispersion '{val}'"))
            })?,
        );
    }
    if v.len() != n_genes {
        return Err(RsomicsError::InvalidInput(format!(
            "{} dispersions for {n_genes} genes",
            v.len()
        )));
    }
    Ok(v)
}

enum DesignRef<'a> {
    Borrowed(&'a Design),
    Owned(Design),
}
impl<'a> DesignRef<'a> {
    fn as_ref(&self) -> &Design {
        match self {
            DesignRef::Borrowed(d) => d,
            DesignRef::Owned(d) => d,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn deviance_zero_at_mle() {
        let y = [10.0, 12.0, 8.0];
        let mu = [10.0, 12.0, 8.0];
        assert!(nb_deviance(&y, &mu, 0.1).abs() < 1e-12);
    }

    #[test]
    fn solve_identity() {
        let mut a = vec![2.0, 0.0, 0.0, 3.0];
        let mut b = [4.0, 9.0];
        let mut x = vec![0.0; 2];
        assert!(solve(&mut a, &mut b, &mut x, 2));
        assert!((x[0] - 2.0).abs() < 1e-12 && (x[1] - 3.0).abs() < 1e-12);
    }

    #[test]
    fn bh_monotone() {
        let p = [0.01, 0.04, 0.03, 0.2];
        let adj = bh_fdr(&p);
        for w in adj.windows(2) {
            let _ = w;
        }
        assert!(adj.iter().all(|&v| (0.0..=1.0).contains(&v)));
    }

    fn design(cols: &[&str], rows: &[&[f64]]) -> Design {
        let n_coef = cols.len();
        let mut data = Vec::new();
        for r in rows {
            data.extend_from_slice(r);
        }
        Design {
            data,
            n_samples: rows.len(),
            n_coef,
            coef_names: cols.iter().map(|s| s.to_string()).collect(),
        }
    }

    #[test]
    fn full_rank_design_estimable() {
        let d = design(
            &["Intercept", "groupb"],
            &[
                &[1.0, 0.0],
                &[1.0, 0.0],
                &[1.0, 0.0],
                &[1.0, 1.0],
                &[1.0, 1.0],
                &[1.0, 1.0],
            ],
        );
        assert!(non_estimable(&d).is_empty());
    }

    #[test]
    fn redundant_column_flags_later_duplicate() {
        // dup is a copy of Intercept; R's qr pivoting flags the later column.
        let d = design(
            &["Intercept", "groupb", "dup"],
            &[
                &[1.0, 0.0, 1.0],
                &[1.0, 0.0, 1.0],
                &[1.0, 0.0, 1.0],
                &[1.0, 1.0, 1.0],
                &[1.0, 1.0, 1.0],
                &[1.0, 1.0, 1.0],
            ],
        );
        assert_eq!(non_estimable(&d), vec!["dup".to_string()]);
    }
}
