//! Differential compat against edgeR's glmFit()+glmLRT()+topTags(sort.by="none").
//!
//! `golden_matches_committed` always runs: our binary against the committed
//! R-captured golden. `live_matches_edger` runs the real R upstream when an
//! r-bioc Rscript is found (RSOMICS_RSCRIPT or ~/miniconda3/envs/r-bioc/bin/
//! Rscript), else loud-skips.
//!
//! Values are compared numerically (not byte-compared): R prints `e-01` where
//! Rust prints `e-1`. LR/PValue/FDR are the statistically meaningful columns and
//! match to ~1e-6; logFC matches edgeR to the slack of edgeR's own tol=1e-6
//! Levenberg (ours converges to the augmented MLE, edgeR stops slightly short).

use std::collections::HashMap;
use std::path::{Path, PathBuf};
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-edger-glm-lrt"))
}

fn parse(text: &str) -> (Vec<String>, HashMap<String, Vec<f64>>) {
    let mut lines = text.lines();
    let header: Vec<String> = lines
        .next()
        .unwrap()
        .split('\t')
        .map(str::to_string)
        .collect();
    let mut rows = HashMap::new();
    for line in lines {
        if line.is_empty() {
            continue;
        }
        let mut f = line.split('\t');
        let gene = f.next().unwrap().to_string();
        let vals: Vec<f64> = f.map(|v| v.parse().unwrap()).collect();
        rows.insert(gene, vals);
    }
    (header, rows)
}

fn compare(expected: &str, got: &str, tol: &[(usize, f64)]) {
    let (he, e) = parse(expected);
    let (hg, g) = parse(got);
    assert_eq!(he, hg, "header mismatch");
    assert_eq!(e.len(), g.len(), "row count mismatch");
    for (gene, ev) in &e {
        let gv = g.get(gene).unwrap_or_else(|| panic!("missing gene {gene}"));
        for &(col, t) in tol {
            let a = ev[col];
            let b = gv[col];
            let err = if a.abs() > 1e-3 {
                (a - b).abs() / a.abs()
            } else {
                (a - b).abs()
            };
            assert!(
                err <= t,
                "{gene} col{col}: expected {a} got {b} (err {err:.2e} > {t:.0e})"
            );
        }
    }
}

#[test]
fn golden_matches_committed() {
    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let out = Command::new(bin())
        .arg(golden.join("counts.tsv"))
        .arg("--design")
        .arg(golden.join("design.tsv"))
        .args(["--coef", "2", "--dispersion", "0.1", "--fdr"])
        .output()
        .unwrap();
    assert!(
        out.status.success(),
        "binary failed: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    let got = String::from_utf8(out.stdout).unwrap();
    let expected = std::fs::read_to_string(golden.join("expected_coef2_disp0.1.tsv")).unwrap();
    // logCPM/LR/PValue/FDR exact to ~1e-5; logFC to print precision.
    compare(
        &expected,
        &got,
        &[(0, 1e-3), (1, 1e-5), (2, 1e-4), (3, 1e-5), (4, 1e-5)],
    );
}

/// The three fail-loud guards mirror edgeR 4.4.0 / limma 3.62.1 errors: a
/// single-column design, a rank-deficient design, and a zero-library column all
/// error in R rather than producing output, so ours must too.
fn run_expect_error(counts: &Path, design: &Path, extra: &[&str], needle: &str) {
    let out = Command::new(bin())
        .arg(counts)
        .arg("--design")
        .arg(design)
        .args(extra)
        .output()
        .unwrap();
    assert!(
        !out.status.success(),
        "expected failure but succeeded; stdout:\n{}",
        String::from_utf8_lossy(&out.stdout)
    );
    let stderr = String::from_utf8_lossy(&out.stderr);
    assert!(
        stderr.contains(needle),
        "stderr missing {needle:?}:\n{stderr}"
    );
}

const COUNTS: &str = "gene\tS1\tS2\tS3\tS4\tS5\tS6\n\
g1\t88\t84\t86\t61\t90\t128\n\
g2\t143\t128\t106\t118\t68\t98\n\
g3\t57\t134\t64\t99\t61\t64\n";

#[test]
fn single_column_design_errors() {
    let dir = tempfile::tempdir().unwrap();
    let counts = dir.path().join("counts.tsv");
    let design = dir.path().join("design.tsv");
    std::fs::write(&counts, COUNTS).unwrap();
    std::fs::write(&design, "Intercept\n1\n1\n1\n1\n1\n1\n").unwrap();
    run_expect_error(
        &counts,
        &design,
        &["--dispersion", "0.1"],
        "at least two columns",
    );
}

#[test]
fn rank_deficient_design_errors() {
    let dir = tempfile::tempdir().unwrap();
    let counts = dir.path().join("counts.tsv");
    let design = dir.path().join("design.tsv");
    std::fs::write(&counts, COUNTS).unwrap();
    std::fs::write(
        &design,
        "Intercept\tgroupb\tdup\n1\t0\t1\n1\t0\t1\n1\t0\t1\n1\t1\t1\n1\t1\t1\n1\t1\t1\n",
    )
    .unwrap();
    run_expect_error(
        &counts,
        &design,
        &["--coef", "2", "--dispersion", "0.1"],
        "not of full rank",
    );
}

#[test]
fn zero_library_column_errors() {
    let dir = tempfile::tempdir().unwrap();
    let counts = dir.path().join("counts.tsv");
    let design = dir.path().join("design.tsv");
    std::fs::write(
        &counts,
        "gene\tS1\tS2\tS3\tS4\tS5\tS6\n\
g1\t0\t84\t86\t61\t90\t128\n\
g2\t0\t128\t106\t118\t68\t98\n\
g3\t0\t134\t64\t99\t61\t64\n",
    )
    .unwrap();
    std::fs::write(
        &design,
        "Intercept\tgroupb\n1\t0\n1\t0\n1\t0\n1\t1\n1\t1\n1\t1\n",
    )
    .unwrap();
    run_expect_error(
        &counts,
        &design,
        &["--coef", "2", "--dispersion", "0.1"],
        "offsets must be finite",
    );
}

fn find_rscript() -> Option<PathBuf> {
    if let Ok(p) = std::env::var("RSOMICS_RSCRIPT") {
        let p = PathBuf::from(p);
        if p.exists() {
            return Some(p);
        }
    }
    let home = std::env::var("HOME").ok()?;
    let cand = PathBuf::from(home).join("miniconda3/envs/r-bioc/bin/Rscript");
    cand.exists().then_some(cand)
}

#[test]
fn live_matches_edger() {
    let Some(rscript) = find_rscript() else {
        eprintln!("SKIP live_matches_edger: no r-bioc Rscript (set RSOMICS_RSCRIPT)");
        return;
    };
    let has_edger = Command::new(&rscript)
        .args(["-e", "suppressMessages(library(edgeR))"])
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false);
    if !has_edger {
        eprintln!("SKIP live_matches_edger: Rscript lacks edgeR");
        return;
    }

    let golden = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/golden");
    let counts = golden.join("counts.tsv");
    let design = golden.join("design.tsv");
    let scratch = std::env::temp_dir().join("rsomics-edger-glm-lrt-compat");
    std::fs::create_dir_all(&scratch).unwrap();
    let r_out = scratch.join("r_de.tsv");

    let script = format!(
        r#"
suppressMessages(library(edgeR))
counts <- as.matrix(read.delim("{c}", row.names=1))
design <- as.matrix(read.delim("{d}"))
y <- DGEList(counts=counts)
fit <- glmFit(y, design, dispersion=0.1)
lrt <- glmLRT(fit, coef=2)
tt <- topTags(lrt, n=nrow(counts), sort.by="none")$table
tt$FDR <- p.adjust(tt$PValue, "BH")
out <- data.frame(gene=rownames(tt), logFC=tt$logFC, logCPM=tt$logCPM,
  LR=tt$LR, PValue=tt$PValue, FDR=tt$FDR)
write.table(out, "{o}", sep="\t", quote=FALSE, row.names=FALSE)
"#,
        c = counts.display(),
        d = design.display(),
        o = r_out.display(),
    );
    let st = Command::new(&rscript)
        .args(["-e", &script])
        .status()
        .unwrap();
    assert!(st.success(), "R edgeR run failed");

    let out = Command::new(bin())
        .arg(&counts)
        .arg("--design")
        .arg(&design)
        .args(["--coef", "2", "--dispersion", "0.1", "--fdr"])
        .output()
        .unwrap();
    assert!(out.status.success());
    let got = String::from_utf8(out.stdout).unwrap();
    let expected = std::fs::read_to_string(&r_out).unwrap();
    compare(
        &expected,
        &got,
        &[(0, 2e-3), (1, 1e-5), (2, 1e-4), (3, 1e-5), (4, 1e-5)],
    );
}
