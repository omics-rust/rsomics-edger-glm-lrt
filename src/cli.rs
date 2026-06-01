use std::path::PathBuf;

use clap::Parser;
use rsomics_common::{CommonFlags, Result, RsomicsError, Tool, ToolMeta};
use rsomics_help::{Example, FlagSpec, HelpSpec, Section};

use rsomics_edger_glm_lrt::{GlmLrtArgs, glm_lrt};

pub const META: ToolMeta = ToolMeta {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
};

#[derive(Parser, Debug)]
#[command(name = "rsomics-edger-glm-lrt", version, about, long_about = None, disable_help_flag = true)]
pub struct Cli {
    pub counts: PathBuf,
    #[arg(long, value_name = "PATH")]
    design: PathBuf,
    #[arg(long, value_name = "N")]
    coef: Option<usize>,
    #[arg(long, value_name = "PATH")]
    contrast: Option<PathBuf>,
    #[arg(long, default_value_t = 0.05)]
    dispersion: f64,
    #[arg(long, value_name = "PATH")]
    dispersion_file: Option<PathBuf>,
    #[arg(long)]
    norm_factors: Option<PathBuf>,
    #[arg(long)]
    fdr: bool,
    #[arg(short = 'o', long, default_value = "-")]
    output: String,
    #[command(flatten)]
    pub common: CommonFlags,
}

impl Tool for Cli {
    fn meta() -> ToolMeta {
        META
    }
    fn common(&self) -> &CommonFlags {
        &self.common
    }

    fn execute(self) -> Result<()> {
        self.common.install_rayon_pool()?;
        let mut out: Box<dyn std::io::Write> = if self.output == "-" {
            Box::new(std::io::stdout().lock())
        } else {
            Box::new(std::fs::File::create(&self.output).map_err(RsomicsError::Io)?)
        };
        let n = glm_lrt(
            &GlmLrtArgs {
                counts: &self.counts,
                design: &self.design,
                norm_factors: self.norm_factors.as_deref(),
                coef: self.coef,
                contrast: self.contrast.as_deref(),
                dispersion: self.dispersion,
                dispersion_file: self.dispersion_file.as_deref(),
                fdr: self.fdr,
            },
            &mut out,
        )?;
        if !self.common.quiet {
            eprintln!("{n} genes tested");
        }
        Ok(())
    }
}

pub static HELP: HelpSpec = HelpSpec {
    name: env!("CARGO_PKG_NAME"),
    version: env!("CARGO_PKG_VERSION"),
    tagline: "Negative-binomial GLM fit + likelihood-ratio test of a coefficient or contrast (edgeR glmFit/glmLRT).",
    origin: None,
    usage_lines: &[
        "<counts.tsv> --design design.tsv [--coef N | --contrast c.tsv] [--dispersion D | --dispersion-file f.tsv] [--norm-factors f.tsv] [--fdr] [-o de.tsv]",
    ],
    sections: &[Section {
        title: "OPTIONS",
        flags: &[
            FlagSpec {
                short: None,
                long: "design",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: true,
                default: None,
                description: "Design matrix TSV: a header of coefficient names then one numeric row per sample.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "coef",
                aliases: &[],
                value: Some("<N>"),
                type_hint: Some("usize"),
                required: false,
                default: Some("last coefficient"),
                description: "1-based design column to test by dropping it from the model.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "contrast",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Contrast vector (one weight per design coefficient) to test instead of --coef.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "dispersion",
                aliases: &[],
                value: Some("<float>"),
                type_hint: Some("f64"),
                required: false,
                default: Some("0.05"),
                description: "Common negative-binomial dispersion shared across genes.",
                why_default: Some("edgeR's fallback when no dispersion is estimated."),
            },
            FlagSpec {
                short: None,
                long: "dispersion-file",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Per-gene dispersions (one per row, gene order), overriding --dispersion.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "norm-factors",
                aliases: &[],
                value: Some("<path>"),
                type_hint: Some("PathBuf"),
                required: false,
                default: None,
                description: "Per-sample normalization factors (TMM etc.); multiplied into library sizes.",
                why_default: None,
            },
            FlagSpec {
                short: None,
                long: "fdr",
                aliases: &[],
                value: None,
                type_hint: Some("flag"),
                required: false,
                default: None,
                description: "Append a Benjamini-Hochberg FDR column.",
                why_default: None,
            },
        ],
    }],
    examples: &[
        Example {
            description: "Test the last design coefficient, common dispersion 0.1",
            command: "rsomics-edger-glm-lrt counts.tsv --design design.tsv --dispersion 0.1 -o de.tsv",
        },
        Example {
            description: "Test coefficient 2 with TMM factors and a BH-FDR column",
            command: "rsomics-edger-glm-lrt counts.tsv --design design.tsv --coef 2 --norm-factors tmm.tsv --fdr -o de.tsv",
        },
        Example {
            description: "Test a contrast vector",
            command: "rsomics-edger-glm-lrt counts.tsv --design design.tsv --contrast contrast.tsv -o de.tsv",
        },
    ],
    json_result_schema_doc: None,
};

#[cfg(test)]
mod tests {
    use super::*;
    use clap::CommandFactory;

    #[test]
    fn cli_debug_assert() {
        Cli::command().debug_assert();
    }
}
