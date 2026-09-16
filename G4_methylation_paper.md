G-richness adjustment attenuates associations between predicted
G-quadruplexes and age-related DNA methylation
================
CCRRBL-20 project team
2026-09-11

# Abstract

G-quadruplexes (G4s) can influence DNA methylation, but their enrichment
in guanine-rich regulatory sequences complicates tests of an association
with methylation aging. We examined whether Quadron-predicted G4
neighborhoods were associated with signed methylation-versus-age slopes
and whether these associations persisted after local G-richness
adjustment. This exploratory analysis used 137 human tissue specimens
from adipose (32), liver (79), and skeletal muscle (26). Sex- and
body-mass-index-adjusted CpG age slopes were regressed on overlapping
stable- and unstable-labelled Quadron indicators and mean methylation,
first without and then with two strand-oriented G fractions. Before
G-richness adjustment, both G4 indicators were positively associated
with signed age slopes in liver and muscle, but not detectably in
adipose. Adjustment reduced or reversed the coefficients. Only the
stable-labelled muscle term retained a Benjamini-Hochberg-adjusted
$q<0.05$: $1.97 \times 10^{-4}$ beta units per decade (approximate 95%
interval $6.19 \times 10^{-5}$ to $3.32 \times 10^{-4}$; $P=0.004$,
$q=0.013$). Direct stable-minus-unstable contrasts were non-significant
in every tissue under both models. Thus, much of the numerical G4
association is sensitive to local sequence composition. The residual
muscle association is small and does not demonstrate a stability-class
difference or a causal G4 effect. Shared subject identifiers, unresolved
donor identity, probe dependence, and incomplete annotation provenance
limit the strength of inference.

**Keywords:** G-quadruplex; Quadron; DNA methylation; aging; guanine
content; skeletal muscle

# Introduction

DNA methylation changes with age, but these changes are not distributed
uniformly across the genome. Age-associated hypermethylation is enriched
at particular regulatory compartments, including bivalent promoters,
while other regions lose methylation ([<span class="nocase">Rakyan et
al.</span> 2010](#ref-rakyan2010)). The reproducibility of some CpG
changes enables methylation-based age prediction ([Horvath
2013](#ref-horvath2013)), but predictive accuracy does not establish
which local molecular processes generate those changes. Identifying
genomic features associated with methylation-versus-age slopes may help
distinguish the sequence contexts in which methylation is maintained
from those in which it is remodeled.

G-quadruplexes are one candidate feature. These structures arise from
stacked guanine tetrads and occur in sequences with the capacity to fold
outside the canonical DNA duplex. G4-seq established a broad landscape
of G4-forming potential in the human genome ([Chambers et al.
2015](#ref-chambers2015)), whereas chromatin mapping showed that
detectable cellular G4s are concentrated in accessible regulatory
regions and are influenced by transcriptional state
([<span class="nocase">Hänsel-Hertsch et al.</span>
2016](#ref-hansel2016)). Consequently, a predicted G4 motif and an
occupied G4 structure in a particular tissue are related but different
biological quantities.

There is a specific mechanistic basis for connecting G4s to methylation.
Mao and colleagues linked G4 formation to CpG-island hypomethylation,
demonstrated high-affinity DNMT1 binding to G4 structures, and showed
that G4 folding could inhibit DNMT1 activity ([<span class="nocase">Mao
et al.</span> 2018](#ref-mao2018)). These findings establish a route by
which DNA structure can affect methylation maintenance. They do not,
however, predict a universal direction of age-associated change: a low
mean methylation level, a positive age slope, and an increased absolute
slope are distinct observations.

Sequence composition is central to this distinction. G4-forming
potential depends on guanine arrangement within a motif and its
surroundings; the Quadron predictor explicitly learns sequence features
related to experimental folding potential ([Sahakyan et al.
2017](#ref-sahakyan2017)). The same G-rich surroundings may also mark
methylation-relevant regulatory architecture. An association between
predicted G4s and age slopes could therefore arise from structure,
correlated sequence context, or both.

Here, we ask whether associations between Quadron-labelled neighborhoods
and signed methylation age slopes persist after adding local
strand-oriented G fractions to a mean-methylation-adjusted model. We
analyze adipose, liver, and skeletal muscle separately using public 450K
methylation data from an obesity-related tissue study
([<span class="nocase">Horvath et al.</span> 2014](#ref-horvath2014)).
Two questions are kept distinct: whether each G4 annotation is
associated with an age-slope shift, and whether the stable- and
unstable-labelled annotations differ directly. This comparison tests
sensitivity to measured G-richness; it is not a genome-wide test of G4
causality or an epigenetic-clock analysis.

# Materials and methods

## Study material and analytical population

Processed Illumina HumanMethylation450 beta matrices and prepared sample
metadata were drawn from GSE61257 (adipose), GSE61258 (liver), and
GSE61259 (muscle), within the GSE61256 study family
([<span class="nocase">Horvath et al.</span> 2014](#ref-horvath2014)).
The prepared analysis inputs contained 32, 79, and 26 specimens,
respectively (Table 1). All had complete age, sex, and BMI metadata. The
available tissues were not matched population samples: in particular,
the muscle series comprised specimens with BMI between 36.9 and 74.9
kg/m². No claim of healthy-population representativeness is made.

Some subject identifiers occur in more than one tissue series. An
identifier also recurs with inconsistent demographics within the adipose
metadata. We therefore report specimen counts rather than a verified
count of unique participants, analyze tissues separately, and do not
interpret the three series as independent replication cohorts. The
prepared beta matrices, rather than newly reprocessed IDATs, were used;
the retained minimal-model outputs do not establish a complete upstream
normalization or probe-QC history.

**Table 1. Prepared tissue specimens and probe sets used in the two
models.** Counts are metadata-eligible specimens, not independent
donors. Individual CpG fits additionally required finite beta
observations. Both second-stage models used the same number of CpGs
within each tissue.

| Tissue          | GEO      | Specimens | Age (years) | Female/male | BMI (kg/m²) | CpGs    |
|:----------------|:---------|----------:|:------------|:------------|:------------|:--------|
| Adipose         | GSE61257 |        32 | 31-79       | 23/9        | 21.4-74.9   | 381,175 |
| Liver           | GSE61258 |        79 | 21-86       | 34/45       | 14.8-66.9   | 381,179 |
| Skeletal muscle | GSE61259 |        26 | 31-64       | 20/6        | 36.9-74.9   | 381,175 |

## CpG-specific methylation age slopes

Within each tissue, beta-matrix columns were matched to metadata by
sample identifier. For CpG $i$, specimen $j$, and tissue $t$, ordinary
least squares estimated

$$
\beta_{ijt}=\alpha_{it}+b_{it}\operatorname{Age}_{jt}
+\gamma_{it}\operatorname{Sex}_{jt}+\eta_{it}\operatorname{BMI}_{jt}+\epsilon_{ijt}.
$$

Each fit required more finite observations than fitted design columns,
ordinarily at least five for the intercept-age-sex-BMI design. The
estimated coefficient $\hat b_{it}$, in beta units per year, was the
second-stage response. Its sign distinguishes a more positive from a
more negative age association; it does not measure stochastic
variability or an individual’s longitudinal change. The probe universe
was restricted by the project’s gene-annotated hg38 manifest and finite
fitted slopes, rather than encompassing every array probe or every
genomic CpG.

Mean methylation, $\bar\beta_{it}$, was calculated for each CpG within
each tissue over the available specimens. This is a same-dataset mean,
not a pre-aging baseline or an independently observed reference. Both
reported annotation models condition on this mean.

## Quadron neighborhood and sequence annotations

The analysis used hg38 Quadron interval files supplied under stable and
unstable labels. Quadron predicts sequence-driven G4-forming propensity
from experimental training data ([Sahakyan et al.
2017](#ref-sahakyan2017)); these input classes are not measurements of
G4 occupancy or thermodynamic stability in the assayed tissues. The
original class-assignment cutoff and input-generation version were not
recoverable from the retained regression artifacts, so labels are used
descriptively rather than assigned an assumed threshold.

Each input interval was padded by 100 bp on both sides. Overlapping or
adjacent intervals were merged, and probes inherited the stable and
unstable annotations present anywhere within their assigned merged
neighborhood. Indicators $S_i$ and $U_i$ could both equal one; they did
not impose mutually exclusive classes. A probe assigned to a mixed
neighborhood need not directly overlap a motif of each class.
Neighborhood annotation was unstranded, and broad merged clusters were
retained in the minimal analyses rather than modeled with a separate
density coefficient.

The sequence correction used hg38 reference G and C fractions, expressed
relative to the annotated probe strand. On a plus-strand probe, the
probe-oriented G fraction is reference G and the opposite-strand G
fraction is reference C; the assignments are reversed on a minus-strand
probe. Their sum is local GC fraction. They represent local base
composition and orientation, not transcriptional strand or two
independently assayed molecular strands. The implementation extracts
sequence around the probe midpoint with a default 1,000-bp flank,
ordinarily a 2,001-bp inclusive window and distinct from the 100-bp G4
padding. The retained coefficient files do not record whether that
default was overridden at execution. In addition, the original
BED/manifest coordinate convention and hg19-derived probe-strand
concordance with hg38 were not independently revalidated for this
manuscript.

## Two nested annotation models

The mean-methylation-adjusted model (M0, without G-richness correction)
was

$$
\hat b_{it}=a_t+\theta_{S,t}S_i+\theta_{U,t}U_i+\lambda_t\bar\beta_{it}+e_{it}.
$$

The G-richness-adjusted model (M1) added only the two strand-oriented
fractions:

$$
\hat b_{it}=a_t+\theta_{S,t}S_i+\theta_{U,t}U_i+\lambda_t\bar\beta_{it}
+\delta_{P,t}G^{\mathrm{probe}}_i+\delta_{O,t}G^{\mathrm{opposite}}_i+e_{it}.
$$

Both were unweighted ordinary least-squares regressions across CpGs,
fitted separately by tissue. Neither included a stable-by-unstable
interaction, motif density, long-cluster status, CpG-island category,
chromatin state, RepeatMasker class, or replication timing. Each G4
coefficient estimates an additive annotation-associated difference in
the signed age slope conditional on the other G4 indicator and the
included covariates. It is not the mean age slope of a purified
stable-only or unstable-only region set. Both models are exploratory
retained analyses, without a claim of prospective preregistration.

## Statistical inference and reporting

For each model, the implementation applied Benjamini-Hochberg adjustment
across the three tissues separately for each coefficient term
([Benjamini and Hochberg 1995](#ref-benjamini1995)). Thus, the reported
$q$ values for stable and unstable terms belong to separate three-test
families, not a combined family covering both annotations, both models,
or the history of model exploration. Direct contrasts
$\Delta_{S-U,t}=\theta_{S,t}-\theta_{U,t}$ used the fitted coefficient
covariance matrix, with a separate three-tissue BH family per model.

Effects were multiplied by ten to express annotation-associated slope
differences per decade. Figure whiskers and table intervals are
approximate model-based 95% intervals,
$(\hat\theta\pm1.96\,\mathrm{SE})\times10$, not simultaneous confidence
intervals. Numerical attenuation between M0 and M1 is descriptive; no
causal fraction or formal test of between-model coefficient change was
estimated. First-stage slope uncertainty, spatial CpG correlation, and
shared-specimen dependence were not propagated through the second-stage
regressions. Their P values and intervals therefore describe the
implemented model, not dependence-robust population evidence.

# Results

## Positive G4-associated age-slope shifts in liver and muscle before G-richness adjustment

The retained models analyzed 381,175 CpGs in adipose and muscle and
381,179 in liver. With mean methylation included but no G-richness
terms, both Quadron indicators were positive in all three tissues
(Figure 1; Table 2). Adipose estimates were small and did not meet the
saved $q<0.05$ criterion. Liver estimates were $2.74 \times 10^{-4}$
beta units per decade for the stable indicator and $3.61 \times 10^{-4}$
for the unstable indicator. The largest coefficients occurred in muscle:
$6.09 \times 10^{-4}$ and $5.50 \times 10^{-4}$, respectively. Both
indicators met the saved BH criterion in liver and muscle.

These positive coefficients indicate a shift toward a more positive
signed age slope conditional on mean methylation and the other
indicator. They do not establish that all annotated CpGs gain
methylation, that absolute age-related change is greater, or that one
tissue differs significantly from another.

![Quadron mean-methylation-adjusted stable and unstable annotation
coefficients in adipose, liver, and
muscle](G4_methylation_paper_files/figures/quadron_minimal_stable_unstable_mean_beta_g4_terms.png)

**Figure 1. Quadron annotation-associated shifts in signed methylation
age slopes without G-richness correction.** The original
`quadron_minimal_stable_unstable_mean_beta_g4_terms` figure is
reproduced. M0 contains stable and unstable neighborhood indicators and
within-tissue mean beta; the underlying age slopes were adjusted for sex
and BMI. Green and orange bars represent the stable- and
unstable-labelled indicator coefficients multiplied by ten. Whiskers are
approximate model-based 95% intervals. Stars refer to each term’s
BH-adjusted P across tissues: \* $q<0.05$, \*\* $q<0.01$, \*\*\*
$q<0.001$; ns denotes $q\geq0.05$. They do not test the difference
between the two bars. All neighborhoods are retained; the labels do not
identify experimentally verified structures in these tissues.

## Local G-richness adjustment reduces the associations

Adding the two local G fractions changed the magnitude and, in some
cases, the sign of the G4 coefficients (Figure 2; Table 2). Both adipose
terms became slightly negative and remained non-significant. In liver,
the stable term changed sign to $-3.56 \times 10^{-5}$ beta units per
decade, while the unstable term decreased to $6.56 \times 10^{-5}$.
Neither liver term was distinguishable from zero by its nominal test or
the saved BH criterion.

Muscle coefficients remained positive but decreased by 67.7% for the
stable indicator and 72.1% for the unstable indicator. The adjusted
stable estimate was $1.97 \times 10^{-4}$ beta units per decade (95%
interval $6.19 \times 10^{-5}$ to $3.32 \times 10^{-4}$; $P=0.004$,
$q=0.013$). This corresponds to an annotation-associated difference of
approximately 0.0197 percentage points of methylation per decade, not an
observed within-person change. The unstable estimate was
$1.53 \times 10^{-4}$ ($P=0.024$, $q=0.073$): nominally significant, but
above the saved BH threshold.

![Quadron mean-methylation and strand-G-richness-adjusted annotation
coefficients and direct class
contrasts](G4_methylation_paper_files/figures/quadron_minimal_strand_g_richness_g4_terms.png)

**Figure 2. Quadron annotation-associated shifts after strand-oriented
G-richness adjustment.** The original
`quadron_minimal_strand_g_richness_g4_terms` figure is reproduced. M1 is
M0 plus probe-oriented and opposite-strand G fractions; no chromatin or
G4-density covariate is added. Bar units, intervals, colors, and
term-star definitions are as in Figure 1. The brackets represent direct
stable-minus-unstable contrasts, tested using coefficient covariance and
BH-adjusted separately across tissues. All bracket contrasts are
non-significant. The y-axis scale differs from Figure 1; numerical
comparisons should use Table 2. A significant stable term and a
non-significant unstable term in muscle do not establish a difference
between classes.

**Table 2. G4 annotation coefficients before and after local G-richness
adjustment.** Effects and approximate 95% intervals are in units of
$10^{-4}$ beta per decade. M0 adjusts for mean methylation; M1
additionally adjusts for the two G fractions. The $q$ values use the
saved within-term, three-tissue BH families.

| Tissue | Annotation | M0 effect \[95% interval\] | M0 q | M1 effect \[95% interval\] | M1 q |
|:---|:---|:---|:---|:---|:---|
| adipose | Stable G4 | 0.82 \[-0.13, 1.77\] | 0.090 | -0.30 \[-1.25, 0.65\] | 0.538 |
| adipose | Unstable G4 | 0.87 \[-0.07, 1.80\] | 0.071 | -0.22 \[-1.16, 0.72\] | 0.651 |
| liver | Stable G4 | 2.74 \[1.94, 3.55\] | 4.43e-11 | -0.36 \[-1.16, 0.45\] | 0.538 |
| liver | Unstable G4 | 3.61 \[2.81, 4.41\] | 2.44e-18 | 0.66 \[-0.14, 1.45\] | 0.159 |
| muscle | Stable G4 | 6.09 \[4.74, 7.45\] | 3.53e-18 | 1.97 \[0.62, 3.32\] | 0.013 |
| muscle | Unstable G4 | 5.50 \[4.16, 6.84\] | 1.28e-15 | 1.53 \[0.20, 2.87\] | 0.073 |

## Direct comparisons do not distinguish stable and unstable annotations

None of the direct stable-minus-unstable contrasts was significant under
either model (Table 3). After G-richness adjustment, nominal contrast P
values were 0.917 in adipose, 0.130 in liver, and 0.696 in muscle. In
particular, the residual stable-labelled muscle term was not detectably
larger than the unstable-labelled term. The comparison does not
establish equivalence, but it prevents interpreting the different
term-significance labels as evidence for a stability-specific effect.

**Table 3. Direct stable-minus-unstable contrasts within the two
retained models.** Effects and intervals are in units of $10^{-4}$ beta
per decade. Intervals use the saved contrast SE, which incorporates
covariance between the two coefficients. Contrast $q$ values were
adjusted across the three tissues separately within each model.

| Model | Tissue  | Stable minus unstable \[95% interval\] | P     | q     |
|:------|:--------|:---------------------------------------|:------|:------|
| M0    | adipose | -0.04 \[-1.59, 1.51\]                  | 0.957 | 0.957 |
| M0    | liver   | -0.87 \[-2.19, 0.45\]                  | 0.197 | 0.592 |
| M0    | muscle  | 0.60 \[-1.61, 2.81\]                   | 0.596 | 0.894 |
| M1    | adipose | -0.08 \[-1.63, 1.47\]                  | 0.917 | 0.917 |
| M1    | liver   | -1.01 \[-2.32, 0.30\]                  | 0.130 | 0.389 |
| M1    | muscle  | 0.44 \[-1.76, 2.64\]                   | 0.696 | 0.917 |

# Discussion

The main finding is the sensitivity of G4-associated methylation age
slopes to local G-richness. A model containing only the two Quadron
indicators and mean methylation yielded positive liver and muscle
associations. Adding two sequence-composition covariates substantially
reduced those estimates, removed detectable liver associations, and left
a small stable-labelled muscle term under the saved multiple-testing
correction. Direct contrasts did not distinguish stable from unstable
annotations in any tissue. Taken together, the results favor an
interpretation centered on sequence context and a limited residual
association, rather than a general stability-dependent methylation-aging
mechanism.

This sensitivity is biologically plausible without identifying its
cause. G-rich sequence helps determine G4-forming propensity, and
G4-prone regions overlap regulatory chromatin ([Sahakyan et al.
2017](#ref-sahakyan2017); [<span class="nocase">Hänsel-Hertsch et
al.</span> 2016](#ref-hansel2016)). Adjustment can therefore absorb
variation attributable to correlated regulatory features, but it can
also condition on sequence properties integral to the proposed
structural mechanism. The coefficient contraction cannot be assigned
wholly to confounding, nor can it be interpreted as a measured
proportion of methylation aging caused by base composition. The models
show which associations survive a particular linear conditioning set.

The residual muscle result requires equally careful interpretation. Its
effect is small on the beta scale, and the stable-minus-unstable
comparison is non-significant. Thus, the result does not support the
claim that stable structures have an effect absent from unstable
structures. It also does not establish muscle specificity: the tissues
were modeled separately, have different age and BMI distributions, and
include shared subject identifiers. A test of whether coefficients
differ across tissues would require a joint, donor-aware design and was
not performed here. Replication in an independent, better-characterized
muscle population is more informative than treating the three tissue
series as interchangeable confirmations.

The association should also be distinguished from the mechanistic
observations that motivated it. G4-dependent DNMT1 inhibition and
CpG-island hypomethylation concern local methylation control
([<span class="nocase">Mao et al.</span> 2018](#ref-mao2018)). Here, the
response is a cross-sectional signed age slope conditional on mean
methylation. A positive G4 coefficient could indicate a greater tendency
toward methylation gain or less pronounced methylation loss relative to
comparable annotated probes. It does not measure G4 occupancy, DNMT1
activity, stochastic drift, or absolute age-slope magnitude. The present
findings neither demonstrate nor refute the proposed biochemical
mechanism.

Several limitations constrain inference. Hundreds of thousands of CpGs
provide many regression observations but do not replace the small number
of biological specimens. Spatial correlation, shared sample-level
estimation error, and heterogeneous precision of CpG slopes were not
modeled in the second stage. The reported SEs and BH values should
therefore not be read as dependency-robust validation. Moreover, mean
beta is estimated from the same specimens as the slopes and entered
linearly; residual nonlinear dependence on mean methylation may remain.
Neither reported model controls for cell composition, CpG-island
architecture, tissue-matched chromatin, or technical batch. The
obesity-related recruitment and restricted muscle age/BMI range further
limit population generalization.

Reproducibility limitations are distinct from these statistical
limitations. The original Quadron class cutoff, execution-time
sequence-window setting, upstream preprocessing details, and
coordinate/strand convention were not completely recoverable from the
retained artifacts. Merged neighborhoods also represent a regional
annotation, not guaranteed direct proximity to every contributing motif.
The manuscript preserves the two result sets as generated rather than
retrospectively changing their definitions. These gaps need resolution
before a coordinate-exact or stability-specific biological claim is
made.

The next decisive analyses should first resolve sample identity and
freeze the annotation and coordinate definitions. An independent tissue
dataset could then test the same estimands with donor resampling and
genomic-block uncertainty, propagating first-stage slope error and
assessing nonlinear mean-methylation effects. Tissue-matched G4
measurements would separate sequence potential from structural
occupancy. Finally, structure-disrupting changes that preserve local
base composition, coupled to methylation measurements, would offer a
more discriminating causal test than another association model. None of
these validations is claimed to have been performed here.

# Conclusion

In these human tissue datasets, local G-richness adjustment markedly
weakens associations between Quadron-predicted G4 neighborhoods and
signed methylation age slopes. A small residual stable-labelled muscle
association remains under the implemented model, but neither retained
analysis demonstrates a stable-versus-unstable difference. The evidence
supports a sequence-context-sensitive association and a focused
replication target, not a general claim that G4 stability determines the
rate of methylation aging.

# Data and code availability

The source series are
[GSE61257](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61257),
[GSE61258](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61258),
and
[GSE61259](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE61259).
[main.nf](main.nf) is the single computational entry point, with input
and resource settings in [nextflow.config](nextflow.config). Its module
invokes the commented R calculation scripts, including
[bin/analyze_g4_methylation_age.R](bin/analyze_g4_methylation_age.R) for
age slopes and
[bin/compare_g4_effects_adjusted.R](bin/compare_g4_effects_adjusted.R)
for M0/M1. The final Nextflow task compares every column of the six
regenerated model tables with the frozen manuscript references. The
[README.md](README.md) documents prepared-input, saved-slope, and
processed-GEO entry points.

The six retained model tables, aggregate cohort summary, and unchanged
original figure PNGs are bundled alongside this manuscript so that the
Markdown does not depend on an ignored results directory.
[G4_methylation_paper_files/source_manifest.csv](G4_methylation_paper_files/source_manifest.csv)
maps each model/figure artifact to its source and MD5 checksum. The full
[M0
coefficients](G4_methylation_paper_files/data/quadron_minimal_stable_unstable_mean_beta_coefficients.csv)
and [M1
coefficients](G4_methylation_paper_files/data/quadron_minimal_strand_g_richness_coefficients.csv)
retain the exact formulas and unrounded results. Rendering checks the
frozen asset hashes without refreshing them. Re-estimating the
scientific models requires the beta matrices or saved slopes and the
relevant annotation checkpoints; Nextflow retains task commands, logs,
and caching metadata separately from the paper.

From the BTEP directory, the Rmd directly regenerates Markdown and
embeds the figure and math resources into HTML. Rendering requires the
documented R packages and internet access for the KaTeX dependency; the
completed HTML can be read offline:

``` r
rmarkdown::render("G4_methylation_paper.Rmd", output_format = "all")
```

# Ethics statement

This work is a secondary analysis of publicly released processed data.
No new human specimens were collected. Recruitment, consent, and
original study oversight are described by the source study
([<span class="nocase">Horvath et al.</span> 2014](#ref-horvath2014)).

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-benjamini1995" class="csl-entry">

Benjamini, Yoav, and Yosef Hochberg. 1995. “Controlling the False
Discovery Rate: A Practical and Powerful Approach to Multiple Testing.”
*Journal of the Royal Statistical Society: Series B (Methodological)* 57
(1): 289–300. <https://doi.org/10.1111/j.2517-6161.1995.tb02031.x>.

</div>

<div id="ref-chambers2015" class="csl-entry">

Chambers, Vicki S., Giovanni Marsico, Jonathan M. Boutell, Marco Di
Antonio, Geoffrey P. Smith, and Shankar Balasubramanian. 2015.
“High-Throughput Sequencing of DNA G-Quadruplex Structures in the Human
Genome.” *Nature Biotechnology* 33 (8): 877–81.
<https://doi.org/10.1038/nbt.3295>.

</div>

<div id="ref-hansel2016" class="csl-entry">

<span class="nocase">Hänsel-Hertsch, Robert, Dario Beraldi, Stefanie V.
Lensing, et al.</span> 2016. “G-Quadruplex Structures Mark Human
Regulatory Chromatin.” *Nature Genetics* 48 (10): 1267–72.
<https://doi.org/10.1038/ng.3662>.

</div>

<div id="ref-horvath2013" class="csl-entry">

Horvath, Steve. 2013. “DNA Methylation Age of Human Tissues and Cell
Types.” *Genome Biology* 14: R115.
<https://doi.org/10.1186/gb-2013-14-10-r115>.

</div>

<div id="ref-horvath2014" class="csl-entry">

<span class="nocase">Horvath, Steve, Wiebke Erhart, Mario Brosch, et
al.</span> 2014. “Obesity Accelerates Epigenetic Aging of Human Liver.”
*Proceedings of the National Academy of Sciences* 111 (43): 15538–43.
<https://doi.org/10.1073/pnas.1412759111>.

</div>

<div id="ref-mao2018" class="csl-entry">

<span class="nocase">Mao, Shi-Qing, Avazeh T. Ghanbarian, Jochen
Spiegel, et al.</span> 2018. “DNA G-Quadruplex Structures Mold the DNA
Methylome.” *Nature Structural & Molecular Biology* 25 (10): 951–57.
<https://doi.org/10.1038/s41594-018-0131-8>.

</div>

<div id="ref-rakyan2010" class="csl-entry">

<span class="nocase">Rakyan, Vardhman K., Thomas A. Down, Siarhei
Maslau, et al.</span> 2010. “Human Aging-Associated DNA Hypermethylation
Occurs Preferentially at Bivalent Chromatin Domains.” *Genome Research*
20 (4): 434–39. <https://doi.org/10.1101/gr.103101.109>.

</div>

<div id="ref-sahakyan2017" class="csl-entry">

Sahakyan, Aleksandr B., Vicki S. Chambers, Giovanni Marsico, Tobias
Santner, Marco Di Antonio, and Shankar Balasubramanian. 2017. “Machine
Learning Model for Sequence-Driven DNA G-Quadruplex Formation.”
*Scientific Reports* 7: 14535.
<https://doi.org/10.1038/s41598-017-14017-4>.

</div>

</div>
