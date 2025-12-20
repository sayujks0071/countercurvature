# Biological Countercurvature Research Report

## Executive Summary
- Unified framework: Normal sagittal curvature and scoliosis-like lateral deviations emerge from a single IEC-Cosserat model operating in different countercurvature regimes.
- Phase diagram: A Riemannian metric derived from information fields yields gravity-dominated, cooperative, and information-dominated regimes, with geodesic deviation as the order parameter.
- Microgravity persistence: Information-driven structure persists as gravitational loading decreases; quantitative deltas require full sweeps (pending).
- Scoliosis bifurcation: Small asymmetries are suppressed in gravity-dominated regimes but amplified into lateral deviations in information-dominated regimes.
- Molecular support: AlphaFold analysis across 36 proteins shows a modest entropy-curvature correlation (r = 0.336, p = 0.0454).

## Abstract
Living systems routinely maintain structure against gravity, from plant stems to vertebrate spines. We propose a quantitative framework of biological countercurvature, where developmental information modifies the effective geometry experienced by a body in a gravitational field. By coupling an Information--Elasticity (IEC) model to Cosserat rod mechanics, we treat the spine as a geodesic in an effective biological spacetime metric $d\ell_{\mathrm{eff}}^2$. We derive a mode selection principle showing that while gravity alone selects a C-shaped sag, the information-coupled system stabilizes an S-shaped counter-curvature mode. A normalized geodesic deviation metric $\widehat{D}_{\mathrm{geo}}$ quantifies this information-driven reshaping. Phase diagrams reveal distinct regimes: gravity-dominated, cooperative, and information-dominated, where the latter predicts the emergence of scoliosis-like lateral deformities as symmetry-broken modes. At the molecular scale, analysis of AlphaFold structures across 36 proteins shows a modest positive correlation between sequence entropy and backbone curvature (r = 0.336, p = 0.0454), providing preliminary support for information-geometry coupling across scales.

## Introduction

### The puzzle of spinal curvature under gravity
Living systems do not simply obey gravity; they negotiate with it. While a passive elastic beam clamped at one end and subject to gravity will sag into a monotonic C-shape, biological structures such as plant stems and vertebrate spines adopt complex, robust geometries that defy this passive tendency. The human spine, in particular, maintains a characteristic S-shaped sagittal profile (cervical and lumbar lordosis, thoracic kyphosis) that is critical for bipedal posture and shock absorption~(white_panjabi_spine, 1990). This shape is not merely a reaction to load but an intrinsic, actively maintained counter-curvature.

### Developmental genetic patterning
The blueprint for this geometry is laid down during embryogenesis. The paraxial mesoderm segments into somites, driven by the segmentation clock and oscillating gene expression (e.g., Notch, Wnt, FGF)~(pourquie2011vertebrate, 2011). These segments acquire distinct identities through the expression of HOX and PAX genes, which specify the morphological characteristics of the resulting vertebrae~(wellik2007hox, 2007). However, the mechanism by which these discrete genetic codes are translated into the continuous, mesoscale geometry of the adult spine remains a fundamental open question.

### Hypothesis: Information--Elasticity Coupling and effective metric
We propose that developmental information acts as a field that modifies the effective geometry experienced by the spine. Drawing an analogy to General Relativity, where matter curves spacetime, we suggest that biological information curves the material manifold of the spine. We formalize this as an Information--Elasticity Coupling (IEC) framework, where a genetic information field $I(s)$ modifies the rest curvature, stiffness, and active moments of the structure. In this view, the spine does not fight gravity; it settles into the geodesic of a curved biological metric $d\ell_{\mathrm{eff}}^2$ shaped by information.

### Contribution and overview
In this work, we:
(i) Define the IEC model and a phenomenological biological metric that maps genetic information to geometric distortions.
(ii) Implement this framework in a 3D Cosserat rod simulation (using PyElastica) to model the spine under gravitational loading.
(iii) Demonstrate that the interplay between gravity and information selects specific spinal modes, shifting the ground state from a passive C-shape to an active S-shape.
(iv) Show that in information-dominated regimes, this same mechanism can amplify small asymmetries into pathological, scoliosis-like deformities.

## Theory: Information--Cosserat Model of Spinal Countercurvature

We propose that the robust S-shaped geometry of the spine arises not from passive mechanical equilibrium under gravity, but from an active counter-curvature mechanism driven by developmental information. We formalize this using an Information--Elasticity Coupling (IEC) framework, where a scalar information field $I(s)$ modifies the effective geometry and energetics of a Cosserat rod.

### Geometry and parameterization
Consider a slender rod parameterized by arc-length $s \in [0, L]$. The configuration is defined by a centerline curve $\mathbf{r}(s) \in \mathbb{R}^3$ and a director frame $\{\mathbf{d}_1, \mathbf{d}_2, \mathbf{d}_3\}(s)$ describing the orientation of cross-sections. The rod deforms under a gravitational field $\mathbf{g} = -g \hat{\mathbf{e}}_z$. In the absence of biological regulation, such a rod would sag into a C-shape (kyphosis) or buckle.

### Information field from developmental patterning
We introduce a scalar field $I(s)$ representing the spatial distribution of developmental identity along the axis. Rather than an arbitrary function, $I(s)$ is grounded in the well-characterized HOX gene expression boundaries that define spinal regionalization~(wellik2007hox). Specifically:
- Cervical-Thoracic Transition: Defined by the anterior expression limit of Hoxc6 (associated with T1).
- Thoracic-Lumbar Transition: Defined by the onset of Hoxc10 and Hoxd10 expression (associated with L1)~(burke1995hox).

We model $I(s)$ as a superposition of these collinear expression domains. The resultant field peaks in the cervical and lumbar regions where counter-curvature (lordosis) is required to resist the flexion moment of gravity, while the thoracic region (dominated by Hoxc6) retains a primary kyphotic curvature.

### Biological metric and effective energy
The central hypothesis of IEC is that the information field modifies the effective metric experienced by the rod. We define a biological metric $d\ell_{\mathrm{eff}}^2$ that dilates or contracts the reference manifold based on information content:

$$
d\ell_{\mathrm{eff}}^2 = g_{\mathrm{eff}}(s)\,ds^2 = \exp\left[2\left(\beta_1 \tilde{I}(s) + \beta_2 \frac{\partial \tilde{I}}{\partial s}\right)\right] ds^2,
$$
where $\tilde{I}$ is the normalized information field and $\beta_{1,2}$ are coupling constants. This metric implies that regions of high information density or gradient have a larger effective length, effectively prescribing a target curvature.

The energetics of the rod are governed by an IEC-modified elastic energy functional. Unlike a passive beam with uniform stiffness $B$, the biological rod minimizes an energy where the bending cost is weighted by the information field:

$$
\mathcal{E} = \int_0^L \frac{1}{2} B(s) \left( \kappa(s) - \kappa_{\mathrm{rest}}(s) \right)^2 w(I(s)) \, ds,
$$
where $\kappa(s)$ is the curvature, $\kappa_{\mathrm{rest}}(s) = \kappa_0 + \chi_\kappa \partial_s I$ is the information-dependent rest curvature, and $w(I) = 1 + \chi_E I(s)$ is a stiffness weighting function. This energy penalizes deviations from the information-prescribed shape more heavily in regions of high $I(s)$.

### Cosserat force and moment balance
The equilibrium configuration is found by minimizing the total potential energy (elastic + gravitational). In the language of Cosserat rod theory, this yields the balance of linear and angular momentum. For a static rod subject to gravity $\mathbf{f}_g = \rho A \mathbf{g}$ and IEC-driven active moments, the equations are:

$$
\begin{aligned}
\mathbf{n}'(s) + \mathbf{f}_g &= \mathbf{0},  \\
\mathbf{m}'(s) + \mathbf{r}'(s) \times \mathbf{n}(s) + \mathbf{m}_{\mathrm{info}}'(s) &= \mathbf{0},
\end{aligned}
$$
where $\mathbf{n}$ is the internal force, $\mathbf{m}$ is the internal moment, and $\mathbf{m}_{\mathrm{info}}$ represents the active couple induced by the information field.

### Mode selection and spinal geometry
The interplay between the gravitational potential (favoring a C-shaped sag) and the IEC energy (favoring an S-shape) can be understood as a mode selection problem. In the linearized planar limit, small deflections $y(s)$ from the vertical satisfy an eigenvalue problem of the form:

$$
\mathcal{L}_{\mathrm{IEC}}[y(s)] = \frac{d^2}{ds^2} \left( B_{\mathrm{eff}}(s) \frac{d^2 y}{ds^2} \right) - \frac{d}{ds} \left( N(s) \frac{dy}{ds} \right) = \lambda_n y_n(s),
$$
where $N(s)$ is the axial tension due to gravity. The information field modifies the operator $\mathcal{L}_{\mathrm{IEC}}$ such that the lowest energy mode $\lambda_0$ shifts from a monotonic C-shape (passive buckling) to a higher-order S-shape (counter-curvature). This spectral shift explains the robustness of the spinal curve: the S-shape becomes the energetic ground state of the information-coupled system.

## Methods

### Deterministic IEC Beam Model
To explore the mode selection mechanism (Eq.~\ref{eq:mode_selection}), we discretize the linearized beam equations using a finite difference scheme on a 1D domain $s \in [0, L]$. The rod is divided into $N=100$ segments. The information field $I(s)$ is mapped to local stiffness $E_i$ and rest curvature $\kappa_{i}$ at each node. We solve the resulting boundary value problem (BVP) using a standard shooting method (or sparse matrix solver for the eigenproblem). This allows rapid exploration of the $(\chi_\kappa, g)$ parameter space to identify regions where S-modes become the ground state.

### 3D Cosserat Rod Implementation (PyElastica)
For full 3D simulations, we utilize PyElastica~(pyelastica_zenodo, 2023; gazzola2018forward, 2018), an open-source Python implementation of Cosserat rod theory. The spine is modeled as a Cosserat rod with the following specifications:

- Discretization: The rod is discretized into $n=50$--$100$ elements.
- IEC Coupling: A custom callback updates the local rest curvature vector $\bm{\kappa}^0(s)$ and bending stiffness matrix $\mathbf{B}(s)$ based on the information field $I(s)$.
- Boundary Conditions: The rod is clamped at the base (sacrum) and free at the top (cranium), simulating a cantilever column under gravity. For specific validation cases, clamped-clamped conditions are used.
- Gravitational Loading: Gravity is applied as a uniform body force $\mathbf{f} = \rho A \mathbf{g}$.
- Damping: To find static equilibrium configurations, external damping ($\nu \sim 0.1$--$1.0$) is applied and the system is integrated until kinetic energy dissipates ($v_{\max} < 10^{-6}$ m/s).

### Parameter Sweeps and Mode Classification
We perform systematic parameter sweeps over coupling strength $\chi_\kappa$ (range $[0, 0.1]$) and gravitational acceleration $g$ (range $[0.01, 1.0]$ $g_{\mathrm{Earth}}$). For each simulation, we compute the equilibrium shape and evaluate:
- Geodesic Deviation $\widehat{D}_{\mathrm{geo}}$ (difference between realized shape and gravity-only geodesic).
- Lateral Deviation $S_{\mathrm{lat}}$ (symmetry breaking in the coronal plane).
- Cobb Angle (clinical proxy for scoliotic curves).

Regimes are classified as gravity-dominated ($\widehat{D}_{\mathrm{geo}} < 0.1$), cooperative ($0.1 < \widehat{D}_{\mathrm{geo}} < 0.3$), or information-dominated ($\widehat{D}_{\mathrm{geo}} > 0.3$).

### AlphaFold Protein Structure Analysis
A BCC protein database spans developmental patterning (HOX, PAX), mechanotransduction, segmentation clock components, longevity/stress response, ECM, and transcription factors (69 proteins total). Structures were retrieved from the AlphaFold Protein Structure Database via the public API using fetch_bcc_structures.py, with PDB files stored in alphafold_analysis/predictions and metadata in alphafold_analysis/metadata.

Quality control excluded XML error responses, files smaller than 100 bytes, and files lacking ATOM/HETATM records. For each retained structure, we computed:
- Sequence entropy (Shannon entropy of amino-acid composition).
- Backbone curvature along the C-alpha trace using a sliding window of 7 residues and a centroid-based curvature proxy (inverse mean radius).
- Flexibility index from bend-angle variability.
- Compactness (radius of gyration and end-to-end distance).
- Mechanical proxies (proline/glycine ratio, instability index, GRAVY).

## Data

### AlphaFold dataset status (current)
- Downloaded PDBs: 37
- Valid PDBs after QC: 36
- Invalid PDBs: 1 (Klotho, XML error response)
- Not found in AlphaFold: 12 proteins

#### Coverage by category
| Category | Targets | Downloaded | Valid | Notes |
|----------|---------|------------|-------|-------|
| HOX | 32 | 22 | 22 | 10 not found |
| PAX | 5 | 0 | 0 | Not downloaded |
| Mechanosensitive | 8 | 1 | 1 | Only YAP1 |
| Segmentation | 11 | 9 | 9 | NOTCH1, FGF4 not found |
| Longevity | 5 | 4 | 3 | Klotho invalid, AMPK missing |
| ECM | 4 | 0 | 0 | Not downloaded |
| Transcription | 4 | 1 | 1 | Only MESP2 |

#### Not found in AlphaFold
HOXA2, HOXA3, HOXA11, HOXB7, HOXB8, HOXC9, HOXD3, HOXD9, HOXD11, HOXD13, NOTCH1, FGF4.

### Data artifacts and locations
- AlphaFold summary report: alphafold_analysis/bcc_analysis_report.md
- AlphaFold data (JSON): alphafold_analysis/bcc_analysis_data.json
- AlphaFold correlation plot: ../alphafold_analysis/figures/entropy_curvature_correlation.png
- IEC figures: outputs/figs/fig_iec_discriminators.png
- Manuscript figures: assets/*.png and assets/*.pdf

### Key numbers box (pending full sweeps)
The following metrics require completion once full sweeps are run:
- Microgravity persistence deltas.
- Phase regime anchor points.
- S-curve shape statistics.
- Scoliosis amplification factors.

## Tables

### AlphaFold summary statistics
- Total proteins analyzed: 36
- Mean sequence length: 452.7 aa
- Mean entropy: 4.024 bits
- Mean curvature: 0.1117
- Entropy-curvature correlation: r = 0.336 (p = 0.0454)

### AlphaFold per-protein metrics
| Protein | Length | Entropy | Mean Curv | Flex Index | Rg | Stiffness | Instability |
|---------|--------|---------|-----------|------------|----|-----------|------------|
| HOXA4 | 320 | 3.956 | 0.1036 | 0.358 | 42.33 | 0.016 | 76.11 |
| HOXB9 | 250 | 4.086 | 0.1079 | 0.334 | 39.59 | 0.010 | 58.40 |
| HOXA5 | 270 | 4.023 | 0.1062 | 0.342 | 41.14 | 0.007 | 65.61 |
| HOXA7 | 230 | 4.069 | 0.1083 | 0.332 | 33.94 | 0.005 | 56.05 |
| HOXD10 | 340 | 4.101 | 0.1027 | 0.363 | 40.67 | 0.008 | 71.70 |
| HOXA1 | 335 | 4.093 | 0.0998 | 0.382 | 44.32 | 0.007 | 59.38 |
| SIRT1 | 747 | 4.086 | 0.1158 | 0.345 | 40.25 | 0.009 | 53.88 |
| HOXD12 | 217 | 3.836 | 0.1008 | 0.294 | 43.57 | 0.013 | 60.19 |
| HOXC8 | 242 | 4.105 | 0.1145 | 0.303 | 34.35 | 0.006 | 63.34 |
| HES1 | 280 | 4.032 | 0.1244 | 0.283 | 40.42 | 0.013 | 52.33 |
| HOXA10 | 410 | 3.922 | 0.0989 | 0.383 | 45.63 | 0.012 | 73.58 |
| HOXD8 | 290 | 4.054 | 0.1112 | 0.315 | 38.03 | 0.013 | 70.57 |
| HOXC10 | 342 | 4.065 | 0.1003 | 0.390 | 42.91 | 0.009 | 57.74 |
| HOXC11 | 304 | 4.076 | 0.1029 | 0.386 | 38.43 | 0.009 | 60.04 |
| HES7 | 225 | 3.760 | 0.1258 | 0.236 | 37.82 | 0.019 | 76.22 |
| YAP1 | 504 | 3.994 | 0.1110 | 0.380 | 44.88 | 0.012 | 65.51 |
| FOXO3 | 673 | 3.974 | 0.0961 | 0.513 | 55.01 | 0.008 | 66.77 |
| DLL3 | 587 | 3.880 | 0.1233 | 0.254 | 37.68 | 0.011 | 55.26 |
| WNT3A | 352 | 4.198 | 0.1353 | 0.212 | 24.81 | 0.005 | 50.65 |
| HOXD1 | 328 | 3.954 | 0.1025 | 0.398 | 46.72 | 0.011 | 66.82 |
| MESP2 | 397 | 3.825 | 0.1014 | 0.336 | 43.54 | 0.012 | 73.92 |
| DLL1 | 247 | 4.089 | 0.1215 | 0.250 | 31.63 | 0.007 | 43.96 |
| HOXD4 | 255 | 4.073 | 0.1106 | 0.312 | 36.73 | 0.013 | 76.31 |
| FGF8 | 233 | 4.101 | 0.1264 | 0.222 | 29.41 | 0.005 | 47.80 |
| HOXB1 | 301 | 3.965 | 0.1023 | 0.399 | 41.58 | 0.013 | 69.02 |
| HOXC4 | 264 | 4.053 | 0.1094 | 0.341 | 37.41 | 0.013 | 81.10 |
| WNT5A | 365 | 4.207 | 0.1383 | 0.206 | 25.41 | 0.002 | 38.85 |
| HOXB2 | 356 | 3.906 | 0.1020 | 0.400 | 45.16 | 0.015 | 85.76 |
| HOXC6 | 153 | 4.036 | 0.1256 | 0.261 | 26.21 | 0.001 | 65.58 |
| PGC1A | 671 | 4.028 | 0.1002 | 0.539 | 49.34 | 0.007 | 78.05 |
| NOTCH2 | 2471 | 4.150 | 0.1289 | 0.325 | 58.07 | 0.007 | 45.58 |
| HOXB6 | 224 | 4.084 | 0.1116 | 0.342 | 34.50 | 0.008 | 74.74 |
| NOTCH3 | 2321 | 4.000 | 0.1291 | 0.329 | 57.16 | 0.011 | 54.79 |
| HOXB4 | 251 | 3.976 | 0.1120 | 0.262 | 36.80 | 0.018 | 90.11 |
| HOXB5 | 269 | 4.005 | 0.1058 | 0.388 | 41.73 | 0.007 | 62.08 |
| HOXA9 | 272 | 4.118 | 0.1046 | 0.364 | 36.80 | 0.009 | 46.63 |

### Key numbers box (template)
Use after full sweep extraction:

```
Microgravity persistence:
- At g = 9.81: Dhat_geo = [TBD], passive energy = [TBD]
- At g = 0.01: Dhat_geo = [TBD], passive energy = [TBD]
- Passive energy collapse: [TBD] percent reduction
- Dhat_geo persistence: [TBD] ratio

Phase diagram regimes:
- Gravity-dominated: chi_kappa = [TBD], g = [TBD], Dhat_geo approx [TBD]
- Cooperative: chi_kappa = [TBD], g = [TBD], Dhat_geo approx [TBD]
- Information-dominated: chi_kappa = [TBD], g = [TBD], Dhat_geo > 0.3

Scoliosis amplification:
- Gravity-dominated: S_lat_asym / S_lat_sym approx [TBD]
- Information-dominated: S_lat_asym / S_lat_sym approx [TBD]
```

## Diagrams

### Figure index (available)
- assets/fig_gene_to_geometry.png: Mapping of HOX/PAX domains to information field and IEC landscape.
- assets/fig_mode_spectrum.png: Gravity-selected vs information-selected eigenmodes.
- assets/fig_countercurvature_panelA.pdf: Passive sag vs IEC-stabilized S-curve.
- assets/fig_countercurvature_panelB.pdf: Microgravity persistence panel.
- assets/fig_countercurvature_panelC.pdf: Additional IEC configuration panel.
- assets/fig_countercurvature_panelD.pdf: Additional IEC configuration panel.
- assets/fig_phase_diagram_scoliosis.png: Phase diagram of countercurvature regimes.
- ../alphafold_analysis/figures/entropy_curvature_correlation.png: AlphaFold entropy-curvature correlation.
- ../outputs/figs/fig_iec_discriminators.png: IEC discriminators (node drift, amplitude modulation, helical threshold map).

### Figure generation notes
- Outputs under outputs/figs are generated by spinalmodes figure utilities (see docs/figures.md).
- AlphaFold correlation plot is generated by alphafold_analysis/analyze_bcc_structures.py.

## Results

### Segmentation-derived information field and IEC landscape
We establish the connection between developmental patterning and the mechanical information field. The mapping from discrete genetic domains (HOX boundaries) to a continuous information field $I(s)$ yields peaks in cervical and lumbar regions. Through the biological metric (Eq.~\ref{eq:biological_metric}), these regions encode the target S-shape into the manifold itself.

### Molecular basis: information-curvature coupling in proteins
We analyzed 36 AlphaFold structures spanning HOX genes (22 proteins), segmentation clock components (9 proteins), longevity factors (3 proteins), and representative mechanosensitive/transcription factors (YAP1, MESP2). We calculated local backbone curvature and compared it to sequence entropy.

The analysis shows a modest positive correlation between sequence information density and structural curvature (r = 0.336, p = 0.0454). The effect is weaker than expected for uniform coupling and likely reflects protein-class heterogeneity, domain architecture diversity, and incomplete coverage for several categories (PAX, ECM, mechanosensitive). These results provide preliminary molecular support for information-geometry coupling while motivating category-stratified analyses as the structure set expands.

### Mode spectrum of the IEC beam in gravity
The passive beam's ground state is a monotonic C-shaped sag. With sufficient IEC coupling ($\chi_\kappa > \chi_{\mathrm{crit}}$), the spectrum shifts and the S-shaped mode becomes the lowest-energy configuration. This confirms that the spinal curve is a gravity-selected mode of the information-modified system.

### 3D Cosserat rod S-curve solutions
3D simulations reproduce cervical and lumbar lordosis and show S-curve persistence in microgravity. The S-curve persists even as $g \to 0$, unlike passive sag, consistent with microgravity observations~(green2018spinal, 2018).

### Phase diagrams and pathology-like modes
Parameter sweeps over $(\chi_\kappa, g)$ reveal gravity-dominated, cooperative, and information-dominated regimes. In the information-dominated regime, small asymmetries (e.g., 5 percent left/right information differences) are amplified into lateral deformities resembling scoliosis.

## Discussion

### Interpretation
The adult spinal shape can be viewed as a standing wave of counter-curvature maintained by developmental information against gravity. The IEC framework provides a quantitative language for this by treating the information field $I(s)$ as a metric modifier that stabilizes the S-curve.

### Links to developmental genetics and evolution
Peaks in the phenomenological $I(s)$ correspond to HOX-coded cervical and lumbar regions, suggesting specific paralogs act as curvature generators. Evolutionarily, bipedalism likely involved tuning this information field to stabilize the S-mode against increased gravitational moments.

### Molecular-scale feedback hypothesis
Across 36 proteins, sequence entropy and backbone curvature are positively correlated (r = 0.336, p = 0.0454), indicating a modest information-geometry signal. We hypothesize that mechanosensitive proteins with curved domains act as geometric sensors that couple tissue-scale mechanics to gene-regulatory feedback. Additional category coverage will clarify whether mechanotransducers show stronger coupling than patterning proteins.

### Comparison to passive pre-stress models
Passive tensegrity or muscle pre-stress models require continuous energy expenditure. In contrast, the IEC framework posits a target shape intrinsic to the material manifold. Muscles then fine-tune deviations rather than enforce an S-shape against a passive baseline.

### Limitations
- The information field is treated as static and phenomenological; real developmental dynamics are time-dependent.
- The spine is modeled as a continuous rod, omitting vertebral and disc anatomy.
- AlphaFold coverage is incomplete (PAX, ECM, mechanosensitive under-sampled), and Klotho was excluded due to an invalid PDB response.

### Future directions
- Patient-specific inference of $I(s)$ from imaging.
- Coupling IEC to growth laws for developmental time-course modeling.
- Integration of sensory feedback (proprioception) as a dynamic information component.

## Conclusion
We have presented a theoretical framework for Biological Countercurvature, positing that spinal geometry is determined by the interaction between a developmental information field and gravity. By coupling a Cosserat rod model with an Information--Elasticity mechanism, the spinal S-curve emerges as a stable equilibrium accessible only when information modifies the effective metric. This framework unifies normal development, microgravity adaptation, and pathological deformity under a single information-geometric mechanism. Molecular-scale analysis provides preliminary support for information-geometry coupling but remains incomplete until additional protein categories are incorporated.

## Legacy and Working Drafts (kept, not deleted)
The following files are retained for provenance but are not merged into the research report:
- `legacy/draft_intro_BCC.md`
- `legacy/draft_methods_BCC.md`
- `legacy/draft_results_BCC.md`
- `legacy/draft_discussion_BCC.md`
- `legacy/draft_conclusion_BCC.md`
- `legacy/draft_theory_BCC.md`
- `legacy/MANUSCRIPT_DRAFT_BCC.md`
- `legacy/biological_countercurvature_manuscript.txt`

## Citations (BibTeX)

```
% Bibliography for Biological Countercurvature of Spacetime
% Minimal but coherent seed bibliography
% Extend as needed with additional references

% ---------------------------------------------------------------------------
% This Paper (Placeholder)
% ---------------------------------------------------------------------------

@article{krishnan2025biological_countercurvature,
  title   = {Biological Countercurvature of Spacetime: An Information--Cosserat Framework for Spinal Geometry},
  author  = {Krishnan, Sayuj and Coauthor, Firstname},
  journal = {To be submitted},
  year    = {2025},
  note    = {Preprint in preparation}
}

% ---------------------------------------------------------------------------
% PyElastica and Cosserat Rod Mechanics
% ---------------------------------------------------------------------------

@software{pyelastica_zenodo,
  author    = {Tekinalp, Arman and Gazzola, Mattia and others},
  title     = {PyElastica: Open-source software for the simulation of assemblies
               of slender, one-dimensional structures using Cosserat rod theory},
  year      = {2023},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.7658872},
  url       = {https://doi.org/10.5281/zenodo.7658872}
}

@article{gazzola2018forward,
  title   = {Forward and inverse problems in the mechanics of soft filaments},
  author  = {Gazzola, Mattia and Dudte, Levi H. and McCormick, Andrew G.
             and Mahadevan, L.},
  journal = {Royal Society Open Science},
  volume  = {5},
  number  = {6},
  pages   = {171628},
  year    = {2018},
  doi     = {10.1098/rsos.171628}
}

@article{zhang2019modeling,
  title   = {Modeling and simulation of complex dynamic musculoskeletal architectures},
  author  = {Zhang, Xiaotian and Chan, Fan K. and Parthasarathy, Tejaswin
             and Gazzola, Mattia},
  journal = {Nature Communications},
  volume  = {10},
  number  = {1},
  pages   = {1--12},
  year    = {2019},
  doi     = {10.1038/s41467-019-12759-5}
}

@article{naughton2021elastica,
  title   = {Elastica: A compliant mechanics environment for soft robotic control},
  author  = {Naughton, Noel and Sun, Jian and Tekinalp, Arman and
             Chowdhary, Girish and Gazzola, Mattia},
  journal = {IEEE Robotics and Automation Letters},
  volume  = {6},
  number  = {2},
  pages   = {3389--3396},
  year    = {2021}
}

@book{antman2005nonlinear,
  title     = {Nonlinear Problems of Elasticity},
  author    = {Antman, Stuart S.},
  edition   = {2},
  year      = {2005},
  publisher = {Springer},
  address   = {New York}
}

@misc{cosseratrods_site,
  title        = {Elastica and Cosserat Rod Theory Case Studies},
  howpublished = {\url{https://www.cosseratrods.org/cosserat_rods/case-studies/}},
  note         = {Accessed 2025-11-16}
}

% ---------------------------------------------------------------------------
% Riemannian Geometry and General Relativity
% ---------------------------------------------------------------------------

@book{lee2018riemannian,
  title     = {Introduction to Riemannian Manifolds},
  author    = {Lee, John M.},
  edition   = {2},
  year      = {2018},
  publisher = {Springer},
  address   = {New York}
}

@article{einstein1916grundlage,
  title   = {Die Grundlage der allgemeinen Relativit{\"a}tstheorie},
  author  = {Einstein, Albert},
  journal = {Annalen der Physik},
  volume  = {49},
  pages   = {769--822},
  year    = {1916},
  doi     = {10.1002/andp.19163540702}
}

@book{wald1984gr,
  title     = {General Relativity},
  author    = {Wald, Robert M.},
  year      = {1984},
  publisher = {University of Chicago Press},
  address   = {Chicago}
}

% ---------------------------------------------------------------------------
% Spinal Biomechanics and Scoliosis
% ---------------------------------------------------------------------------

@book{white_panjabi_spine,
  title     = {Clinical Biomechanics of the Spine},
  author    = {White, Augustus A. and Panjabi, Manohar P.},
  edition   = {2},
  year      = {1990},
  publisher = {Lippincott}
}

@article{weinstein2008adolescent,
  title   = {Adolescent idiopathic scoliosis},
  author  = {Weinstein, Stuart L. and Dolan, Lori A. and Cheng, Jack C. Y.
             and Danielsson, Aina and Morcuende, Jose A.},
  journal = {The Lancet},
  volume  = {371},
  number  = {9623},
  pages   = {1527--1537},
  year    = {2008},
  doi     = {10.1016/S0140-6736(08)60658-3}
}

@article{deBrito2014sitting,
  title   = {Ability to sit and rise from the floor as a predictor of all-cause mortality},
  author  = {de Brito, Leonardo Barbosa and Ricardo, Denise Rodrigues and
             Ara{\'u}jo, Danielle da Silva and Ramos, Plinio Santos and
             Myers, Jonathan and Ara{\'u}jo, Claudio Gil Soares de},
  journal = {European Journal of Preventive Cardiology},
  volume  = {21},
  number  = {7},
  pages   = {892--898},
  year    = {2014},
  doi     = {10.1177/2047487312471759},
  note    = {First published online December 2012}
}

% ---------------------------------------------------------------------------
% Developmental Biology and Genetics
% ---------------------------------------------------------------------------

@article{grimes2016zebrafish,
  title   = {Zebrafish model of idiopathic scoliosis link cerebrospinal fluid flow to defects in spine curvature},
  author  = {Grimes, D. T. and Boswell, C. W. and Morante, N. F. C. and others},
  journal = {Science},
  volume  = {352},
  number  = {6291},
  pages   = {1341--1344},
  year    = {2016},
  doi     = {10.1126/science.aaf6419},
  note    = {Ciliary flow and left-right asymmetry in spinal development}
}

@article{wellik2007hox,
  title   = {Hox genes and regional patterning of the vertebrate body plan},
  author  = {Wellik, D. M.},
  journal = {Developmental Biology},
  volume  = {306},
  number  = {2},
  pages   = {359--372},
  year    = {2007},
  doi     = {10.1016/j.ydbio.2007.02.028},
  note    = {HOX gene expression and segmental identity}
}

@article{pourquie2011vertebrate,
  title   = {Vertebrate segmentation: from cyclic gene networks to scoliosis},
  author  = {Pourqui{\'e}, Olivier},
  journal = {Cell},
  volume  = {145},
  number  = {5},
  pages   = {650--663},
  year    = {2011},
  doi     = {10.1016/j.cell.2011.05.011},
  note    = {Segmentation clock and somitogenesis in vertebrate development}
}

% ---------------------------------------------------------------------------
% Microgravity and Spinal Health
% ---------------------------------------------------------------------------

@article{green2018spinal,
  title   = {Spinal Health during Unloading and Reloading Associated with Spaceflight},
  author  = {Green, David A. and Scott, Jonathan P. R.},
  journal = {Frontiers in Physiology},
  volume  = {8},
  pages   = {1126},
  year    = {2018},
  doi     = {10.3389/fphys.2017.01126}
}

@article{marfia2023microgravity,
  title   = {Microgravity and the intervertebral disc: The impact of spaceflight and
             simulated microgravity on disc degeneration},
  author  = {Marfia, Giacomo and others},
  journal = {Frontiers in Physiology},
  volume  = {14},
  pages   = {1124991},
  year    = {2023},
  doi     = {10.3389/fphys.2023.1124991}
}

% ---------------------------------------------------------------------------
% Counterbend Mechanics (PNAS)
% ---------------------------------------------------------------------------

@article{gadelha2013counterbend,
  title   = {Counterbend phenomena in flagellar mechanics},
  author  = {Gadelha, Hermes and others},
  journal = {Proceedings of the National Academy of Sciences},
  volume  = {110},
  number  = {30},
  pages   = {12180--12185},
  year    = {2013},
  doi     = {10.1073/pnas.1306988110},
  note    = {Nonlocal sliding parameters $\mu$ and $\gamma$ in flagellar mechanics}
}

% ---------------------------------------------------------------------------
% Additional Biomechanics and Rod Theory
% ---------------------------------------------------------------------------

@book{cowin2007tissue,
  title     = {Tissue Mechanics},
  author    = {Cowin, Stephen C. and Doty, Stephen B.},
  year      = {2007},
  publisher = {Springer},
  address   = {New York},
  note      = {Constitutive theory for biological tissues}
}

@article{oreilly2017modeling,
  title   = {Modeling nonlinear problems in the mechanics of strings and rods},
  author  = {O'Reilly, Oliver M.},
  journal = {Springer International Publishing},
  year    = {2017},
  doi     = {10.1007/978-3-319-50598-5},
  note    = {Advanced Cosserat rod theory and numerical methods}
}
```
