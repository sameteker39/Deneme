# Filtering or Modeling Disengagement? An Empirical Comparison of Effort-Moderated IRT and Item Response Tree Estimates in the First Digital TIMSS

## Abstract

The transition of the Trends in International Mathematics and Science Study (TIMSS) to a fully digital administration in 2023 has, for the first time, placed item-level response times at the center of routine large-scale calibration, and with them the long-standing problem of rapid guessing. Two families of models dominate current attempts to protect item and ability estimates from non-effortful responding: the effort-moderated item response theory (EM-IRT) model, which removes the influence of flagged responses, and item response tree (IRTree) models, which retain those responses and represent disengagement as a latent process alongside the trait of interest. The two paradigms rest on different assumptions and offer different by-products, yet the evidence bearing on their comparison is almost entirely simulation-based. Drawing on the TIMSS 2023 eighth-grade mathematics data from Türkiye (N = 4,925; 43 multiple-choice items), we estimated item difficulty, item discrimination, and ability under a standard two-parameter logistic model, an EM-IRT model, and a two-node IRTree, using the normative-threshold (NT10) rule to flag rapid guesses and four alternative rules—including a mixture-based rule—for a sensitivity analysis. Rapid guessing was infrequent but non-trivial (6.2% of timed responses) and behaved as theory predicts. Once rapid guesses were accounted for, EM-IRT and the IRTree accuracy node produced item and ability estimates that were, for practical purposes, indistinguishable (correlations ≥ .99), and both departed only modestly from the naïve model. The two paradigms parted company not in the parameters they share but in what the IRTree adds: a disengagement dimension that correlated negatively with mathematics ability. These conclusions were stable across threshold rules; only the disengagement–ability correlation was sensitive to how rapid guessing was operationalized. The findings clarify when the added machinery of an IRTree is worth its cost and when a filtering approach suffices.

*Keywords:* rapid guessing; test-taking effort; effort-moderated IRT; item response tree models; response time; TIMSS 2023

---

## 1. Introduction

Low-stakes achievement tests ask examinees to do something the testing situation gives them little reason to do: work hard on every item. When students see no personal consequence attached to their performance, some stop trying, and a portion of their responses ceases to reflect what they know. The most conspicuous form this withdrawal takes is rapid guessing—selecting an option so quickly that the examinee cannot plausibly have read the stem. Because rapid guesses are essentially random with respect to the measured trait, they introduce construct-irrelevant variance that can distort the very quantities test developers care about most: item parameters, ability estimates, and the group comparisons built on them (Wise, 2017; Wise & DeMars, 2006).

International assessments have absorbed this problem into their operational reality. Across the 2018 PISA science assessment in 71 countries, Rios et al. (2022) found that more than half of examinees rapid-guessed on roughly one in ten items, with the mean number of rapid guesses per examinee differing by as much as 500% between countries. Such differential disengagement is not a nuisance to be set aside; where it varies systematically across the groups a study exists to compare, it threatens the comparability that gives international assessment its purpose. The concern is sharpened because disengagement is rarely random: lower-performing examinees tend to rapid-guess more (Michaelides et al., 2020), so the responses any correction removes are concentrated among precisely those students whose scores are already most fragile.

Until recently, this line of inquiry was confined to computer-based programs that happened to record response latencies. TIMSS changed that calculus: with its 2023 cycle the study completed a transition to fully digital administration (von Davier et al., 2024), and item-level timing became available for a program that anchors the mathematics and science achievement of dozens of educational systems. For the first time, the methods developed to accommodate rapid guessing can be brought to bear on TIMSS at scale—and the choice among them acquires practical weight, because whichever model a national center or the international consortium adopts will shape the numbers that enter public reporting.

Two broad strategies for keeping rapid guessing from contaminating parameter estimates embody genuinely different philosophies. The first, the effort-moderated IRT (EM-IRT) model of Wise and DeMars (2006), treats a flagged response as uninformative and neutralizes it: the response function for rapid guesses is a flat line, so these responses contribute nothing to the likelihood and, in effect, drop out of estimation. The second declines to discard anything. Item response tree (IRTree) models decompose each observed response into a sequence of latent decisions—first whether the examinee engaged, then, conditional on engagement, whether the answer was correct—and estimate the trait governing each decision (Leventhal et al., 2023). Disengagement is not filtered out but modeled, as a latent propensity that lives alongside ability and may correlate with it.

These are not two implementations of one idea. EM-IRT buys parsimony at the price of a strong assumption (that flagged responses are ignorable once flagged) and returns nothing about disengagement itself. The IRTree relaxes that assumption and yields a person-level disengagement measure, but at the cost of a larger, multidimensional model whose behavior on operational data is less well charted. A practitioner choosing between them would want to know two things: whether the paradigms actually produce different item and ability parameters on real data or converge, and how much the answer depends on the analyst's unavoidably arbitrary choice of a response-time threshold.

The literature that could answer these questions is dominated by simulation. Careful studies have mapped how EM-IRT behaves under misclassification (Rios, 2021), under violations of its ignorability assumptions (Rios et al., 2021b), and against multidimensional variants that model the disengagement–ability link (Wang et al., 2024); others have traced how ignoring rapid guessing biases item parameters in specifiable directions (Alahmadi et al., 2025). This body of work is indispensable, but simulations answer only the question the analyst poses, under the model the analyst writes down. What remains comparatively thin is empirical evidence: side-by-side estimates from competing paradigms on the same operational responses, showing where they agree, where they diverge, and how fragile the difference is. The present study supplies such evidence for the setting where the question has newly become live—TIMSS 2023 mathematics in Türkiye, a context in which rapid guessing has been documented to intensify as items grow more difficult (Yılmaz Koğar et al., 2023).

We pursue three questions: how prevalent and valid rapid guessing is in these data; how item and ability parameters differ across a naïve two-parameter logistic (2PL) model, an EM-IRT model, and a two-node IRTree, and where the two disengagement-aware paradigms converge or diverge; and how sensitive those comparisons are to the response-time threshold. In answering them we aim less to crown a winner than to specify the conditions under which the choice between filtering and modeling matters, and those under which it does not.

## 2. Background

### 2.1 Rapid guessing and the integrity of measurement

Rapid guessing occupies a distinctive place among threats to validity because it is behavioral rather than dispositional and local rather than global. An examinee is not simply a "guesser"; the same person may reason carefully through one item and click through the next (Wise, 2017). This intra-individual variability is what makes response time such a powerful diagnostic: a response produced faster than any genuine solution process would allow betrays disengagement on that item, whatever the examinee did elsewhere. And because rapid guesses on multiple-choice items succeed only at chance, they compress the ability–correctness association that discrimination captures and distort difficulty in directions that depend on which examinees guessed and how often. Alahmadi et al. (2025) reconciled a body of conflicting results by showing that the direction and magnitude of this bias hinge jointly on item-level factors (true parameter values and the number of rapid guesses an item attracts) and the ability gap between rapid guessers and the rest. The distortion reaches beyond a single calibration: because item parameters propagate into adaptive item selection, equating, and standard setting, unaddressed rapid guessing can cascade through a whole testing operation (Alahmadi et al., 2025) and introduces linking error when rapid-guessing rates differ across groups placed on a common scale (Rios & Soland, 2021).

### 2.2 Identifying rapid guessing through response time

Turning response latencies into a usable indicator requires a rule that separates rapid guessing from genuine solution behavior. The dominant approach sets, for each item, a response-time threshold below which responses are flagged. Among the many procedures proposed, the normative threshold method of Wise and Ma (2012) has become a de facto standard: the threshold is a fixed fraction (commonly 10%, yielding the NT10 rule) of the item's mean response time, capped to prevent implausibly long thresholds on slow items. It is data-driven yet transparent, and it scales effortlessly to the hundreds of items in a large assessment.

Its very convenience raises the question of arbitrariness. Thresholds can also be set by inspecting the bimodality of response-time distributions, by cumulative-proportion arguments, by information-based criteria (Wise, 2019), or by fitting an explicit mixture of "rapid" and "solution" components to log-transformed times and cutting at a model-derived point—the mixture response-time quantile (MRTQ) logic. These procedures disagree: they flag markedly different proportions of responses, and all misclassify to some degree because the latency distributions of engaged and disengaged behavior overlap (Holopainen et al., 2026). Reassuringly, however, a meta-analysis of studies applying two or more procedures to the same test found that although threshold typology changed how many responses were flagged, it was largely unrelated to aggregate measurement properties and test performance; the act of filtering mattered more than the rule chosen (Rios et al., 2021a). Whether that reassurance transfers from filtering-based scoring to the model-based paradigms compared here is an open question we address directly.

### 2.3 The effort-moderated IRT model

The EM-IRT model of Wise and DeMars (2006) formalizes filtering within the item response framework. For solution behavior, correctness follows a standard IRT model; for rapid guesses, the probability of a correct answer is fixed at chance and is therefore flat in ability. Because a flat function has zero derivative, rapid guesses do not shift the peak of the likelihood, and the model estimates parameters as though those responses had been set aside—so that EM-IRT is closely allied, conceptually and computationally, to treating flagged responses as missing and calibrating on what remains—an approach that, among the available recodings, has proven a reliable way to recover item and ability parameters (Rios et al., 2023).

The model's appeal is its parsimony, and a decade of evidence has clarified its boundary conditions. When its assumptions hold, EM-IRT recovers item and ability parameters more accurately than models that ignore rapid guessing, and it is comparatively robust to the assumption that rapid guessing occurs randomly across items (Rios et al., 2021b). Its vulnerability lies elsewhere. When rapid-guessing propensity is related to ability, so that flagged responses are not missing at random with respect to the trait, bias re-enters, and although EM-IRT still outperforms the 2PL, its item estimates can shift by more than a fifth of a standard deviation under extreme conditions (Rios et al., 2021b); misclassification compounds the problem, with underclassification generally more damaging than overclassification (Rios, 2021). These findings motivated extensions that relax the ignorability assumption, including a two-step estimator and multidimensional variants that model the association between rapid-guessing propensity and ability (Wang et al., 2024). The same literature shows, however, that when rapid guessing is modest the elaborate and the simple approaches converge: effort-moderated scoring and a multidimensional competitor produced negligibly different ability estimates when rapid guesses made up 12% or less of responses (Rios et al., 2024). This convergence-under-modesty is a thread we pick up empirically.

### 2.4 IRTree models and the modeling of disengagement

Item response tree models approach the same data from the opposite direction. Rather than removing the influence of aberrant responses, they enlarge the model to represent the processes that generate them. An IRTree decomposes each response into a sequence of latent binary decisions arranged as a tree; each internal node is governed by its own IRT sub-model, and the traits underlying the nodes are estimated jointly (De Boeck & Partchev, 2012; Jeon & De Boeck, 2016). First used to disentangle content from response style in questionnaire data, separating whether a respondent agrees from how strongly (Meiser et al., 2019; Plieninger, 2020), the framework has since been extended to multidimensional and ordinal decision nodes (Meiser et al., 2019), to dynamic and longitudinal formulations that let response processes drift across a test (Ames et al., 2021; Merhof et al., 2023), and to diagnostics for parameter heterogeneity across respondents (Debelak et al., 2024), and it is now presented as a general tool for modeling response processes that traditional IRT cannot accommodate (Kim et al., 2025).

Applied to rapid guessing, the tree has a natural two-node form. The first node asks whether the examinee engaged with the item at all; the second, entered only when the first answer is "engaged," asks whether the response was correct. Leventhal et al. (2023) developed exactly this structure and showed what it buys. Because the model neither removes examinees nor assumes independence between disengagement and the trait, it yields less biased ability estimates for students who rapid-guess and, uniquely among the common corrections, returns a latent disengagement trait that can be studied in its own right. That second node is where the IRTree and EM-IRT meet: both condition the accuracy model on engaged responses, so the tree's accuracy node and the retained-response calibration of EM-IRT draw on the same information. What the tree adds is the first node and the correlation it permits between the two traits. Whether this structure repays its cost in the currency practitioners care about (item and ability parameters) rather than only in interpretive richness is, again, a question operational data can answer and simulation cannot fully settle.

### 2.5 The present study

The strands above converge on a gap. Simulation evidence robustly shows that disengagement-aware models outperform naïve ones, that EM-IRT and multidimensional alternatives converge when rapid guessing is modest, and that threshold choice matters little for aggregate properties under filtering. What is scarce is empirical evidence that places EM-IRT and an IRTree side by side on operational large-scale data, reports the item and ability parameters each yields, and asks where they agree and diverge—and, to our knowledge, none exists for the first digital TIMSS, in a system where rapid guessing is known to respond to item difficulty (Yılmaz Koğar et al., 2023). We address this gap with three research questions:

- **RQ1.** How prevalent is rapid guessing in the TIMSS 2023 mathematics multiple-choice items for Türkiye, and does the flagged behavior display the accuracy signature the construct requires?
- **RQ2.** How do item discrimination, item difficulty, and ability estimates compare across a naïve 2PL model, an EM-IRT model, and a two-node IRTree, and where do the two disengagement-aware paradigms converge or diverge?
- **RQ3.** How robust are these comparisons to the response-time threshold rule used to flag rapid guessing, from conservative normative rules to a more liberal mixture-based rule?

## 3. Method

### 3.1 Data and sample

Data came from the eighth-grade mathematics assessment of TIMSS 2023, the first cycle administered entirely on digital devices (von Davier et al., 2024). We analyzed the Türkiye sample of 4,925 students. TIMSS employs a rotated booklet design in which each student receives only a subset of the item pool, so the person-by-item matrix is sparse by design: each student answered a handful of the items analyzed here, while each item was answered by roughly 700 students. This structure suits item calibration, which pools information across all respondents to an item, and is handled natively by marginal-maximum-likelihood estimation, which conditions each person's likelihood on the items that person saw; it does, however, limit the precision of individual ability estimates, a point we revisit below.

### 3.2 Instruments

We restricted the analysis to multiple-choice mathematics items, for which the chance probability of a correct guess is well defined (one in four for four-option items) and the accuracy signature of rapid guessing is therefore interpretable. Forty-three items met this criterion and carried the timing needed to flag rapid guessing. Each response recorded the selected option, scored dichotomously against the answer key, and the first-response time in seconds—the interval from item presentation to the initial answer, the conventional latency for rapid-guessing research because it captures first engagement rather than later revisiting. Timing records reserved for not-reached or omitted responses were coded as missing.

### 3.3 Flagging rapid guessing

Rapid guesses were identified with the normative threshold method (Wise & Ma, 2012). For each item the threshold was set at 10% of the mean first-response time, capped at 10 seconds (the NT10 rule); responses at or below the item threshold were flagged as rapid guesses and the remainder as solution behavior. Before computing each item mean, we winsorized the upper tail of the response-time distribution at two standard deviations, so that a small number of very long latencies could not inflate the mean and, through it, the threshold. Items with fewer than 30 timed responses were not assigned a threshold.

Because any single rule is to some degree arbitrary, we repeated the entire analysis under four alternatives: NT10 without winsorizing; a stricter normative rule using 5% of the mean (NT5); a fixed five-second threshold; and a mixture-based rule (MRTQ) that fits a two-component normal mixture to each item's log response times and cuts at the quantile corresponding to the weight of the lower-mean, rapid-guessing component (Holopainen et al., 2026). These rules span conservative, data-free cutoffs to liberal, distribution-based ones, bracketing the decisions a reasonable analyst might make.

### 3.4 Models and estimation

Three models were fitted to the same responses. The **standard 2PL** calibrated all scored responses, ignoring rapid guessing, and served as the baseline against which the cost of doing nothing could be read. The **EM-IRT** model implemented effort moderation: flagged responses were rendered uninformative and the 2PL was estimated from the retained, solution-behavior responses. Fixing the rapid-guess response function at chance and treating flagged responses as missing are equivalent for item and ability estimation, because the flat function contributes a constant to the likelihood that does not depend on the parameters; we note this equivalence explicitly, as it is the reason the EM-IRT and IRTree accuracy estimates can be expected to lie close together. The **IRTree** represented each timed response as two pseudo-items: a first node coding whether the response was a rapid guess and a second coding correctness, defined only for engaged responses. The nodes loaded on two correlated latent dimensions, disengagement propensity and mathematics ability, estimated jointly with a freely estimated correlation.

All models were estimated by marginal maximum likelihood with the expectation–maximization algorithm; person parameters were obtained as expected a posteriori (EAP) estimates, and the common two-parameter logistic metric placed all discrimination, difficulty, and ability estimates on comparable scales. Estimation used the mirt package in R (Chalmers, 2012). Because the IRTree nodes are between-item multidimensional, with each pseudo-item loading on one dimension, the accuracy node's item parameters are governed by the engaged responses, exactly as in EM-IRT, while the correlation between dimensions is informed by the joint pattern across nodes.

### 3.5 Comparison and analysis plan

We compared the models on the parameters that enter operational use. For discrimination and difficulty (across all 43 items) and for ability (across the sample), we computed Pearson and Spearman correlations across models together with mean differences and the root mean squared difference (RMSD), which captures the typical magnitude of disagreement in the parameter's own metric. The disengagement dimension, available only from the IRTree, was summarized by its correlation with ability. The sensitivity analysis re-estimated the EM-IRT and IRTree models under each of the five threshold rules, reporting for every rule the rapid-guessing rate, the accuracy of flagged and unflagged responses, the mean discrimination of the accuracy node, the disengagement–ability correlation, and the correlation of each rule's item and ability estimates with those from the primary analysis.

## 4. Results

### 4.1 Prevalence and validity of rapid guessing

Under the primary NT10 rule, 6.2% of timed responses were flagged as rapid guesses—modest yet far from negligible, and within the range reported for multiple-choice items in comparable international settings, where per-item rates below 6% to 10% are typical (Michaelides et al., 2020; Rios et al., 2022). More important than the rate is its behavioral signature. Flagged responses were correct 42% of the time, against 54% for solution behavior. That flagged responses were substantially less accurate is the pattern the construct demands: rapid guesses should approach chance, unflagged responses should exceed it, and the ordering should be unambiguous (Wise, 2017). Flagged accuracy did not descend all the way to the 25% chance level, a residue consistent with the fact that thresholds misclassify at the margins where engaged and disengaged latency distributions overlap (Holopainen et al., 2026); some hurried but genuine solution behavior is inevitably swept in. The direction and size of the gap nonetheless support the validity of the classification.

### 4.2 Item parameters across models

The central question is what accounting for rapid guessing does to the item parameters, and how much the two disengagement-aware paradigms differ in doing it. Table 1 reports the comparisons, and two results stand out. First, accounting for rapid guessing moved the item parameters in the expected direction but only modestly. Mean discrimination rose from 1.43 under the naïve 2PL to 1.45 under EM-IRT and 1.46 under the IRTree accuracy node—the attenuation caused by chance-level guesses is relieved once those guesses are removed, but the shift is small at this prevalence—while mean difficulty was similarly stable (−0.11, −0.13, and −0.11). Correlations between the naïve and effort-moderated item parameters were high (r = .99 for both parameters) and the typical disagreement small (RMSD = 0.10 for discrimination, 0.06 for difficulty). Yet aggregate stability masked movement for particular items: the hardest items, which attract proportionally more rapid guesses and on which lucky guesses most inflate the apparent proportion correct, tended to become harder once those guesses were discounted, and the item whose difficulty changed most rose from 1.80 to 2.07 logits.

The second, more striking result concerns the two disengagement-aware models, which produced essentially the same item parameters: the correlation was 1.00 for difficulty and .998 for discrimination, with RMSDs of 0.03 and 0.04. This near-identity follows from the models' structure. Both condition the accuracy model on the same engaged responses, and at a prevalence of 6% the information the IRTree draws from modeling disengagement jointly does little to move the accuracy node's item parameters. The empirical convergence thus confirms, on operational data, the pattern simulation had led us to expect when rapid guessing is modest (Rios et al., 2024): for item calibration, filtering and modeling arrive at the same place.

**Table 1.** *Cross-model comparison of item and ability parameters (primary NT10 analysis).*

| Comparison | Parameter | r | Mean difference | RMSD |
|---|---|---|---|---|
| Standard → EM-IRT | Difficulty (b) | .998 | −0.017 | 0.062 |
| Standard → IRTree (ACC) | Difficulty (b) | .998 | 0.001 | 0.070 |
| EM-IRT → IRTree (ACC) | Difficulty (b) | 1.000 | 0.018 | 0.025 |
| Standard → EM-IRT | Discrimination (a) | .991 | 0.019 | 0.096 |
| Standard → IRTree (ACC) | Discrimination (a) | .992 | 0.027 | 0.090 |
| EM-IRT → IRTree (ACC) | Discrimination (a) | .998 | 0.008 | 0.037 |
| Standard → EM-IRT | Ability (θ) | .987 | −0.002 | 0.125 |
| Standard → IRTree (ACC) | Ability (θ) | .980 | 0.002 | 0.156 |
| EM-IRT → IRTree (ACC) | Ability (θ) | .994 | 0.007 | 0.082 |

*Note.* ACC = accuracy node of the IRTree. RMSD = root mean squared difference. Item parameters compared across 43 items; ability across the students estimable under each model pair (n = 4,843–4,885).

### 4.3 Ability estimates

The picture for ability mirrors that for the items. Estimates from the three models were strongly correlated (.987 between the naïve and EM-IRT models and .994 between EM-IRT and the IRTree), yet the agreement was not perfect, and the disagreement was concentrated where theory says it should be. Because rapid guessing is unevenly distributed across students, correcting for it leaves most estimates untouched and adjusts those of the students who rapid-guessed. The mean difference between models was near zero, as the common metric requires, but the RMSD of 0.08 to 0.16 logits reflects real movement for the disengaged minority—small for a program reporting population distributions and group means, but not negligible for inferences about the individual students for whom the disengagement-aware models were designed (Leventhal et al., 2023).

### 4.4 The disengagement dimension

If the two paradigms agree so closely on item and ability parameters, what does the IRTree offer that EM-IRT does not? The answer lies in the parameter they do not share. The IRTree estimated a latent disengagement propensity whose correlation with mathematics ability was −0.33: students more prone to rapid guessing tended to be lower in ability—the relationship EM-IRT assumes away and that, when strong, is precisely the condition under which effort-moderated filtering begins to bias its own estimates (Rios et al., 2021b). Here the relationship is moderate rather than strong, consistent with the paradigms' close agreement on the shared parameters: the ignorability assumption is violated, but not enough at this prevalence to drive the models apart. The disengagement dimension is therefore best understood not as a competitor to the ability estimate but as an additional, substantively meaningful trait, a description of who withdraws and how that withdrawal relates to proficiency, that filtering discards by construction.

### 4.5 Sensitivity to the threshold rule

The robustness analysis asked whether these conclusions depend on how rapid guessing was operationalized. They do not, with one instructive exception. Across the five rules the rapid-guessing rate varied fourfold, from 3.3% under the strict NT5 rule to 14.1% under the liberal MRTQ rule (Table 2)—the expected consequence of moving from conservative to liberal cutoffs, echoing the meta-analytic finding that procedures flag materially different proportions of responses (Rios et al., 2021a). Yet the item and ability parameters were almost entirely insensitive to this variation: across every rule, the EM-IRT discrimination, difficulty, and ability estimates correlated with the primary analysis at .99 or above. The one parameter that responded was the disengagement–ability correlation, which strengthened from −0.16 under NT5 to −0.43 under MRTQ. This is not a defect but a signature: how strongly disengagement relates to ability depends on how much disengaged behavior a rule captures, so a liberal rule naturally recovers a stronger relationship.

The practical message is twofold. For the parameters that enter reporting, the threshold rule is close to immaterial—an empirical extension, to the model-based paradigms, of the conclusion that accounting for rapid guessing matters more than the rule chosen (Rios et al., 2021a). For the disengagement dimension that gives the IRTree its distinctive value, the rule is consequential, and analysts who interpret that dimension should treat their threshold choice as a substantive decision rather than a technical default.

**Table 2.** *Sensitivity of key quantities to the response-time threshold rule.*

| Threshold rule | RG rate | Accuracy \| RG | Accuracy \| solution | Mean *a* (EM-IRT) | Disengagement–ability *r* | *a* corr. with NT10 | θ corr. with NT10 |
|---|---|---|---|---|---|---|---|
| NT10 (winsorized; primary) | .062 | .42 | .54 | 1.45 | −0.33 | 1.00 | 1.00 |
| NT10 (unwinsorized) | .064 | .42 | .54 | 1.46 | −0.34 | .999 | .999 |
| NT5 (winsorized) | .033 | .48 | .54 | 1.44 | −0.16 | .996 | .994 |
| Fixed 5 s | .071 | .42 | .54 | 1.46 | −0.34 | .998 | .996 |
| MRTQ (mixture) | .141 | .38 | .56 | 1.48 | −0.43 | .992 | .978 |

*Note.* RG = rapid guessing. Correlations in the final two columns are between each rule's EM-IRT estimates and those of the primary NT10 analysis.

## 5. Discussion

### 5.1 Principal findings

This study asked whether two philosophically distinct responses to rapid guessing—filtering it out with an effort-moderated model, or modeling it as a latent process with an item response tree—produce different item and ability parameters on operational data, and how much any difference depends on the threshold choice. The answers are clarifying. Rapid guessing in the TIMSS 2023 mathematics items for Türkiye was modest but real and behaved as the construct requires. Once it was accounted for, the two paradigms produced item and ability estimates that were, for every practical purpose, the same, and both departed only slightly from a model that ignored rapid guessing altogether. Where they diverged was not in the parameters they share but in what the tree adds: a disengagement dimension, negatively related to ability, that filtering cannot yield. These conclusions held across a wide range of threshold rules; only the disengagement–ability correlation moved with the rule, in the direction its logic predicts.

### 5.2 Interpretation and connection to prior work

The convergence of EM-IRT and the IRTree accuracy node is neither surprising nor trivial once its source is understood. Both models estimate the accuracy structure from engaged responses, and at the prevalence observed here the extra information the tree draws from its disengagement node does little to reshape the accuracy parameters. The result gives empirical flesh to a pattern simulation had suggested in the abstract: effort-moderated scoring and its multidimensional relatives coincide when rapid guessing is modest (Rios et al., 2024; Wang et al., 2024). Our contribution is to show that the coincidence survives contact with operational data, with all the messiness (matrix sampling, sparse ability information, imperfect classification) that such data bring.

The near-invariance of item and ability parameters to the threshold rule extends a known result into new territory. Rios et al. (2021a) established, for filtering-based scoring, that threshold typology has little bearing on aggregate measurement properties even as it changes how many responses are flagged; we find the same insensitivity for the model-based paradigms, and for both item and ability estimates. The disengagement dimension is the exception that proves the rule: it is precisely the quantity that depends on how much disengaged behavior a threshold captures, and so precisely the quantity that moves when the threshold changes. Its sensitivity is a reason to interpret it cautiously, not to distrust the stable parameters around it.

Together, these results locate the choice between the paradigms where it belongs. The two are not rivals estimating the same quantity with different accuracy; on the shared quantities they are, at modest prevalence, interchangeable. They differ in scope. EM-IRT answers a narrow, well-posed question (what are the item and ability parameters, cleansed of rapid guessing?) with minimal machinery and a strong assumption. The IRTree answers it equally well while also answering one EM-IRT cannot pose: who disengages, and how does disengagement relate to what we are measuring? That question is not idle, for disengagement is itself of interest to systems that wish to reduce the withdrawal of students from low-stakes assessment (Leventhal et al., 2023; Nagy et al., 2022), and its relationship to ability bears on the fairness of comparisons among groups that disengage at different rates (Rios et al., 2022).

### 5.3 Implications for practice

For a program whose sole aim is to protect item calibration and ability estimation from rapid guessing, the evidence counsels reassurance and parsimony. At the prevalence typical of well-administered international mathematics assessments, an effort-moderated model suffices; the more elaborate tree yields the same parameters, and the threshold rule can be chosen for transparency rather than agonized over, because the reported parameters scarcely notice the difference—a practically consequential simplification for centers that must calibrate many items under time pressure.

For a program that wants more than clean parameters, whether to monitor disengagement, study its correlates, or reason about its implications for comparability, the calculus changes. The IRTree delivers the same protected parameters and, at little additional cost, a disengagement trait that can be modeled, predicted, and compared across groups. Such a program should treat the threshold rule as a substantive choice, however, because the disengagement dimension is the one quantity that rule genuinely affects; the prudent course is to report the disengagement–ability relationship across a range of rules, as we have done, rather than rest an interpretation on a single default.

### 5.4 Limitations

Several limitations mark the boundaries of these conclusions. The analysis was confined to multiple-choice items, for which chance performance is defined and the accuracy signature of rapid guessing is interpretable; constructed-response items, where disengagement manifests as omission or perfunctory answers, require different treatment and may not converge in the same way. The rotated booklet design that makes TIMSS efficient renders individual ability estimates imprecise, because each student answers few items, so our conclusions about ability agreement are most secure at the level of distributions and group comparisons; the small individual-level differences between models deserve closer study in designs with denser person information. Rapid guessing here was modest, and the convergence we document should weaken as prevalence rises and the disengagement–ability relationship strengthens—the very conditions under which EM-IRT's ignorability assumption begins to bite (Rios et al., 2021b) and the tree's willingness to model that relationship should matter more. Finally, all thresholds misclassify: the residual above-chance accuracy of flagged responses shows that our "rapid guesses" contain some engaged behavior, and approaches that classify at the response-by-examinee level or fold response time into the model rather than dichotomizing it first (Deribo et al., 2021; Ulitzsch et al., 2019) may sharpen the picture. Whether the convergence generalizes to systems with heavier disengagement (Rios et al., 2022) is an empirical question a single system and grade cannot settle.

### 5.5 Future directions

The most direct extension is to trace the boundary at which the paradigms separate. Because their convergence is a function of prevalence and of the disengagement–ability relationship, a design spanning systems, grades, or subjects that differ in these quantities could map where filtering ceases to suffice and modeling begins to pay. The comparison would be strengthened by moving beyond dichotomize-then-model logic toward models that treat response time as continuous information about engagement (Deribo et al., 2021; Lu et al., 2020; Ulitzsch et al., 2019), and by folding item features known to drive disengagement—difficulty, cognitive demand, position (Yılmaz Koğar et al., 2023)—into the disengagement node itself, so that the tree explains not only who disengages but where. As TIMSS accumulates digital cycles, the disengagement dimension could also be tracked over time and related to the contextual variables the study already collects, turning a psychometric correction into a substantive account of test-taking motivation.

## 6. Conclusion

Faced with rapid guessing, an analyst can filter it away or model it in. On the first digital TIMSS mathematics data from Türkiye, these two courses led to the same item and ability parameters—parameters a threshold rule could scarcely disturb. The effort-moderated model and the item response tree are therefore not competitors on the ground they share but instruments of different scope. When the goal is clean calibration at modest prevalence, the simpler instrument suffices and the choice of threshold is a detail; when the goal also includes understanding disengagement, the tree earns its added structure by measuring something the filter throws away. Knowing which goal is in view, and that the shared parameters are safe either way, lets practitioners spend their modeling effort where it actually changes the answer.

---

## References

Alahmadi, S., et al. (2025). From item estimates to test operations: The cascading effect of rapid guessing. *Journal of Educational Measurement.* https://doi.org/10.1111/jedm.70010

Ames, A. J., et al. (2021). Modeling changes in response style with longitudinal IRTree models. *Multivariate Behavioral Research.*

Chalmers, R. P. (2012). mirt: A multidimensional item response theory package for the R environment. *Journal of Statistical Software, 48*(6), 1–29.

Debelak, R., et al. (2024). Investigating heterogeneity in IRTree models for multiple response processes with score-based partitioning. *British Journal of Mathematical and Statistical Psychology.*

De Boeck, P., & Partchev, I. (2012). IRTrees: Tree-based item response models of the GLMM family. *Journal of Statistical Software, 48*(1), 1–28.

Deribo, T., Kroehne, U., & Goldhammer, F. (2021). Model-based treatment of rapid guessing. *Journal of Educational Measurement, 58*(2), 281–303.

Holopainen, S., et al. (2026). Misclassification produced by rapid-guessing identification methods and their suitability under various conditions. *Educational and Psychological Measurement.*

Jeon, M., & De Boeck, P. (2016). A generalized item response tree model for psychological assessments. *Behavior Research Methods, 48*(3), 1070–1085.

Kim, N., et al. (2025). Digital Module 37: Introduction to item response tree (IRTree) models. *Educational Measurement: Issues and Practice.*

Leventhal, B. C., et al. (2023). An illustration of an IRTree model for disengagement. *Educational and Psychological Measurement.*

Lu, J., et al. (2020). A mixture model for responses and response times with a higher-order ability structure to detect rapid guessing behaviour. *British Journal of Mathematical and Statistical Psychology.*

Meiser, T., et al. (2019). IRTree models with ordinal and multidimensional decision nodes for response styles and trait-based rating responses. *British Journal of Mathematical and Statistical Psychology.*

Merhof, V., et al. (2023). Dynamic response strategies: Accounting for response process heterogeneity in IRTree decision nodes. *Psychometrika.*

Michaelides, M. P., et al. (2020). The relationship between response-time effort and accuracy in PISA science multiple-choice items. *International Journal of Testing.*

Nagy, G., et al. (2022). The role of rapid guessing and test-taking persistence in modelling test-taking engagement. *Journal of Computer Assisted Learning.*

Plieninger, H. (2020). Developing and applying IR-tree models: Guidelines, caveats, and an extension to multiple groups. *Organizational Research Methods.*

Rios, J. A. (2021). Assessing the accuracy of parameter estimates in the presence of rapid guessing misclassifications. *Educational and Psychological Measurement.*

Rios, J. A., et al. (2021a). Does the choice of response time threshold procedure substantially affect inferences concerning the identification and exclusion of rapid guessing responses? A meta-analysis. *Large-scale Assessments in Education, 9*(1).

Rios, J. A., et al. (2021b). Parameter estimation accuracy of the effort-moderated item response theory model under multiple assumption violations. *Educational and Psychological Measurement.*

Rios, J. A., et al. (2022). An investigation of item, examinee, and country correlates of rapid guessing in PISA. *International Journal of Testing.*

Rios, J. A., & Soland, J. (2021). Investigating the impact of noneffortful responses on individual-level scores: Can the effort-moderated IRT model serve as a solution? *Educational and Psychological Measurement.*

Rios, J. A., et al. (2023). A comparison of response time threshold scoring procedures in mitigating bias from rapid guessing behavior. *Educational and Psychological Measurement.*

Rios, J. A., et al. (2024). Is effort-moderated scoring robust to multidimensional rapid guessing? *Educational and Psychological Measurement.*

Ulitzsch, E., et al. (2019). A hierarchical latent response model for inferences about examinee engagement in terms of guessing and item-level non-response. *British Journal of Mathematical and Statistical Psychology.*

von Davier, M., Kennedy, A., Reynolds, K., Fishbein, B., Khorramdel, L., Aldrich, C., Bookbinder, A., Bezirhan, U., & Yin, L. (2024). *TIMSS 2023 international results in mathematics and science.* TIMSS & PIRLS International Study Center, Boston College.

Wang, B., Huggins-Manley, A. C., Kuang, H., & Xiong, J. (2024). Enhancing effort-moderated item response theory models by evaluating a two-step estimation method and multidimensional variations on the model. *Educational and Psychological Measurement.*

Wise, S. L. (2017). Rapid-guessing behavior: Its identification, interpretation, and implications. *Educational Measurement: Issues and Practice, 36*(4), 52–61.

Wise, S. L. (2019). An information-based approach to identifying rapid-guessing thresholds. *Applied Measurement in Education.*

Wise, S. L., & DeMars, C. E. (2006). An application of item response time: The effort-moderated IRT model. *Journal of Educational Measurement, 43*(1), 19–38.

Wise, S. L., & Ma, L. (2012, April). *Setting response time thresholds for a CAT item pool: The normative threshold method.* Paper presented at the annual meeting of the National Council on Measurement in Education, Vancouver, Canada.

Yılmaz Koğar, E., et al. (2023). Examination of response time effort in TIMSS 2019: Comparison of Singapore and Türkiye. *International Journal of Assessment Tools in Education.*
