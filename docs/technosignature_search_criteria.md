# Technosignature search criteria from an archaeological framework

Oct 13, 2026 · @Dylan Foley

## Purpose

This brief gives a small technosignature survey a set of search criteria drawn from archaeology. Some will look counterintuitive to astronomers, because they start from what technology is physically, not from what our own technology looks like.

The criteria are under Search criteria. The sections before them give the reasoning, and Limits and open questions says what the framework cannot do. The criteria are proposals that follow from the framework, not results that have been tested on survey data.

## The shared problem

Neither SETI nor archaeology has a computable definition of technology, yet both exist to detect it. Each substitutes a list. SETI uses a proxy list: narrowband emission, waste heat, atmospheric disequilibrium, megastructure transits. Archaeology uses typology, a list of things already agreed to be tools.

Every item on SETI's list is something our own technology does, projected outward. That is anthropocentric enumeration standing in for a criterion. Both fields recognise instances of a category that neither has defined.

## The proposed definition

Technology is structured, controlled motion at a boundary that transforms the agent's effective environment. Its physical quantity is Θ(μ): the change in force, distance and speed that a motion pattern μ produces.

```latex
\Theta(\mu) = [\Delta F,\ \Delta d,\ \Delta v]
```

Θ is a dimensional physical quantity, not a measure of improbability. Because it names a mechanism, it predicts which signatures should appear, instead of flagging anomalies after the fact.

## What archaeology contributes

Archaeology is the only available labelled calibration set for technosignature discrimination. SETI has no confirmed technological generators beyond our own, so it cannot validate an origin criterion on the sky. The terrestrial record holds a very large number of confirmed technological traces with independent ground truth.

The two fields have complementary bottlenecks. Archaeology is data-rich and integration-poor. SETI is integration-adequate, because it inherits astronomy's coordinate and time frames, but positive-example-poor.

| Field | Time coordinate | Evidence | Generators |
| --- | --- | --- | --- |
| Astronomy | Inherited by propagation | Remote observation | Natural |
| Archaeology | Must be constructed | Contact traces | Confirmed technological |
| SETI | Inherited by propagation | Remote, no contact possible | Unconfirmed |

The NASA volume edited by Vakoch, Archaeology, Anthropology, and Interstellar Communication (2014), is a standing call from the SETI community for archaeological input to technosignature interpretation.

## The key restriction

Only deposited-regime technology is observable at interstellar distance. The framework separates two regimes. In the entrained regime the transformation happens during the motion, like a door while it is being pushed. In the deposited regime a standing configuration keeps transforming its environment with no ongoing motion from the maker.

Nobody resolves a hand axe in use at four parsecs, so entrained technology is out of reach. What remains is megastructures, altered atmospheres and waste heat from persistent infrastructure. SETI is therefore archaeology restricted to the deposited regime at planetary scale.

This derivation explains the shape of the technosignature literature. The proxy list is deposited-regime because the observing geometry permits nothing else.

## Position against assembly theory

Assembly theory is a detector without a mechanism, and Θ is a mechanism that generates detectors. Marshall, Cronin and Walker built assembly theory to tell biological from non-biological molecules in exoplanet atmospheres, and Jaeger's commentary extends it to technology: monitor complexity and abundance, then infer selection. It is computable from remote data such as mass spectra, which is a real practical advantage.

Its neutrality is also its ceiling. It shows that something biased the outcome and, by design, assumes nothing about what did. The wider family of measures (entropy, Kolmogorov complexity, Bennett's logical depth, assembly index) all quantify departure from randomness. Geological and chemical processes also produce long histories, so a substrate-neutral measure cannot separate a technological generator from a slow mineralogical one.

Θ is not in that family. It specifies a physical mechanism, so the two are orthogonal rather than competing: one could compute the assembly index of a Θ-trace.

## Search criteria

Seven criteria follow from the framework, and the first five differ from standard practice. Each is stated with what it changes.

1. **Rank signatures by persistence, not by resemblance to human technology.** Standard practice enumerates what our technology does. Persistence decides what a remote observer can see at all, so the useful output is a persistence-ranked taxonomy of remotely observable proxies.
2. **Look for ratcheting, not just symmetry or complexity.** Natural processes make symmetry. Only life makes symmetry that ratchets: irreversible accumulation, recombination and chain-deepening. This is a candidate origin criterion grounded in the physical definition, and it is not yet tested on survey data.
3. **Treat complexity-plus-abundance flags as necessary, not sufficient.** Require a proposed mechanism that predicts a second, independent observable, so that an anomaly becomes a testable hypothesis.
4. **Build the natural null from slow generators.** The null model must include geological and chemical processes with long histories, because those are what substrate-neutral measures cannot exclude.
5. **Do not assume early technology looks industrial.** The proxy list models only industrial civilisation. A proposed library of technosignature spectra would characterise earlier regimes (pre-lithic, Oldowan, Acheulean, composite tools, agricultural transformation, industrial) by stability duration, transition geometry and relationship to environmental perturbation.
6. **State the blind spot.** A null result says nothing about civilisations whose technology is entrained and leaves no standing configuration. Report that as a limit of the survey, not as evidence of absence.
7. **Validate on labelled terrestrial positives before using sky data.** Where a discriminator can be tested on archaeological or terrestrial cases with known origin, test it there first. Fit the simplest generator first and let the residual decide whether to escalate to a richer model; this method comes from the framework's trace-inversion pipeline and has not yet been adapted to survey data.

## Limits and open questions

Θ cannot be measured at interstellar distance. The framework supplies a theory of what technology is and a way to generate and rank remote proxies. It is not a drop-in remote detector, and assembly index keeps the practical advantage there.

- The exoplanet side has no validation domain. The hominin record is a single case of the process, and exoplanets are none, so claims about the sky are hypotheses within a project this size.
- The ratcheting criterion and the spectra library are proposals. Neither has been tested against survey-type data.
- Open: how to word "deposited regime" for actively running emitters such as narrowband transmitters, which a specialist may argue are not standing configurations.
- Open: how to rank persistence quantitatively for this survey's instrument and target set.

## Provenance and references to verify

Nothing here was looked up. The content comes from the author's earlier discussions of the framework, and the references below need checking before they are cited to the survey team.

- The author's conference paper and the manuscript defining technology physically: add full citations and current status.
- Vakoch (ed.), Archaeology, Anthropology, and Interstellar Communication, NASA, 2014.
- Assembly theory: Marshall, Cronin and Walker; Sharma et al. (2023); Jaeger's commentary (2024); Zenil's critique. Check for formal replies since these appeared.
- Bennett (1988) on logical depth.
