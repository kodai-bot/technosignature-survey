# Research programme: an archaeological technosignature survey with TESS

This turns the criteria in `technosignature_search_criteria.md` (the "brief") into work this repository can carry out.

**Instrument:** TESS photometry. It records each star's brightness every 2–10 minutes for about 27 days per sector. Many stars have been observed in several sectors since 2018.

**Deposited regime only:** the brief restricts remote detection to standing configurations. For this instrument that means **anything that persistently blocks, re-emits or modulates starlight**. Waste heat and altered atmospheres need other instruments, so here they appear only as cross-checks (WP4).

Status of the pieces built so far (Oct 2026):

| Piece | State |
|---|---|
| `technosig/scan.py` | short-dip event catalogue plus explain-away classifier; 1,095 Sector 59 stars swept |
| `technosig/longdips.py` | long-dip (Boyajian-type) pass; works, but the SAP route needs proper systematics correction |
| `followup/locate_dip_source.py` | locates which star a dip is on, using pixel-level data |
| `results/FIELD_NOTES.md` | first sweep's results; BD+67 274 eclipse traced through two wrong hosts |

---

## How each criterion becomes something we can run

### C1. Rank by persistence → *what survives long enough for TESS to see it*

For a photometric survey, persistence reduces to three measurable properties of a signal:

- **Lifetime** of the configuration: how long it exists at all.
- **Duty cycle**: the fraction of the time it is in front of the star, or otherwise visible.
- **Repeatability**: whether it comes back, so that a second sector can confirm it.

**Work:** write a short persistence table for candidate deposited proxies on this instrument:
- orbiting occulters (single bodies, swarms, rings)
- slowly growing occultation (an accumulating structure)
- persistent shape anomalies in transits
- secular dimming

Rate each one on those three properties and on detectability per sector. That table is the answer to the brief's open question "how to rank persistence quantitatively for this instrument".

**It also tells us how to schedule observations:** low-duty-cycle signals need many sectors, so multi-sector stars rank above single-sector ones.

### C2. Ratcheting → *irreversible, accumulating change across epochs*

**The photometric form of a ratchet:**
- the system moves from one state to another and does not return;
- each new state builds on the previous one;
- the steps get deeper over time.

On a light curve that could show up as occultation that grows sector by sector, or a transit shape that changes one way only.

**Warning: monotonic change alone is not a ratchet.** Plenty of natural processes accumulate and don't reverse:
- stars evolve in one direction;
- dust builds up;
- dust can clear and not come back.

So the operational test has to require **accumulation with recombination or deepening**: later states include and extend earlier ones. A simple one-way drift isn't enough.

**Work:**
- (a) Define a ratchet statistic on a time series. Candidates:
  - time-reversal asymmetry;
  - step detection with a no-return test;
  - the depth of nesting of successive states.
- (b) Validate it on terrestrial series of known origin **before** using any sky data (see C7).
- (c) Then apply it to stars with long baselines, combining TESS multi-sector data with archive photometry (ASAS-SN, ZTF, and DASCH photographic plates going back about 100 years). Boyajian's star's long-term fading is the obvious test case.

### C3. A second, independent observable → *every anomaly must predict something else*

For a dip, the standard independent observables are:

- **Chromaticity.** Dust makes dips deeper in blue light than in red; a solid opaque body dims all colours equally. Needs multi-band follow-up or archive colours.
- **Infrared excess.** A structure that absorbs starlight must re-radiate it as heat. Check by cross-matching with WISE and Gaia.
- **Repeat timing.** A mechanism with an orbit predicts when the dip comes back.
- **Pixel location.** Already done: `locate_dip_source.py`. BD+67 274 shows why this test comes first.

**Work:** add these as columns in the candidate table. A candidate is promoted only when a proposed mechanism names which of these it predicts.

### C4. Natural null from slow generators → *build the natural assemblage before looking for the artificial one*

The classifier's current mundane explanations are only the instrument and common stellar behaviour. The brief asks for slow natural generators too. For dimming events these are the astronomical equivalent of geology:

- dust from colliding planetesimals, and disintegrating rocky planets;
- exocomets;
- "dipper" young stars still surrounded by disks;
- long-period binaries with a disk (the ε Aurigae type);
- ring systems (the J1407 type);
- spots and starspot evolution.

**Work:**
- Build a reference catalogue of known examples of each class, observed by TESS where possible.
- Extend `classify()` so each class has a description of its expected shape, duration and colour.
- Report residue only after all of them have had their turn.

This is the reference collection a finds specialist works from.

### C5. Don't assume early technology looks industrial → *state what a non-industrial deposited signature would look like here*

On TESS photometry, most pre-industrial regimes are out of reach. Land use and agriculture change a planet's reflectivity, far below this instrument's sensitivity.

**Work:** record each regime from the brief's spectra library (stability duration, transition geometry, relationship to environmental disturbance), mark whether this instrument could observe it, and carry the unobservable ones into the blind-spot statement (C6). The library itself belongs to the theory side of the project, not this repo.

### C6. State the blind spot → *injection–recovery tests, plus a regime statement*

A null result is only meaningful with a known detection efficiency.

**Work:**
- **Injection–recovery.** Inject synthetic signals into the real light curves:
  - box, V and asymmetric dips;
  - multi-component dips;
  - slow fades;
  - a range of depths, durations and shapes.

  Run the full pipeline and map the recovery fraction.
- **What a null result then says.** Something like: "no occultation deeper than X% lasting Y–Z hours in N stars, at R% completeness".
- **Regime statement.** Separately, list the regimes the survey cannot see at all: entrained technology, and everything in C5 marked unobservable.

### C7. Validate on labelled positives; fit the simplest generator first

**Two validation tracks:**

- **Discriminators on terrestrial series.** Test the ratchet statistic (C2) on time series of known origin:
  - technological: the Paige & Perreault procedural-unit sequence;
  - natural slow processes: sediment, geochemical and climate-proxy series;
  - astronomical natural: long-term variable-star records.

  It has to separate these before it touches the sky.
- **The pipeline on astronomical labelled sets.** Already started:
  - TOIs and known planets were recovered as periodic;
  - TOI-5552.01 was recovered as an unexplained single transit;
  - TOI-1685 b was missed (now a known limit).

**Simplest generator first:** this maps directly onto transit-shape inversion.
- Fit the simplest model, a spherical body (planet or star) with standard limb darkening.
- Escalate only on structured residuals: oblate body, ring, disk, swarm, then non-natural geometry.
- This follows the published idea that a transit's shape carries information about the shape of the object in front.
- It is the trace-inversion pipeline adapted to this survey, and it can be built now.

---

## Target selection: archaeological priors

These lists choose *where* to apply the criteria. The work starts from the one confirmed site, us, and works outward.

| List | Source | Rationale | Approx. size |
|---|---|---|---|
| Solar-sibling candidates | Gaia chemistry plus orbit back-integration; published candidates | shared birth cluster, so shared early material | tens to low hundreds |
| Close-passage stars | Gaia astrometry back-integrated | physically plausible contact or probe paths | hundreds |
| Earth transit zone | published list (Kaltenegger & Faherty 2021) | stars that could have seen *our* transit | ~1,700 |
| Control sample | random stars matched on brightness and type | the same pipeline on non-prior targets | matched |

**The control sample is required.** Without it, the priors can't be shown to matter.

---

## Order of work

1. **WP0. Housekeeping.**
   - Move the 2023 notebooks to `msc-2023/`.
   - Move the brief and this programme into `docs/`.
   - Write the README.
   - First commit.
2. **WP1. Null library (C4).**
   - Catalogue the natural slow generators.
   - Extend the classifier to use them.
3. **WP2. Injection–recovery (C6).**
   - Measure detection efficiency on the Sector 59 sample.
   - Write the blind-spot statement.
4. **WP3. Targets.**
   - Build the three prior lists and the control sample.
   - Fetch their TESS light curves.
5. **WP4. Second observables (C3).**
   - Add the IR-excess and pixel-location columns.
   - Add chromaticity where archive data allows.
6. **WP5. Shape inversion (C7).**
   - Fit the simplest generator to every surviving dip.
   - Keep the residuals as the escalation signal.
7. **WP6. Ratchet (C2).**
   - Design the statistic.
   - Validate it on terrestrial series.
   - Then apply it to long-baseline stars.
8. **WP7. Report.**
   - The persistence table (C1).
   - The results and the blind spots.

WP1–WP3 can run alongside each other. WP6 depends on the theory side, which has the terrestrial validation series.

## Notes on the brief

- **Criterion 2 needs an operational guard against simple accumulation.** Monotonic natural processes are common, so the photometric test must require accumulation with recombination or deepening, not just one-way change. See C2.
- **On the open question about transmitters.** One possible wording: the *configuration* is deposited (it persists without the maker's motion), while its *emission* is driven. That keeps transmitters inside the deposited regime.
- **The references in the brief still need checking.** So do the ones this document mentions: Kaltenegger & Faherty 2021, the transit-shape papers, and DASCH. None has been verified here.
