# Methods reporting checklist for quantitative ISH

There is no reporting standard for quantitative in situ hybridisation. The
nearest usable substitutes are the Schmied 2024 image-analysis checklist,
QUAREP-LiMi and ARRIVE 2.0. This checklist follows those.

`write_methods_stub()` fills most of it from the config into
`METHODS_STUB.md` in the results folder. Anything it cannot fill appears as
`TODO` rather than being quietly omitted.

---

## 1. Tissue and pre-analytics

- [ ] Species, strain, sex, age, n per group
- [ ] Fixation: fixative, duration, temperature
- [ ] FFPE or frozen; time to fixation
- [ ] **Section thickness in um.** A first-order multiplier on puncta per cell,
      and not corrected for anywhere in this field. Published work spans a
      four- to fivefold range, so counts are not comparable across studies or
      sample formats at different thickness, and no Abercrombie-type or disector
      correction is applied
- [ ] Storage time and conditions of cut sections
- [ ] RNA quality if assessed

## 2. Assay

- [ ] Kit, catalogue number, **lot**
- [ ] Probe catalogue numbers and target region
- [ ] For a compound probe: **which strand** (siRNA passenger versus guide; ASO)
- [ ] Protease type, concentration, duration; retrieval duration and
      temperature. Where a protocol number is load-bearing, confirm it against
      the manual for the kit lot in use, as these differ between revisions
- [ ] Chromogen or fluorophore per channel, amplification dilution, incubation
      time
- [ ] **Positive control probe** and the result
- [ ] **Negative control probe** and the result
- [ ] For co-detection: antibody, clone, dilution, and whether it survived the
      ISH protocol. Antibody survival through a protease-based protocol is poor
      for most antibodies, and membrane epitopes fare better than soluble ones.
      State whether an IHC-only serial section was run
- [ ] For fluorescence in autofluorescent tissue: how autofluorescence was
      handled (quenched, masked as its own channel, unmixed, or not), and which
      channel each probe occupies. Lipofuscin is worst in green, so a rare
      target should never be placed there in aged tissue

## 3. Acquisition

- [ ] Instrument and objective. **Single-dot resolution needs 40-63x**; 20x
      gives regional context and intensity, not countable dots
- [ ] Pixel size in um, taken from the image metadata rather than assumed
- [ ] **z handling**: single plane, z-stack (number of planes, step size), or
      extended-focus projection. Modality matters most where measurement is
      hardest: confocal and widefield can disagree on the lowest-expressing
      target while agreeing on higher ones. Thin FFPE sections are largely
      adequate in a single plane; thick frozen sections undercount, and
      out-of-focus dots inflate background. Extended-focus projections make dots
      visible but do not preserve axial separation, so they undercount in dense
      cells
- [ ] Exposure and gain **fixed once and recorded**, not re-autoexposed per slide
- [ ] For fluorescence: unmixing library composition. It needs a DAPI-only slide
      **and an unstained sample from the same study**, so autofluorescence
      becomes its own endmember
- [ ] Batch structure. If the study spans batches, was the order randomised (one
      control, one treated, alternating)? Was a control animal included in every
      slide batch to anchor the background threshold?

## 4. Image analysis: the free parameters

Every one of these silently determines the result.

- [ ] Software and **version**
- [ ] Segmentation method. Built-in cell detection in general-purpose tools
      underperforms dedicated segmentation models on demanding tissue
- [ ] **Cell expansion radius** used to assign puncta to nuclei
- [ ] **Minimum spot size**
- [ ] **Spot intensity threshold**
- [ ] **Probe copy intensity threshold**, where applicable
- [ ] **The single-dot calibration used for cluster deconvolution.** An
      "estimated spots" column is cluster area divided by an assumed spot area;
      a "copies" column is cluster intensity divided by a single-copy intensity
      calibrated on isolated dots. Without the divisor the numbers are not
      comparable between slides
- [ ] Colour deconvolution or unmixing vectors: fixed across the batch, or
      estimated per image? Per-image estimation moves the optical-density scale
      from slide to slide
- [ ] **Whether ROIs were operator-selected, and whether the selector was
      blinded**
- [ ] Atlas registration if used

Export the per-cell table and do all thresholding, classification, normalisation
and inference in code under version control. The GUI is a feature extractor. The
commonest reproducibility failure in this field is that the decisive parameter
tuning happened inside a GUI and was never written down.

## 5. Readout

- [ ] **Which expression regime the data are in**, per stratum. Sparse, up to
      roughly 10-15 puncta per cell, supports a genuine count. Dense supports an
      index only. Saturated supports integrated intensity or per cent labelled
      area
- [ ] Primary readout, declared before analysis
- [ ] **Positivity threshold and where it came from.** The rule is mean plus one
      SD of puncta on negative-control tissue, counting only non-zero values
- [ ] Scoring scale **defined explicitly**. Scales for this assay are not
      standardised and published vendor documentation is internally inconsistent
      about them, so the definition has to travel with the result
- [ ] H-score range and how the bins were defined. Bin boundaries are roughly
      geometric, so the score behaves like a log summary: differences are not
      proportional to abundance, and bin 0 contributes nothing, so two tissues
      with very different fractions of expressing cells can score identically

## 6. Denominator

- [ ] Which denominator, and why
  - **Per cell** — for cell-autonomous questions. Depends entirely on
    segmentation; cell-density differences between groups become apparent
    expression differences
  - **Per mm2** — for neuropil and unsegmentable tissue. Confounded by
    cellularity: oedema, infiltration and atrophy all move it. **Always report
    cells/mm2 alongside**
  - **Per nuclear area** — rarely right. Nuclear area is not cytoplasmic volume
    and the signal is largely cytoplasmic
- [ ] Cells and cells/mm2 per animal per region
- [ ] **Do not normalise to a housekeeping probe.** Degradation hits
      high-expressing housekeepers hardest and low-to-moderate expressors least,
      so the control's response to a pre-analytical insult is not proportional to
      the target's and the ratio adds noise instead of cancelling artefact. Use
      one as a sample-level inclusion gate instead

## 7. Statistics

- [ ] **What one row of the analysed data frame is.** State it explicitly
- [ ] Random effect for animal, and section within animal
- [ ] Distribution assumed, and the comparison that justified it
- [ ] Both hurdle components reported: fraction of positive cells, and counts
      among positive cells
- [ ] `offset(log(exposure))`, never a pre-computed ratio
- [ ] **Intraclass correlation reported**
- [ ] No `log(x + 1)` on counts; a log link on untransformed counts
- [ ] Pseudobulk sanity check
- [ ] Multiplicity: BH-FDR across the declared grid, fixed before analysis
- [ ] Exclusions listed with reasons that are image-based and predate looking at
      the numbers

## 8. Interpretation limits for compound-detection work

- [ ] Presence and pattern reported, not amount
- [ ] No claim that a negative region means the compound is absent; no limit of
      detection has been published in any units
- [ ] No concentration language; no correlation of ISH signal against LC-MS or
      hybridisation-ELISA tissue concentration has been published
- [ ] No claim to distinguish full-length compound from n-1, n-2 or
      chain-shortened metabolites
- [ ] Detected compound is not active compound: two uptake pathways exist, one
      ending in the endolysosome where the oligonucleotide is inert
