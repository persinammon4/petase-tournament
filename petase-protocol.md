# PETase Expression and Activity Protocol

## Cloning

Golden gate assembly was used to clone DNA encoding CDS of each PETase sequence into PET28a (KanR) plasmid in the insert region as annotated in this genebank file. Sequence was verified using whole plasmid sequencing (Quintara Biosciences).

## Expression

BL21(DE3) cells harboring plasmid DNA with each ORF were pre-cultured in 1 mL LB (Teknova) containing 50 µg/mL kanamycin and 0.5% glucose overnight at 37°C with shaking at 300 rpm in a 96-deep well block. The next day, 15 µL of cells per construct were sub-cultured in 1 mL of TB medium (Teknova) containing 50 µg/mL kanamycin, were grown at 37°C while shaking at 300 rpm until mid-log phase (approximately 2-3 hours, OD 0.6-0.8), induced with 0.5 mM IPTG. Post-induction, the growth temperature was lowered to 18°C overnight. Cells were harvested by centrifugation at 4,000 x g and pellets were frozen at -80°C before use.

## Purification

Harvested cell pellets were incubated for one hour with gentle shaking at room temperature with 1 mL of lysis buffer (20 mM sodium phosphate pH 7.4, 500 mM NaCl, 20 mM imidazole, 1 U/mL Benzonase, PI tablets, and 0.1 mg/mL lysozyme) and lysates are clarified by centrifugation at 4,000 x g for 15 min at 4°C.

For small-scale purification, 1 mL of soluble lysate was used for IMAC pulldown with KingFisher APEX. Magnetic resin (50 µL) was washed 3 × 1 mL with 20 mM sodium phosphate pH 7.4, 500 mM NaCl, and 20 mM imidazole, and then protein was eluted in 100 µL of 20 mM sodium phosphate pH 7.4, 500 mM NaCl, and 300 mM imidazole.

Protein concentration of purified proteins was determined by the average of two replicates in an A660 assay (Pierce) using a BSA standard curve ranging from 50-2,000 µg/mL, and was measured spectrophotometrically on a BMG LabTech SPECTROstar Nano.

## Specific Activity Assay

### Sample Preparation

Purified protein in the storage buffer containing 20 mM sodium phosphate pH 7.4, 500 mM NaCl, and 300 mM imidazole was placed on wet ice until completely thawed. 4 µL of purified protein was added to a 96-well Costar v-bottom plate containing 196 µL of 5 mg/mL powdered PET substrate (Goodfellow ES30-PD-000132) resuspended in buffer using a Perkin-Elmer Janus Mini pipetting station.

### Reaction Conditions

| Condition | Buffer | pH | Temperature |
|-----------|--------|-----|-------------|
| 1 | 50 mM Citrate | 5.5 | 30°C |
| 2 | 50 mM Glycine | 9.0 | 30°C |

### Incubation

Plates were sealed and shaken at 450 rpm within a temperature-controlled plate shaker (Allsheng Thermo Shaker) at 30°C for two, six, and eight hours.

### Reaction Termination

Reaction was stopped by incubating on wet ice for 5 minutes followed by centrifugation at 4,000 rpm for 5 minutes at 4°C to pellet the PET substrate. 50 µL of supernatant was transferred to a 96-well Costar v-bottom plate containing 100 µL of crash solution (0.1% formic acid in acetonitrile) containing 300 µM d4-TPA (internal mass spectroscopy standard). Plates were sealed with pierceable plate seal and shaken at 450 rpm at room temperature for 5 minutes followed by centrifugation at 4,000 rpm for 5 minutes at room temperature.

### Mass Spectrometry Analysis

Enzyme product was determined using an Applied Biosystems RapidFire / Sciex 4000 QQQ Mass Spectrometer. A 10 µL sample of each well was injected into the mass spectrometer and the peak areas of TPA, MHET, and the internal standard were determined. Included with each plate of samples was a blank, negative, and positive control. Additionally, a standard curve of TPA and MHET was included with each set of plates. Samples were tested in duplicate.

## Data Analysis

Peak areas were normalized using the internal standard d4-TPA to correct for any differences incurred during sampling. Normalized peak area for each product was interpolated from a standard curve of TPA and MHET to determine the molar equivalent generated of the enzyme product. Products generated at each time point were analyzed by linear regression to calculate a reaction rate. The amount of protein in each sample (previously determined using a standard BCA assay) and the calculated reaction rate were used to determine the specific activity in **µmoles product/min/mg enzyme**.
