# A genome scale model for Rhodococcus opacus

This repository contains notebooks that construct a genome scale model for R. opacus PD630 that includes 1773 genes, 2384 reactions, and 1581 metabolites. The model is then tested using constraint based reconstruction and analysis (COBRA) methods including flux balance analysis (FBA), parsimonious FBA (pFBA), E-Flux2, and SPOT (Simplified Pearson cOrrelation with Transcriptomic data). The model's predictions were tested against 13C-metabolic flux analysis data and measured growth rates.

- [System Requirements](#system-requirements)
- [Instructions for use](#instructions-for-use)
- [Summary of notebooks](#summary-of-notebooks)
- [Reference](#reference)
- [License](#license)

## System Requirements

The code was written using python 3.7

## Instructions for use

To run the code in this repository use the following commands:

<ol>
  <li>git clone https://github.com/garrettroell/OpacusBiodesignGR.git</li>
  <li>python3 -m venv venv</li>
  <li>source venv/bin/activate</li>
  <li>pip install -r requirements.txt</li>
</ol>

## Summary of notebooks

- Notebook A: Adds uptake reactions to the genome scale model
- Notebook B: Adds custom biomass reactions to the model
- Notebook C: Annotates and curates the model for MEMOTE testing
- Notebook D: Calculates experimental growth parameters
- Notebook E: Runs FBA and pFBA for both carbon sources
- Notebook F: Runs E-Flux2 with transcripts from glucose conditions
- Notebook G: Runs SPOT with transcripts from glucose conditions
- Notebook H: Runs E-Flux2 with transcripts from phenol conditions
- Notebook I: Runs SPOT with transcripts from phenol conditions
- Notebook J: Calculates ATP maintance using FBA
- Notebook K: Calculates ATP maintance using FBA

## Reference
This work has not been been peer reviewed.

## License

This code is distributed under the 3-Clause BSD license specified in the [license][1] file. It is open source and commercially usable.

---

[1]: license
