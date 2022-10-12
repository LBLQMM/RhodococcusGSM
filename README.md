# A genome scale model for *Rhodococcus opacus*

This repository contains notebooks that construct a genome scale model for *R. opacus* PD630 that includes 1773 genes, 2384 reactions, and 1581 metabolites. The model is then tested using constraint based reconstruction and analysis (COBRA) methods including flux balance analysis (FBA), parsimonious FBA (pFBA), E-Flux2, and SPOT (Simplified Pearson cOrrelation with Transcriptomic data). The model's predictions were tested against 13C-metabolic flux analysis data and measured growth rates.

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
- Notebook E: Determines FBA and pFBA solutions for phenol and glucose
- Notebook F: Calculates E-Flux2 solutions for glucose
- Notebook G: Calculates SPOT solutions for glucose
- Notebook H: Calculates E-Flux2 solutions for phenol
- Notebook I: Calculates SPOT solutions for phenol
- Notebook J: Calculates ATP maintance using FBA
- Under development -- Notebook K: Uses OptForce to predict gene knockouts

## Instructions for producing initial metabolic reconstruction

<ol>
  <li>download genome FASTA file from NCBI and save as r_opacus_bologna.faa from https://www.ncbi.nlm.nih.gov/assembly/GCF_020542785.1</li>
  <li>Run command to build model: carve r_opacus_bologna.faa -u grampos -o r_opacus_bologna_raw.xml</li>
  <li>Run command to gapfill model: gapfill r_opacus_bologna.xml -m M9,LB -o r_opacus_bologna_gapfilled.xml</li>
</ol>

The following versions were used to make [r_opacus_bologna_gapfilled.xml][1]:

<ul>
  <li>CarveMe 1.5.1</li>
  <li>Diamond 0.9.14</li>
  <li>CPLEX 12.10.0.0</li>
</ul>

## Reference
This work has not been been peer reviewed.

## License

This code is distributed under the 3-Clause BSD license specified in the [license][2] file. It is open source and commercially usable.

[1]: models/r_opacus_bologna_gapfilled.xml
[2]: license
