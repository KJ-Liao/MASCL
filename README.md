# MASCL (Molecular Assembly Simulation in Crystal Lattice)

## Introduction

The MASCL Protocol and AAI-PatchBag framework were developed to facilitate accurate analysis of protein structural features and provide deeper insights into protein–protein interactions. This project demonstrates how to:

- Use the PIPER docking program, along with MATLAB scripts that handle symmetry, to predict crystal packing for a given protein.
- Assess potential protein crystallization conditions by analyzing similarities of protein packing interfaces using AAI-PatchBag.

## Contents

Detailed explanations of each part (including the working principles, dataset descriptions, and script execution steps) are available in the respective README.txt files.

- ### Dataset
  Stores the datasets used in this study, including PDB entries and associated structural information.

- ### Crystallization_Condition_Parsing
  Contains the datasets, scripts, and related files for parsing crystallization conditions, along with preliminary parsing results.

- ### MASCL
  Explains how to generate pairwise assemblies using AF3 or AF2 + C2-DIPER, then perform additional docking to obtain tetrameric complexes, check symmetry, and ultimately derive the minimal protein crystal packing model. Includes scripts and examples.

- ### AAI-PatchBag
  Describes how to extract and segment molecular surfaces, build a PatchBag Library, and expand it into an AAI-PatchBag by applying physicochemical properties (AAI-index). Provides the necessary scripts and examples.

- ### Xtal_CondiRanker
  Offers a simple script that combines sequence similarity, RMSD, 3DZD, molecular weight, and other conventional protein features with AAI-PatchBag descriptors. This script compares the target protein against a pre-built database to rank similar PDB entries, which can inform crystallization strategies.

## Required External Tools

- PIPER    : physical-based docking program.                                <br>
   available in: https://cluspro.org/downloads.php
- EDTSurf  : rapid macromolecular surface calculation.                      <br>
   available in: https://zhanggroup.org/EDTSurf/
- generate_3dzd.py  : molecular shape in ply format and 3DZD computation.   <br>
   available in: https://github.com/kiharalab/3d-af_surfer
- MakeShape, Shape2Zernike	: 3DZD computation.                              <br>
   available in: https://github.com/jerhoud/zernike3d/tree/main
