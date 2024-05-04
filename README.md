# Trend Filtering for Functional Data


## Overview
This repository contains the R implementation of locally adaptive smoothing methods for functional data, as proposed in the following paper:

**Wakayama, T., & Sugasawa, S. (2023). [Trend filtering for functional data.](https://doi.org/10.1002/sta4.590)**

## Description
The repository includes the following files:

- `functions/` : Scripts implementing the proposed methods and related functions
  - `FTF-function.R`  : Functions for functional trend filtering (FTF)
  - `FHP-function.R`  : Functions for functional HP filtering
  - `FTFG-function.R` : Functions for functional trend filtering on graphs
  - `FHPG-function.R` : Functions for functional HP filtering on graphs
  - `FDPC-select.R`   : Function for tuning parameter selection

- `simulation.R` : Script for applying the proposed methods to generated datasets 

- `COVID19_cases/` : Script for applying the proposed methods to COVID-19 datasets
  - `AdjacencyMatrix_JPN.csv` : Adjacency matrix of prefectures of Japan
  - `CovidPerPopulation.csv`  : COVID-19 cases per one million population by prefecture
  - `Covid_Trend_Estimate.R`  : Trend estimation of transition of COVID-19 cases

## Demo
The 2nd order FTF is capable of localizing its estimates around strong inhomogeneous spikes, indicating its effectiveness in detecting events or spots of interest.

![covid_locally0-1](https://user-images.githubusercontent.com/44727480/127317873-f1d9c418-548b-426e-9d6d-c10aad01e9be.png)

