# Locally Adaptive Smoothing For Functional Data


## Overview
This repository provides R code implementing functional trend filtering and other smooothing methods, as proposed by the following paper. 

Wakayama, T. and Sugasawa, S. (2021). Locally Adaptive Smoothing for Functional Data. https://arxiv.org/abs/2104.02456


## Description
The repository includes the following files.

- functions : Script implementing the proposed methods and the related
  - FTF-function.R  : functions for functional trend iltering, including FTF(method) and FTF_select(parameter tuning)
  - FHP-function.R  : functions for functional HP filtering, including FHP(method) and FHP_select(parameter tuning)
  - FTFG-function.R : functions for functional trend Filtering on graph
  - FHPG-function.R : functions for functional HP filtering on graph
  - FDPC-select.R   : function for tuning parameter of FDPC

- simulation.R : Script for applying the propsed time-series methods to generated datasets 

- COVID19_cases : Script for applying the proposed spatial methods to the COVID-19 datasets
  -  AdjacencyMatrix_JPN.csv : Adjacency matrix of prefectures of Japan
  -  CovidPerPopulation.csv  : COVID-19 cases per one million population by prefecture
  -  Covid_Trend_Estimate.R  : trend estimation of transition of COVID-19 cases


## Demo
2nd order FTF (functional trend filtering) is able to localize its estimates around strong inhomogeneous spikes, which implies that it is able to detect the event or spot of interest.


![covid_locally0-1](https://user-images.githubusercontent.com/44727480/127317873-f1d9c418-548b-426e-9d6d-c10aad01e9be.png)

