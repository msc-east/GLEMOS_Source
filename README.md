## Description

This repository is a part of the GLEMOS WorkSpace and contains the sourse code for the GLEMOS model including scientific routines for modeling mercury (Hg), particle-bound heavy metals (Cd, Pb) and persistent organic pollutants (POPs). The current stable release is 2.2.2.

## The GLEMOS User Manual

A copy of the GLEMOS User Manual can be found in the root directory of each of the four GLEMOS repositories.

## About GLEMOS

Global modeling framework GLEMOS is a multi-scale multi-pollutant simulation platform developed for operational and research applications within the EMEP programme [Tarrason and Gusev, 2008; Travnikov et al., 2009, Jonson and Travnikov, 2010, Travnikov and Jonson, 2011]. The framework allows simulations of dispersion and cycling of different classes of pollutants (e.g. heavy metals and persistent organic pollutants) in the environment with a flexible choice of the simulation domain (from global to local scale) and spatial resolution. In addition, GLEMOS supports multi-media description of the pollutants cycling in the environment. A modular architecture of the modeling system allows flexible configuration of the model set-up for particular research tasks and pollutant properties. More infromation about the model can be found on the MSC-E website https://www.msceast.org/j-stuff/glemos.

## Model evaluation

The GLEMOS modelling system was extensively evaluated in a number of numerical experiments and multi-model studies within the Task Force on Hemispheric Transport of Air Pollution (TF HTAP). The validation program included testing the atmospheric transport, evaluation of model performance against observations and assessment of source attribution abilities on a global scale. In addition, the atmospheric transport module of GLEMOS was recently tested in a numerical experiment based on dispersion of radioactive isotopes from the Fukushima-1 accident. The model performance in simulation of Hg pollution on a global scale was tested in the multi-model assessments within the Global Mercury Observation System (GMOS) project [Travnikov et al., 2017] and the Global Mercury Assessment 2018 [AMAP/UN Environment, 2019].

## System Requirement 

The model has been verified to run on Linux Ubuntu 20.04 with HDF5 v1.10.5, netcdf-c v4.7.2, and netcdf-fortran v4.5.2 libraries preinstalled. The model uses from up to 34 GB of RAM, depending on the pollutant and the type of calculation. The matrix calculations require much more RAM than field ones, and Hg uses more RAM than POPs, while POPs uses more than HM). Currently the model is capable only with the f95 compiler that is provided with Oracle Developer Studio 12.6. Besides, the model does not support multithreading (works only in a single thread mode).

## Downloading the input data

The model input data are separated into several storages. Each piece of data is reachable by corresponding link:

* Emissions - https://drive.google.com/file/d/1UlbO98wB_J9sKhq2kuIPt5iDHwY40Jkm/view?usp=sharing
* MeteoData - https://drive.google.com/file/d/1K06uLzNB0ZHmt0VbkzDoxuAcP5lyblwy/view?usp=sharing
* Dust - https://drive.google.com/file/d/1OM82xTY_q6YIqNu6g0uuAKQ5kjbgZsaH/view?usp=sharing
* InitCond - https://drive.google.com/file/d/1HPA7rM-iD0qUD7CHJj4Gavp_4u32xp_L/view?usp=sharing
* LandCover - https://drive.google.com/file/d/1Kpp62GMpdVXWgFJ85JRRPirfy50t4dLC/view?usp=sharing
* ReactData - https://drive.google.com/file/d/1gsn5z7J77ReN_GfAaL-h5IjqIJjeymKL/view?usp=sharing

This dataset covers a short time period (January 2020) and is an example of input data for the model. It includes data for monthly modeling of Hg, Cd, Pb and BaP on global 1°x1° and on regional 04°x04° scales (the regional domain covers the EMEP region). Detailed description of these data as well as instructions for organizing local storage can be found in the chapter 5 of the GLEMOS User Manual. 

## Configuring of model runs

Simulation management and configuration as well as the desctiption of main features can be found in sections 6 and 1 of the GLEMOS User Manual, respectively.

## The output structure

The model supports output in netcdf4 and txt formats. More information can be found in the GLEMOS User Manual.
