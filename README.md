# Modeling-based strategies for ICU-bed management during the COVID-19 pandemic

Geunsoo Jang, Su Min Park, Yongin Choi, Sung June Choi, and Hyojung Lee

## Abstract
During the COVID-19 pandemic, fluctuations in daily cases and severe infections posed substantial challenges to healthcare capacity, particularly in intensive care units (ICUs). This period underscored the need for data-driven decision-making systems to efficiently manage the limited healthcare resources. Therefore, to mitigate ICU overload during pandemic surges, we evaluated strategies to tighten non-pharmaceutical interventions (NPI) and expand ICU-bed capacity. We developed an age-structured transmission model incorporating ICU occupancy rate (IOR)-based activation thresholds, decision-to-action delays, and combined interventions targeting both demand (NPI intensity) and supply (capacity expansion). Our results indicated that during a seven-day delay in implementing control measures, maintaining the IOR below a high-occupancy benchmark of 0.8 requires escalation to at least NPI Level 3 when activation occurs before the IOR exceeds 0.5. When NPIs are strengthened with ICU-bed expansion, maintaining an IOR below 0.8 becomes feasible under milder NPIs (e.g., NPI Level 1), provided activation remains timely (before the IOR exceeds 0.6) and capacity increases by at least 10%. Escalation to NPI Level 3 further reduces the required expansion and yields control that is more stable across conditions. These findings provide a data-driven basis for ICU resource planning and policy design during health crises. The framework can inform the assessment of adaptive intervention strategies and support evidence-based healthcare resource planning in future health emergencies beyond the COVID-19 context

Link to the paper: 

## Code Description
This repository contains the code and data for the paper "Modeling-based strategies for ICU-bed management during the COVID-19 pandemic".

### 1. Data and Parameters Overview

The project relies on specific data and parameters organized as follows:

* **`data/` folder**: This directory contains sample input datasets for the simulations and key parameters used in the model.
    * `covid-19 cases & death`: Historical data on COVID-19 cases and deaths.
    * `population`: Population data relevant to the simulation.
    * `vaccine`: Vaccination data.
    * `beds`: Total bed data.
    * `current icu`: ICU data.
    * `severe_rate`: Fatality rate for severe cases.
    * `initial`: Initial state.
    * `unreported_rate`: Unreported rate by age-group.
  
### 2. Fitting Execution

The calibration and model fitting are executed via scripts located in the **`Fitting/` folder**:

* **`Fitting/fit_run.m`**: This script runs the fitting procedure and calibrates model parameters using observed ICU data/COVID-19 cases/COVID-19 death.
* **`Fitting/odef.m`**: This script contains the main ODE system
* **`Fitting/odeoperation.m`**: Implements the simulation workflow, performing operations such as parameter loading, preprocessing, and running the ODE solver for the fitting procedure
* **`Fitting/odeoperation_ode.m`**: A lower-level module that defines the ODE operations used inside the solver

### 3. Scenario Execution

The simulations are executed via scripts located in the **`scenario/` folder**:

* **`scenario/S1.m`**: This script runs the **NPI scenario** simulations.
* **`scenario/S2.m`**: This script runs the **comprehensive scenario** simulations.

### 4. Results Storage

All simulation outputs and results are saved in the **`result/` folder**.

