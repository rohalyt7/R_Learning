# **R Framework for Human Subjects Research**
**Author:** Thomas Rohaly

**Role:** Computational Research Psychologist

**Status:** Active Development

## **Overview**
This repository provides a modular ETL (Extract, Transform, Load) and Analysis pipeline designed for the behavioral and physiological sciences. The core philosophy is to provide a standardized toolkit that guides researchers from raw data ingestion and cleaning through to advanced statistical modeling and reproducible result generation.

## Core Philosophy: Simulation-Driven Rigor
Unlike traditional tutorials using static datasets, this project emphasizes Data Simulation as a prerequisite for robust analysis.

- **Tool-Based Simulation:** Early modules utilize the *simstudy* package and Base R to introduce fundamental data generation.

- **Custom Research Frameworks:** Advanced modules leverage custom-built functions to generate highly controlled, clustered, and hierarchical data. This allows for the testing of models against known ground truths before applying them to "noisy" real-world physiological streams.

For the source code behind these custom simulation frameworks, see the MultilevelPowerSim Repository.

## **Current Modules**
**1. Data Manipulation & Basic Simulation**
- **Workflow Comparison:** A side-by-side evaluation of Base R vs. Tidyverse syntax for data cleaning and transformation.

- **ETL Essentials:** Techniques for data reshaping (Wide vs. Long formats) and cleaning to prepare raw inputs for the analytic pipeline.

- **Basic Simulation:** Introduction to generating synthetic research data using simstudy.

**2. Exploratory Data Analysis (EDA)**
- Visualizing distributions, trends, and outliers using ggplot2 to ensure data integrity prior to modeling.

**3. OLS Assumption Diagnostics**
A comprehensive suite for verifying the mathematical foundations of linear modeling, including:

- **Normality & Homoscedasticity:** Visual and statistical checks for residual distribution.

- **Independence of Observations:** Critical for non-clustered OLS applications.

- **Multicollinearity:** Diagnostic checks (VIF/Correlation matrices) to ensure feature independence and model stability.

## Roadmap
- **Inference:** T-tests and ANOVA within the Tidyverse framework.
- **Regression Modeling Suite:** Simple/Multiple linear, Quadratic, and Logistic models, with a focus on robust estimation methods.
- **Multilevel Modeling:** Linear Mixed-Effects Models (LMMs) for clustered/longitudinal data.
- **Bayesian Frameworks:** Transitioning to brms for robust hierarchical modeling.
- **Bio-Signal ETL:** Specialized pipelines for handling missing data and sensor "dropouts" in wearable research.
- **Advanced Data Simulation:** A deep dive into building custom simulation frameworks for clustered and longitudinal data to validate model power and sensitivity before data collection.

## **License**
This project is licensed under the MIT License - see the LICENSE file for details.
