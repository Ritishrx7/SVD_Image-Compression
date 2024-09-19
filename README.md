# Simulation in R  
*Master 1 in Modelisation Statistique Project*  
Université Franche-Comté, Besançon (2022–2023)  

## **Project Overview**  
This project simulates the random walk of a turtle on a 2D grid using **Monte Carlo methods**. The goal is to estimate the distribution of the number of revisits \((N_n)\) to previously visited points for different step sizes \((n = 100)\), \((n = 1000)\), and \((n = 10000)\).

### **Authors**:
- A.V Hurtado Quiceno
- Narain Ritish

## **Problem Description**  
The turtle starts at the origin \((0,0)\) and moves randomly in one of four directions: \((0,1)\), \((1,0)\), \((0,-1)\), or \((-1,0)\). After each step, we record whether the turtle revisits a previously visited point.

For each value of \(n\), we calculate \(N_n\), the number of revisits, and use Monte Carlo simulations to analyze the distribution of \(N_n\).

## **Key Steps**:
1. **Simulate the random walk** for different values of \(n\).
2. **Use Monte Carlo methods** to estimate the distribution of revisits.
3. **Visualize the results** using histograms and plots.

## **Key Aspects**:
- **Objective**: To estimate the distribution of revisit events using Monte Carlo simulations.
- **Simulation**: Conducted for \(n = 100\), \(n = 1000\), and \(n = 10000\) steps.
- **Monte Carlo Methods**: Used to run multiple simulations and generate a graphical representation of \(N_n\).
- **Graphical Analysis**: Histograms and density plots were created to understand the distribution of revisit events.

## **R Packages Used**:
- `ggplot2`: For creating visualizations like histograms and density plots.
- `dplyr`: For data manipulation and summarization.
- `purrr`: For mapping functions and efficiently running multiple simulations.
- `tidyr`: To organize and reshape the simulation data.
- `tibble`: For creating and managing tidy data frames.
- `gridExtra`: For arranging multiple plots in a grid format for comparison.

## **Project Files**:
- `SimulationProject_Ritish_Andrea.Rmd`: R Markdown file containing the full analysis and code.
- `index.html`: Rendered HTML report showing results of the analysis.
