# ABM-Schelling
A MATLAB implementation of Schelling's Segregation Model

# Agent-Based Modeling of Schelling's Segregation Model

MATLAB implementation of Schelling's segregation model as an agent-based simulation on a fixed-boundary lattice, updated sequentially each iteration.

## Overview

Two agent groups ("Red" and "Green") occupy a square lattice; each agent checks its Moore neighborhood, and dissatisfied agents (below a desired percent-similar threshold) relocate to empty cells until the population converges or a 1000-iteration cap is reached. Default parameters replicate NetLogo v6.2's segregation demo [web:73].

## Requirements

- MATLAB R2020b
- No add-on toolboxes required beyond base MATLAB (`conv2`, `imagesc`, `spy`)

## Usage

```matlab
ABMSchellingSegregation
```

On first run, the script prompts for:

| Prompt | Description |
|---|---|
| ScrptPath | Directory where the script is located |
| ResPath | Directory where result figures and `.mat` files are saved |
| start | Starting simulation number (e.g., 1) |
| finish | Final simulation number to run (e.g., 10) |
| spyflag | Enter 1 to suppress visualization of unhappy agents, 0 to visualize |

The script re-invokes itself automatically to run sequential simulations from `start` to `finish`, saving results and clearing variables between runs.

## Model Parameters

| Parameter | Default | Description |
|---|---|---|
| N | 25 | Lattice side length (N x N grid) |
| density | 0.95 | Fraction of lattice occupied by agents |
| des_pct_sim | 0.50 | Desired minimum fraction of similar neighbors for satisfaction |
| msz | 8 | Marker size for plotting unhappy agents |
| Neighborhood | 1st-order Moore (8 neighbors) | Defined via a 3x3 convolution kernel |

## Output

For each simulation run, the script generates:
- PNG snapshots of the lattice at each iteration (`<SimNum>_SchellingSim_Iter<t>.png`)
- A plot of unhappy-agent count over time (`<SimNum>_EvolutionOfUnhappiness.png`)
- A `.mat` file with agent states, iteration count, and unhappiness history (`<SimNum>_SchellingSim.mat`)

## License

Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0).

## Author

© 2021 Joseph Blochberger — jblochb2@jhu.edu.

## Citations

> Blochberger, J. 2021. *Agent-Based Modeling of Schelling's Segregation Model.* https://github.com/joeblochberger/ABM-Schelling/ABMSchellingSegregation.m, GitHub. Retrieved Month Day, Year.

> Schelling, T., 1971. "Dynamic Models of Segregation." *Journal of Mathematical Sociology*, 1(2), pp. 143–186 [web:64].

> Wilensky, U. 1999. *NetLogo.* http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

![ABM Schelling Realization](https://github.com/joeblochberger/ABM-Schelling/blob/main/github_pdf.png)
