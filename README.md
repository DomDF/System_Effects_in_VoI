# System Effects in Identifying Risk-Optimal Data Requirements for Digital Twins of Structures

This repository is associated with a research paper that has been submitted to the Journal of Reliability Engineering & System Safety.
The preprint is available here: link.

Within the "analysis" folder, the file "system_effects_in_VoI.jl" is a Julia file that contains all of the calculations used to produce the results in the paper.
Within the "results_files" folder, are the simulations used to estimate the expected value of various combinations of data collection activities. This is the raw data used to produce the visualisations in the paper.
 
Abstract:

Structural Health Monitoring (SHM) technologies offer much promise to the risk management of the built environment, and they are justifiably an active area of research. SHM provides opportunities to collect large amounts of data about the local environment, including measurements of strains at key locations. However, this data alone (even when coupled with the state of the art in statistical analysis) is only related to one aspect of structural integrity management. Information regarding material properties, such as toughness and strength, cannot presently be monitored using sensors, but are measured in destructive lab tests. Similarly, the presence or emergence of geometrical anomalies require inspection for detection and sizing.

Value of Information (VoI) analysis is an established statistical framework for quantifying the utility of a prospective data collection activity in the context of solving an underlying decision problem. VoI can be used to identify how much engineering teams should be willing to pay for material testing, inspection and SHM to support the risk management of the built environment. In this paper, the expected value of various combinations of these different data types are quantified, demonstrating that dependencies between different elements of a structural integrity management system dictate expected optimal data collection planning. In summary, system-level decision making, requires system-level models.

Calculations were completed to provide decision support for the integrity management of bridges. Multi-variate Bayesian models were used with influence diagrams, and the associating Julia code has been made available. Although specific results are case-dependent (on the prior uncertainties in the system, and utility functions), sensitivity analysis and wider discussion considers the implications more broadly.
