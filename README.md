Overview
This Julia project analyzes Gottesman–Kitaev–Preskill (GKP) states entangled as a two-mode cluster state. It performs homodyne measurements in rotated phase space and studies the resulting state of the second mode. The analysis includes computing the quadrature q probability distribution, reconstructing the state on the Bloch sphere, determining its 
θ and ϕ parameters, and evaluating their corresponding probability distributions. The project simulates two-mode output states and produces comprehensive visualizations to support detailed quantum-state analysis.

Project Structure
run_01.jl - Main simulation script that orchestrates the full analysis pipeline
run_02_avg.jl - Secondary analysis script for averaging results
2mode_output_state.jl - Core functions
pdf_sampling.jl - Probability density function (PDF) calculations and sampling routines
utils.jl - Utility functions for plotting and data analysis
outputs - Generated simulation results and visualization outputs


-Usage
  Run the main simulation: julia run_01.jl
    This executes simulations for the parameter combinations specified and generates outputs in outputs/{timestamp}/ with plots and analysis data.
    --Output Files
      For each simulation (defined by rotation angle θ_r and β parameter), the script generates:
        *_pdf_q.png - Probability density function of q_homodyne measurement outcome
        *_pdf_theta.png - Probability density of θ
        *_pdf_phi.png - Probability density of ϕ
        *_pdf_2d.png - Joint 2D probability distribution
        *_fidelity_q.png - Fidelity vs. q
        *_maxima.txt - Local maxima coordinates and theoretical predictions
        *_params.txt - Simulation parameters
        
  Run the second simulation: julia run_02_avg.jl
    This executes simulation for a fix β parameter and averaging the joint probablity distribution over the specified range of rotation angle θ_r.
    
-Dependencies
  Julia (with plotting capabilities)
  Plots.jl - Visualization
  Dates.jl - Timestamping
