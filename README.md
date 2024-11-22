# Python and MATLAB code and data for statistical physics analysis of social behavior in Eco-HAB mice.

Original codes for the manuscript:\
Modelling collective behavior in groups of mice housed under semi-naturalistic conditions \
Chen, Winiarski, Puścian, Knapska, Mora, Walczak (2023)


For original data and intermediate steps: \
The folder data contains original data.\
The folder processed_data contains data processed by the Jupyter notebook in python_code.\
The folder output contains output / results from the maximum entropy model learning, which are used to generate figures.


In the python_code folder, the Jupyter notebook shows how Fig.1 is generated.



In matlab_code folder:

1. The codes that generate intermediate results (e.g. inferred parameters of the maximum entropy model) from processed data is in the maxent subfolder. 
To fully operate the code, download UGM package (Mark Schmidt 2007) from https://www.cs.ubc.ca/~schmidtm/Software/UGM.html.

2. The codes that generates figures from intermediate results:


plot_fig2.mlx
generates all non-schematic panels of Figure 2, including Fig. 2C, 2D, 2E,
and the supplementary materials Fig.S2, S3, S4, S5.

plot_fig2f_fig2sup8.mlx
generates panel F of Figure 2, and its linked Supplementary Figure 8.

plot_fig3.mlx
generates all non-schematic panels of Figure 3, and its linked Supplementary Figures.

plot_fig4.mlx
generates all non-schematic panels of Figure 4, and its linked Supplementary Figures.

plot_figs1.mlx
generates supplementary figure Figure 1 - figure supplement 1.

plot_figs7.mlx
generates supplementary figure Fig.S7.

plot_figs8_s9.mlx
generate supplementary figures Fig.S8 and S9.

plot_figs12_mi.mlx
generates supplementary figure Fig.S12.

plot_figs13.mlx
generates supplementary figure Fig.S13.
