# 3rd-Year-Project
Backup for my 3rd year project

Brief Overview:

The Wall shear stress exerted on an artery due to blood flow has been shown to be an early indicator of cardiovascular disease. However, it is difficult to measure the WSS on an artery, current methods are too invasive and unreliable. The best way currently to measure the WSS on an artery is to use a CFD solver (after imaging the vessel and establishing boundary conditions). CFD solvers take a non-trivial amount of time to solve a flow case, the meshing process in particular can take some time and compute resources. The idea behind this project is to train a deep learning to model to predict the WSS on an artery based on a number of features (I have yet to determine which features are best, but currently they are (X,Y,Z) coordinates, velocity, tangent and normal surface vectors, distance from flow source, curvature). 

One significant problem that I am solving first is the availability of data, it is quite a time consuming task to image the required coronary arteries, as a result there are no large datasets that are publicly available that have the required data. However, I was given a small dataset of around 50 geometries as part of my 3rd Year Project, from this I developed a Statistical Shape Model in MATLAB. I can now generate statistically valid geometries. I also have code for a batch solver on MATLAB that calls ANSYS CFX.

The next step is to retrain an existing U-Net architecture that I have separately.

I am also looking at other types of architectures such as PINNs and Mesh Graph Nets that may increase accuracy.
