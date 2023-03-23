# Cellular-Level Computational Models of the Mouse Heart

These complex models describe the electrophysiology of mouse cardiac cells, and simulate heart beats at the cellular level

## Monorepo containing: 
- Atrial Myocyte Model, driven by an external stimulus
- Atrial Myocyte Model coupled to Sinoatrial Pacemaker Cell model

## Brief Technical Introduction

Central to most, developed forms of life lies the heart; modelling its function, to better understand the root causes of dysfunction, have become growing areas of research in recent years. Electrical processes within the heart drive mechanical contraction, and cardiac electrophysiology is the study of such processes. Electrophysiological models of the heart describe the passage of an action potential (AP) through various cardiac cells, for which the potential of a cell membrane is rapidly raised (depolarized) and sequentially lowered (repolarized) in relation to extracellular space. This is driven by a number of transmembrane ion currents, which vary in amplitude throughout the AP and depend on dynamic physiological cell parameters, such as localised ion concentrations and the membrane potential. These parameters are known as state variables; computational models aim to numerically solve the coupled, non-linear, typically stiff ODEs that describe the variation of state variables over time. At each integration step, individual ion currents are calculated, enabling a simulated cardiac AP to be traced.

Developed within my Master's project for the University of Manchester Physics Department


