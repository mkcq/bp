# bp
Computer Science Bachelor Project 2023/2024.
Implementing Peachy assignment(s) in the Julia language.

## Peachy assignments
Peachy assignments can be found here:
- https://tcpp.cs.gsu.edu/curriculum/?q=peachy

The assignments that I have looked at are:
1. Introductory Python-based MPI Assignment (Mandelbrot) - EduPar 2023
2. Wind Tunnel Simulation - EduHPC 2021
3. Using MPI For Distributed Hyper-Parameter Optimization And Uncertainty Evaluation - EduHPC 2023
4. Agent-based simulation of fire extinguishing - EduHPC 2019

## Introductory Python-based MPI Assignment (Mandelbrot) - EduPar 2023
The assignment by Alina Lazar is mentioned in this paper:
- https://tcpp.cs.gsu.edu/curriculum/sites/default/files/peachy-edupar23.pdf
- https://tcpp.cs.gsu.edu/curriculum/?q=system/files/Lazar%20peachy.pdf

It's a nice introduction to MPI with a "trivially parallel" problem.

## Wind Tunnel Simulation - EduHPC 2021
The assignment by Arturo Gonzalez-Escribano and Yuri Torres is mentioned in this paper:
- http://tcpp.cs.gsu.edu/curriculum/sites/default/files/PeachyEduHPC21_0.pdf
- https://trasgo.infor.uva.es/sdm_downloads/wind-tunnel-peachy-assignment/

It looks like a good next step to MPI with a problem that resembles Successive Over Relaxation (SOR). Other keywoards such as load balancing, latency hiding and as the paper mentions "wave fronts" are also relevant.

## Using MPI For Distributed Hyper-Parameter Optimization And Uncertainty Evaluation - EduHPC 2023
The assignment by Erik Pautsch, John Li, Silvio Rizzi, George K. Thiruvathukal, and Maria Pantoja is mentioned in this paper:
- https://tcpp.cs.gsu.edu/curriculum/?q=system/files/PeachyEduHPC23.pdf
- https://figshare.com/articles/conference_contribution/Using_MPI_For_Distributed_Hyper-Parameter_Optimization_And_Uncertainty_Evaluation/24549673/1

Didn't continue due to the fact that the paper mentions that "The main Distributed Computing concept covered is how to divide a problem into chunks that can be easily assigned as independent tasks to nodes in MPI". Despite the interesting topic it unfortunately resembles the first Peachy assignment as a "trivially parallel" problem.

## Agent-based simulation of fire extinguishing - EduHPC 2019
The assignment by Arturo Gonzalez-Escribano and Jorge Fernandez-Fabeiro is mentioned in this paper:
- https://tcpp.cs.gsu.edu/curriculum/?q=system/files/Peachy_Eduhpc_19.pdf
- https://tcpp.cs.gsu.edu/curriculum/?q=system/files/fireAgents.zip (download)

Although an interesting topic, unfortunately at the time of writing (25-04-2024) all the necessary materials for the assignment were not included/could not be found.

## Citations/References
- Simon Byrne, Lucas C. Wilcox, and Valentin Churavy (2021) "MPI.jl: Julia bindings for the Message Passing Interface". JuliaCon Proceedings, 1(1), 68, doi: 10.21105/jcon.00068
- H. Bal, et al.,"A Medium-Scale Distributed System for Computer Science Research: Infrastructure for the Long Term" in Computer, vol. 49, no. 05, pp. 54-63, 2016.
doi: 10.1109/MC.2016.127
keywords: {programming;protocols;supercomputers;peer-to-peer computing;optical fiber networks;image processing}
Abstract: The Dutch Advanced School for Computing and Imaging has built five generations of a 200-node distributed system over nearly two decades while remaining aligned with the shifting computer science research agenda. The system has supported years of award-winning research, underlining the benefits of investing in a smaller-scale, tailored design.
url: https://doi.ieeecomputersociety.org/10.1109/MC.2016.127
- https://www.francescverdugo.com/XM_40017/dev/