# Cell-Cell Interaction in Prostate Cancer 

A computational-mathematical multiscale model. The model requires a set of partial differential equations (PDE), a set life-cycle flowcharts and a cellular automaton framework to integrate them (1).

Reference:
1. A Unique s•	Stromal reactivity differentially drives tumour cell evolution and prostate cancer progression. Frankenstein Z, Basanta D, Franco OE, Gao Y, Javier RA, Strand DW, Lee M, Hayward SW, Ayala G, Anderson A. Nature Ecology & Evolution. 2020 May 11.
 
Abstract

Prostate cancer (PCa) progression is a complex eco-evolutionary process driven by the feedback between evolving tumour cell phenotypes and microenvironmentally driven selection. To better understand this relationship, we used a multiscale mathematical model that integrates data from biology and pathology on the microenvironmental regulation of PCa cell behaviour. Our data indicate that the interactions between tumour cells and their environment shape the evolutionary dynamics of PCa cells and explain overall tumour aggressiveness. A key environmental determinant of this aggressiveness is the stromal ecology, which can be either inhibitory, highly reactive (supportive) or non-reactive (neutral). Our results show that stromal ecology correlates directly with tumour growth but inversely modulates tumour evolution. This suggests that aggressive, environmentally independent PCa may be a result of poor stromal ecology, supporting the concept that purely tumour epithelium-centric metrics of aggressiveness may be incomplete and that incorporating markers of stromal ecology would improve prognosis.

![image](https://github.com/user-attachments/assets/69476264-6e65-4f3f-ab78-7f506cb9aa87)

a, Interaction network of key model variables. Interactions between cells (coloured nodes) and microenvironmental variables (lilac nodes) are represented as either green (positive) or red (negative) connections. Multicoloured connectivity represents the spectrum of possible tumour phenotypes with different levels of growth factor and MMP production. Bicoloured connectivity represents two different degrees of stromal reactivity. b–d, In silico reconstruction of the normal prostate peripheral zone tissue. b, Histopathological slide of the whole normal prostate, highlighting the peripheral zone, filled with epithelial acini surrounded by stroma (magenta). c, In silico representation of the complete peripheral zone, including ductal structures and cellular densities that mimic normal anatomy. This constitutes the domain where all simulations were performed. The inset on the bottom left is an example of a sample simulation d, Representation of a single reconstructed duct and the surrounding stroma, as well as the total number of cell types. e, Cell decision flow charts for each cell type in the model. The phenotypic behaviour of an individual cell is based on the interaction between the cell and the local microenvironment. GF, growth factor.



Method

Multiscale prostate peripheral zone model

The model we developed builds on an HCA model. The definition of an HCA model requires a set of partial differential equations that characterize the physical microenvironment, a set of life-cycle flow charts that characterize the behaviour of the cells under microenvironmental constraints and a cellular automaton framework to integrate them. The following system of nonlinear partial differential equations define growth factor (G), MMP (E) and ECM/basement membrane (M) as the three key continuous microenvironmental variables:

![image](https://github.com/user-attachments/assets/805849fe-52d9-46f4-b0f7-465b6c5d1e3b)


δG and δE are the growth factor and MMP diffusion coefficients, respectively; x,y represents the x and y dimensions of the two-dimensional domain; t is time; m0, αB, γ, χRS, ρRS, βS, μE, ηL, φ, ζ, κ, νB, τI and σ are positive constants with biologically significant values as based on Basanta28. Then, a discretized form of these equations is solved numerically on a two-dimensional lattice that represents a small slice of prostate tissue. All cell types (tumour, basal, luminal, stromal and inflammatory) are modulated by the microenvironment on this lattice and can migrate, proliferate, die and mutate according to the life-cycle flow charts. 

Tumour cells in the model have two continuously variable phenotypes, namely growth factor (γ) and MMP (ζ) production. These traits are passed from a parent tumour cell to its two daughter cells with some small variation, chosen at random from an interval equally weighted in both directions to avoid biased drift. The model is agnostic with respect to specific biological mechanisms that underlie this drift, which could include gradual accumulation of mutations, regulation of gene transcription by epigenetics or aneuploidy, or changes in the number or structure of organelles, for example. The evolution of these phenotypes in time and space is an important consideration of this work.

The switch between stroma (S) and reactive stroma (RS) phenotypes is driven by the growth factor stimulus. Reactive stromal cells are activated if the level of growth factor (G) is above the threshold GRS and are deactivated if the growth factor level G falls below the threshold GRS.
