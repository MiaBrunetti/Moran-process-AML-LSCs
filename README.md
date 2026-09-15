# Moran-process-AML-LSCs
This repository accompanies the article "Mathematical modelling of clonal reduction therapeutic strategies in acute myeloid leukemia". It includes the original cytarabine cell viability data and the necessary code for the mathematical analysis.

Authors: Mia Brunetti<sup>1,2</sup>; Isabella A. Iasenza<sup>3,4</sup>; Adrianne L. Jenner<sup>5</sup>; Noël J-M Raynal<sup>2,6</sup>; Kolja Eppert<sup>4,7</sup>; Morgan Craig<sup>1,2</sup><br>
<sup>1</sup>Département de Mathématiques et de Statistiques, Université de Montréal, Montréal, Canada<br>
<sup>2</sup>Sainte-Justine University Hospital Research Center, Montréal, Canada<br>
<sup>3</sup>Division of Experimental Medicine, Department of Medicine, McGill University, Montréal, Canada<br>
<sup>4</sup>Research Institute of the McGill University Health Centre, Canada<br>
<sup>5</sup>School of Mathematical Sciences, Queensland University of Technology, Brisbane, Australia<br>
<sup>6</sup>Département de Pharmacologie et Physiologie, Université de Montréal, Montréal, Canada<br>
<sup>7</sup>Department of Pediatrics, McGill University, Montréal, Canada<br>

## Repository structure
| Folder | Description |
|:-----|:------------|
| Cell viability | Cell viability data and commands for IC50 curve fitting. |
| Cytarabine viability data | Original viability data of HSCs and LSCs under cytarabine. |
| Figure realization | Commands and functions for creating figures. |
| Moran process | Commands for running the Moran process. |
| PKPD responses | Commands for creating PKPD treatement models and for merging MATLAB structures. |
| Pharmacokinetics fitting | Commands for fitting the drugs' pharmacokinetics models. |
| Toxicity fitting | Cardiac glycoside toxicity data and commands for IC50 curve fitting. |

## Requirements
- MATLAB R2025b or later.

## Workflow
### 1. Pharmacokinetics Fitting
Run "Commands_PK_fitting.m" and save work in the [PKPD responses](./PKPD%20response) folder.

### 2. Cell Viability
Run "Commands_Cell_Viability_fitting.m" and save work in the [PKPD responses](./PKPD%20response) folder.

### 3. Toxicity Fitting
Run "Commands_Toxicity_fitting.m" and save work in the [PKPD responses](./PKPD%20response) folder.

### 4. PKPD Responses
Run "Commands_merge_structure.m" and save work in [PKPD responses](./PKPD%20response) folder. Then, run "Commands_Treatment_Model.m" and save work in [Moran process](./Moran%20process) folder.

### 5. Moran Process
Run "Commands_MoranProcessAML.m" and save work in [Figure realization](./Figure%20realization) folder.

### 6. Figure Realization
Run "Commands_Figures.m" to display the graphs necessary to recreate Figures 2-5 and Supplementary Figure S1.

## Citations
If you use any of the data or this code, please cite the following publication: <br>
Brunetti M, Iasenza IA, Jenner AL, Raynal NJM, Eppert K, Craig M. Mathematical modelling of clonal reduction therapeutic strategies in acute myeloid leukemia. Leukemia Research. 2024;140:107485. [https://doi.org/10.1016/j.leukres.2024.107485.](https://github.com/user-attachments/assets/11db3f3d-7ed7-4619-96e9-bb09682bc8db).

