# Public Release: Postselection technique for optical Quantum Key Distribution with improved de Finetti reductions

This is a public release of the code used in *Postselection technique for optical Quantum Key Distribution with improved de Finetti reductions* \[[arXiv](https://arxiv.org/abs/2403.11851)]. This was built for [v2.0.1](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/releases/tag/v2.0.1) of the Open QKD Security package.

This is the code (and data) required to generate the plot shown in Fig. 2 of the paper, demonstrating the performance of the postselection technique for the three-state protocol.


## Install instructions
> [!CAUTION]
> This repository is for archival and transparency purposes, we do not guarantee compatibility with other versions of the Open QKD Security package beyond the ones listed above.

### As zip
1. Download the linked version of the code from above, and follow all [installation instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/bb1c6490c6bffb0661cef52f6b48de41b5e78027).
2. Also follow the additional Mosek install instructions if you want an exact match.
3. Download the latest release of this code on the side bar, unzip in your preferred directory, and add this folder to the Matlab Path.
4. Run `mainThreeStateFinite.m` to regenerate plots. The original plot and data is available in the folder "ResultsForPaper".

### With git
1. To clone this repository and its exact submodules, navigate to your desired directory and run,
```
git clone --recurse-submodules https://github.com/Optical-Quantum-Communication-Theory/PostselectionTechniqueForOpticalQuantumKeyDistributionWithImprovedDeFinettiReductions
```
2. Follow all further [install instructions](https://github.com/Optical-Quantum-Communication-Theory/openQKDsecurity/tree/bb1c6490c6bffb0661cef52f6b48de41b5e78027).
3. Also follow the additional Mosek install instructions if you want an exact match.
4. Add this repository's folder to the Matlab path.
5. Run `mainThreeStateFinite.m` to regenerate plots. The original plot and data is available in the folder "ResultsForPaper".
