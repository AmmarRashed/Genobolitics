# Genobolitics

Genobolitics is a bioinformatics tool that extends [Metabolitics](https://metabolitics.readthedocs.io/en/latest/) to include gene expression data for pathway-level analysis of diseases. This project aims to bridge the gap between genome-scale modeling and pathway-level analysis by integrating metabolic data with gene expression data.

Metabolitics is described in the paper:
> Cakmak A, Celik MH. Personalized Metabolic Analysis of Diseases. IEEE/ACM Trans Comput Biol Bioinform. 2021 May-Jun;18(3):1014-1025. doi: 10.1109/TCBB.2020.3008196. [PubMed](https://pubmed.ncbi.nlm.nih.gov/32750887/)

## Project Overview

Genobolitics uses gene expression data to predict reactions' flux distribution and analyze pathway-level alterations in diseased cells. The approach includes:

- Flux Variability Analysis (FVA) with dynamically built linear programming models
- Translation of boolean relations among genes to reaction coefficients
- Calculation of Diff values at reaction and pathway levels
- Machine Learning-based classification for disease association

The project has been tested with breast cancer and lung cancer datasets, achieving classification accuracy between 82-90%.

## Repository Structure

```
.
├── classes/
│   ├── ParsingDataset.pyc
│   ├── geno_classifier.py       # Contains classification algorithms
│   ├── geno_utils.py            # Utility functions for gene processing
│   └── genobolitics.py          # Core implementation
├── dataset/                     # Contains the datasets used
├── notebooks/                   # Jupyter notebooks with analysis
├── results/                     # Analysis results
├── requirements.txt             # Project dependencies
└── final_presentation.pdf       # Project presentation
```

## Installation

```bash
git clone https://github.com/yourusername/Genobolitics.git
cd Genobolitics
pip install -r requirements.txt
```

## Dependencies

The project depends on several Python libraries including:
- NumPy
- Pandas
- scikit-learn
- COBRApy
- JobLib
- GEOparse
- PyHGNC
- [Metabolitics](https://metabolitics.readthedocs.io/en/latest/)

See `requirements.txt` for the complete list.

## Features

- Integration of gene expression data with metabolic models
- Dynamic Linear Programming for flux analysis
- Pathway-level differential analysis with Diff values
- Machine learning classification for disease prediction
- Hierarchical clustering for disease ontology construction
- Statistical significance analysis with ANOVA

## How to Contribute

Contributions to this project are welcome. Please feel free to submit a pull request or open an issue for bugs, questions, or feature requests.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgements

Thanks to Prof. Ali Çakmak and Muhammed Hasan Çelik for their support throughout the research cycle.
