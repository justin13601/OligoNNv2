# OligoNN

A program to train a neural network using gene expression datasets GSE6269 and GSE63990 and use the resulting parameters to generate DNA sequences for a molecular computer modelled after the neural network. Contains files to clean the datasets, train the neural network using PyTorch, store the parameters of a trained model, generate DNA sequences for the molecular model using the trained parameters and store them in a dataframe.

***University of Toronto Institute of Biomedical Engineering: ChouLab - Chou, L.***

*Credits:*
-	Ryan C. Lee
-	Ariel Corsano
-	Justin Xu

## Dependencies
- NumPy
- Pandas
- Biopython
- Matplotlib
- PyTorch
- Torchsummary
- Torchvision
- Json
- Random
- Csv
- Import_ipynb
- Scikit-learn
- NUPACK (not pip installable)

## Package File Tree

```
|--------core
|       |------__init__.py
|       |------check_hamming.py
|       |------energy_matrix.py
|       |------nupack_wrap.py
|       |------oligo_database.py
|       |------oligo_gen.py
|       |------oligo_gen_string.py
|       |------probe_seqs.py
|       |------quant_energy.py
|       |------weight_to_seq_pipeline.py
|------data
|       |------GSE6269
|       |       |------normalized_6269_GPL96.txt
|       |       |------normalized_6269_GPL570.txt
|       |------GSE63990
|       |       |------normalized_63990.txt
|       |       |------GSE63990_series_matrix.txt
|       |------UCDKhan
|       |       |------khan_train.csv
|       |       |------khan_test.csv
|       |       |------khan.rda
|       |------HG-U133A_2.probe_tab
|------FASTA
|       |------// FASTA Files
|------models
|       |------json
|       |       |------pytorch_model_12-24-12-3.json
|       |       |------pytorch_model_5-4-3.json
|       |       |------pytorch_model_3-3-3-3.json
|       |------matlab
|       |       |------exp_geo.m
|       |       |------mdl_batch.m
|       |------pytorch
|       |       |------model_MulticlassANN_bs16_lr0.001_epoch295_12-24-12-3
|       |       |------model_ANN_bs16_lr0.001_epoch1333_5-4-3
|       |       |------model_ANN_bs16_lr0.001_epoch866_3-3-3-3
|------notebooks
|       |------data_cleaner.ipynb
|       |------gene_selection.ipynb
|       |------PyTorch_GSE63990_GSE6269.ipynb
|       |------other
|       |       |------PyTorch_UCDKhan.ipynb
|------weight_to_seq_outputs
|       |-----// Output Directory with .csv Files
|------__init__.py
|------__main__.py
|------README.md
|------README.docx
|------requirements.txt
|------setup.py
```

## Testing
To run `__main__.py`, execute:

    $python3 -m OligoNN

in project directory containing OligoNN.

## Files
.core/check_hamming.py
-	Hamming distance filter implemented in oligo_database
-	Authored: Corsano, A.
-	Last Edited: August 20, 2020

.core/energy_matrix.py
-	Functions to produce and plot energy matrices, as well as calculate binding energy
-	Authored: Chou, L; Edited: Corsano, A.
-	Last Edited: August 20, 2020

.core/nupack_wrap.py
-	Nupack wrapper to convert Nupack functions to python 
-	Authored: Pierce lab, Caltech; Edited: Corsano, A
-	Last Edited: May 28, 2020

.core/oligo_database.py
-	Inserts oligos into database class and can export as csv
-	Authored: Xu, J; Edited: Corsano, A.
-	Last Edited: August 20, 2020

.core/oligo_gen.py
-	Generates DNA strands stored as Bio.Seq objects in class oligo and contains functions to manipulate strands
-	Authored: Corsano, A. 
-	Last Edited: August 13, 2020

.core/oligo_gen_string.py
-	Generates DNA strands stored as strings in class oligo and contains functions to manipulate strands
-	Authored: Xu, J. 
-	Last Edited: August 18, 2020

.core/probe_seqs.py
-	Stores Affymetrix probe sequences (and optional FASTA sequences) of top genes in a dataframe
-	Authored: Corsano, A.
-	Last Edited: July 22, 2020

.core/quant_energy.py
-	Contains functions to convert binding energies to weights and find sequences with specific binding energies
-	Authored: Corsano, A.
-	Last Edited: August 20, 2020

.core/weight_to_seq_pipeline.py
-	Generates DNA sequences for molecular model 
-	Authored: Corsano, A.
-	Last Edited: August 20, 2020

.notebooks/gene_selection.ipynb
-	Gene input selection via sklearn’s chi-2 scoring feature selection
-	Authored: Lee, R; Edited: Corsano, A., Xu, J.
-	Last Edited: Aug 06, 2020

.notebooks/data_cleaner.ipynb
-	Data series matrix processing
-	Authored: Lee, R; Edited: Corsano, A., Xu, J.
-	Last Edited: Aug 08, 2020

.notebooks/PyTorch_GSE63990_GSE6269.ipynb
-	PyTorch neural network testing and training notebook for GSE63990 & GSE6269 data
-	Authored: Xu, J
-	Last Edited: Aug 19, 2020

.notebooks/other/PyTorch_UCDKhan.ipynb
-	PyTorch neural network testing and training notebook for UCDKhan data
-	Authored: Xu, J
-	Last Edited: Aug 02, 2020

.models/matlab/exp_geo.m
-	Host gene expression classifiers diagnose acute respiratory illness etiology
-	Sci Transl Med 2016; PMID: 26791949
-	Authored: Henao, R.; Edited: Xu, J
-	Last Edited: Aug 19, 2020

.models/matlab/mdl_batch.m
-	Batch preparation MATLAB script
-	Authored: Henao, R.
-	Last Edited: Feb 09, 2016

## Contributors
-	Henao, R.
    -	MATLAB exp_geo.m & mdl_batch.m script
-	Pierce lab, California Institute of Technology
    -	nupack_wrap.py
-	Affymetrix Human Genome U133A 2.0 Array - Support Materials
    -	HG-U133A_2.probe_tab

## References
Lopez, R., Wang, R. & Seelig, G. A molecular multi-gene classifier for disease diagnostics. Nature Chem 10, 746–754 (2018). https://doi.org/10.1038/s41557-018-0056-1

Tsalik, E. L., Henao, R., Nichols, M., Burke, T., … Woods, C. W. (2016). Host gene expression classifiers diagnose acute respiratory illness etiology. Science Translational Medicine, 8(322), 322ra11-322ra11. https://doi.org/10.1126/scitranslmed.aad6873

Khan, J., Wei, J. S., Ringnér, M., Saal, L. H. (2001). Classification and diagnostic prediction of cancers using gene expression profiling and artificial neural networks. Nature Medicine, 7(6), 673–679. https://doi.org/10.1038/89044

FASTA Gene Sequences from GeneCards.
