# DLPAlign

A Deep Learning based Progressive Alignment for Multiple Protein Sequences
## Build

```
cd DLPAlign/
make
```

## Training

```
cd Classifier/
python pair_classifier.py
```

or

```
jupyter-notebook pair_classifier.ipynb
```


## Models

The models could be downloaded from `https://drive.google.com/file/d/16xPQSXi6kk1p31xdlF4P3_PToeVcqA2r/view?usp=sharing`. (Please put the unziped dir `Classifier`in the project home dir.)

The best model is `pairs_classification_cnn_bilstm_final4`.

## Run

`python DLPAlign.py <bench_dir>`