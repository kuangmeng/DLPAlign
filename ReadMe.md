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

## Citation
```
@inproceedings{kuang2020dlpalign,
  title={DLPAlign: A Deep Learning based Progressive Alignment Method for Multiple Protein Sequences},
  author={Kuang, Mengmeng and Liu, Yong and Gao, Lufei},
  booktitle={CSBio'20: Proceedings of the Eleventh International Conference on Computational Systems-Biology and Bioinformatics},
  pages={83--92},
  year={2020}
}
```
