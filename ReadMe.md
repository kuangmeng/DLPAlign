# DLPAlign

A Deep Learning based Progressive Alignment for Multiple Protein Sequences

Sequential modeling-based two-stage Progressive Multiple Sequence Alignment

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


## Models

The old models could be downloaded from `https://drive.google.com/file/d/16xPQSXi6kk1p31xdlF4P3_PToeVcqA2r/view?usp=sharing`. (Please put the unzipped dir `Classifier`in the project home dir.) and the best model is `pairs_classification_cnn_bilstm_final4`. 

The new best model can be downloaded from `https://connecthkuhk-my.sharepoint.com/:u:/g/personal/mmkuang_connect_hku_hk/EZZ-7tIaeXVDirnslWN3RLwBkxVnfTnlSRWZGmysMTLzMA?e=hQzQ2p`.

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
