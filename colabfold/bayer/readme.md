Usage of the AlphaPulldown utils:

Note, this workflow is still computing many features more than once. The AlphaPulldown does a way better job at this,
for small scale experiments this is fine and still manageable to compute.

We start with a bait.fasta and a candidates.fasta.

```
>bait
ELVIS
head candidates.fasta
>c1
LIVES
>c2
ISKING
```


1) generate the paired sequences (one vs. many)

```
bash
python ../colabfold/bayer/colab_prepare_pulldown.py --bait bait.fasta --candidates candidates.fasta --output 01_pulldown_experiments.fasta
# result
>bait_and_c1
ELVIS:LIVES
>bait_and_c2
ELVIS:ISKING
```

2) generate the MSA
```
colabfold_search 01_pulldown_experiments.fasta XXXXXX/mmseqs_dbs/ 02_pulldown_alignments --thread 125 --db1 uniref30_2103_db

# result
candidate_and_c1.a3m
candidate_and_c2.a3m
```

3) split the MSA again by number of gpus

```
python ../colabfold/bayer/colab_split_dir.py --input_dir 02_pulldown_alignments --output_dir 03_pulldown_alignments_splitted --n_splits 2
```

4) structure predition

```
...
```

5) summary

...
