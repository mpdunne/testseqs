![alt text](http://empede.co.uk/imgrepos/testseqs_head.png? "TestSeqs - randomly generated sample sequences")

What does TestSeqs do?
======================

TestSeqs generates simple random sets of amino acid sequences, providing all-purpose simulated gene families.

https://github.com/mpdunne/testseqs

Usage
=====

TestSeqs is written in python, and requires no installation (though you may wish to put the download directory in your bash $PATH). To run TestSeqs with no options whatsoever, simply type ```testseqs```. This will give 10 random amino acid sequences, each derived from an initial sequence of length 40.

Sequence options
================

Each sequence outputted by TestSeqs is derived from an initial sequence, and subsequently mutated according to user preference. Mutations occur at a rate given by the ```-e``` or ```--errorprob``` option, the default value for which is ```0.1```. If a mutation occurs, it is a gap with probability given by ```-g``` or ```--gapprob```, the default value for which is ```0.2```. Increasing either of these values will decrease the similarity between the outputted sequences.

You can change the length of the randomised initial sequence using the ```-l``` or ```--length``` options (default is ```40```), and you can change the number of outputted sequences using ```-n``` or ```--numseqs``` (default is ```10```). Alternatively you can specify your own initial sequence using ```-i``` or ```--initial```.

If you have a specific set of amino acids you wish the sequences to be drawn from, you may use the ```-a``` or ```--aas``` option. For example, ```testseqs --aas "AGHMLPQWR"```.

By default, TestSeqs starts every sequence with an ```M``` and ends each sequence with a ```*```. To turn this off, use the options ```-m false``` and ```-x false``` respectively. To specify an alternative start amino acid, use ```-c``` or ```--altstart```, followed by your amino acid of choice.

Amino acid distributions
========================
By default, TestSeqs chooses a substitute amino acid uniformly at random. Weights can be applied to the amino acid choices in one of two ways. Firstly, a weights file can be specified using ```w``` or ```--weights```, followed by the file name. A weights file is a tab-separated CSV file, the first column of which should contain amino acids, and the second column a weight. It doesn't matter what scale the weights are. For example, a weights file could start:

```
A    3
C    4
D    1
E    29
...
```

Alternatively, weights can be learned from an existing set of amino acid sequences, in FASTA format, specified using the ```-f``` or ```--fasta``` options.

Tree mode
=========
A species tree can be used to generate more realistic sequence relations. You can specify a tree in Newick format using the ```-t``` or ```--tree``` options. The tree must have branch lengths. The ```errorprob``` value means something slightly different here: mutations are modelled as a Poisson process along the branch lengths, with parameter λ=e/(1-e), where e is the specified error probability. 

Gene duplication and loss events can also be included in tree mode. The gene loss parameter is specified using ```-z``` or ```--lossprob```; the gene duplication parameter is specified using ```-y``` or ```--dupprob```. Again these options specify parameters in a Poisson model of gene loss and gene duplication.
