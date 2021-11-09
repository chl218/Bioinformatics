# Bioinformatics


## Finding Hidden Messages in DNA

In the first half of the course, we investigate DNA replication, and ask the question, where in the genome does DNA replication begin?  We will see that we can answer this question for many bacteria using only some straightforward algorithms to look for hidden messages in the genome.

In the second half of the course, we examine a different biological question, when we ask which DNA patterns play the role of molecular clocks.  The cells in your body manage to maintain a circadian rhythm, but how is this achieved on the level of DNA?  Once again, we will see that by knowing which hidden messages to look for, we can start to understand the amazingly complex language of DNA.  Perhaps surprisingly, we will apply randomized algorithms, which roll dice and flip coins in order to solve problems.

Finally, you will get your hands dirty and apply existing software tools to find recurring biological motifs within genes that are responsible for helping Mycobacterium tuberculosis go "dormant" within a host for many years before causing an active infection.

## Genome Sequencing

You may have heard a lot about genome sequencing and its potential to usher in an era of personalized medicine, but what does it mean to sequence a genome?

Biologists still cannot read the nucleotides of an entire genome as you would read a book from beginning to end. However, they can read short pieces of DNA. In this course, we will see how graph theory can be used to assemble genomes from these short pieces. We will further learn about brute force algorithms and apply them to sequencing mini-proteins called antibiotics.

In the first half of the course, we will see that biologists cannot read the 3 billion nucleotides of a human genome as you would read a book from beginning to end.  However, they can read shorter fragments of DNA. In this course, we will see how graph theory can be used to assemble genomes from these short pieces in what amounts to the largest jigsaw puzzle ever put together.

In the second half of the course, we will discuss antibiotics, a topic of great relevance as antimicrobial-resistant bacteria like MRSA are on the rise.  You know antibiotics as drugs, but on the molecular level they are short mini-proteins that have been engineered by bacteria to kill their enemies.  Determining the sequence of amino acids making up one of these antibiotics is an important research problem, and one that is similar to that of sequencing a genome by assembling tiny fragments of DNA.  We will see how brute force algorithms that try every possible solution are able to identify naturally occurring antibiotics so that they can be synthesized in a lab.

Finally, you will learn how to apply popular bioinformatics software tools to sequence the genome of a deadly Staphylococcus bacterium that has acquired antibiotics resistance.

## Comparing Genes, Proteins, and Genomes

Once we have sequenced genomes in the previous course, we would like to compare them to determine how species have evolved and what makes them different.

In the first half of the course, we will compare two short biological sequences, such as genes (i.e., short sequences of DNA) or proteins.  We will encounter a powerful algorithmic tool called dynamic programming that will help us determine the number of mutations that have separated the two genes/proteins.

In the second half of the course, we will "zoom out" to compare entire genomes, where we see large scale mutations called genome rearrangements, seismic events that have heaved around large blocks of DNA over millions of years of evolution.  Looking at the human and mouse genomes, we will ask ourselves: just as earthquakes are much more likely to occur along fault lines, are there locations in our genome that are "fragile" and more susceptible to be broken as part of genome rearrangements?  We will see how combinatorial algorithms will help us answer this question.

Finally, you will learn how to apply popular bioinformatics software tools to solve problems in sequence alignment, including BLAST.

## Molecular Evolution

In the previous course in the Specialization, we learned how to compare genes, proteins, and genomes.  One way we can use these methods is in order to construct a "Tree of Life" showing how a large collection of related organisms have evolved over time.

In the first half of the course, we will discuss approaches for evolutionary tree construction that have been the subject of some of the most cited scientific papers of all time, and show how they can resolve quandaries from finding the origin of a deadly virus to locating the birthplace of modern humans.

In the second half of the course, we will shift gears and examine the old claim that birds evolved from dinosaurs.  How can we prove this?  In particular, we will examine a result that claimed that peptides harvested from a T. rex fossil closely matched peptides found in chickens. In particular, we will use methods from computational proteomics to ask how we could assess whether this result is valid or due to some form of contamination.

Finally, you will learn how to apply popular bioinformatics software tools to reconstruct an evolutionary tree of ebolaviruses and identify the source of the recent Ebola epidemic that caused global headlines.

## Genomic Data Science and Clustering

How do we infer which genes orchestrate various processes in the cell?  How did humans migrate out of Africa and spread around the world? In this class, we will see that these two seemingly different questions can be addressed using similar algorithmic and machine learning techniques arising from the general problem of dividing data points into distinct clusters.

In the first half of the course, we will introduce algorithms for clustering a group of objects into a collection of clusters based on their similarity, a classic problem in data science, and see how these algorithms can be applied to gene expression data.

In the second half of the course, we will introduce another classic tool in data science called principal components analysis that can be used to preprocess multidimensional data before clustering in an effort to greatly reduce the number dimensions without losing much of the "signal" in the data.

Finally, you will learn how to apply popular bioinformatics software tools to solve a real problem in clustering.

## Finding Mutations in DNA and Proteins

In previous courses in the Specialization, we have discussed how to sequence and compare genomes. This course will cover advanced topics in finding mutations lurking within DNA and proteins.

In the first half of the course, we would like to ask how an individual's genome differs from the "reference genome" of the species. Our goal is to take small fragments of DNA from the individual and "map" them to the reference genome.  We will see that the combinatorial pattern matching algorithms solving this problem are elegant and extremely efficient, requiring a surprisingly small amount of runtime and memory.

In the second half of the course, we will learn how to identify the function of a protein even if it has been bombarded by so many mutations compared to similar proteins with known functions that it has become barely recognizable.  This is the case, for example, in HIV studies, since the virus often mutates so quickly that researchers can struggle to study it.  The approach we will use is based on a powerful machine learning tool called a hidden Markov model.

Finally, you will learn how to apply popular bioinformatics software tools applying hidden Markov models to compare a protein against a related family of proteins.

## Big Data in Biology

In this course, you will learn how to use the BaseSpace cloud platform developed by Illumina (our industry partner) to apply several standard bioinformatics software approaches to real biological data.

In particular, in a series of Application Challenges will see how genome assembly can be used to track the source of a food poisoning outbreak, how RNA-Sequencing can help us analyze gene expression data on the tissue level, and compare the pros and cons of whole genome vs. whole exome sequencing for finding potentially harmful mutations in a human sample.

Plus, hacker track students will have the option to build their own genome assembler and apply it to real data!


