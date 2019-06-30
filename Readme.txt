Introduction
=============

InPhenotype is a computational approach for combining gene-expression
data together with genotyping data to predict quantitative physiological
traits. The major advantage of InPhenotype is its ability to exploit the 
continuous nature of molecular data while modeling gene-locus interactions

Usage
=====

InPhenotype incorporates a prediction mode, which yields quantitative trait
predictions for unseen individuals. For that purpose InPhenotype takes as 
input a training-set file and a testing-set file.


In the /data/ directory we have a toy dataset which was generated from a
response of 100 features: 50 genes and 50 SNPs across 50 individuals.

Train.txt
This file specifies the training data for all individuals.
The format of this file:
<row> = <phenotype/SNP/Gene name> <individuals 1 value> <individual 2 value>...
The first row must define the phenotype, followed by SNP defining rows and gene
defining rows (in that order).

Test.txt
This file specifies the testing data for new individuals.
The format of this file is similar to the train.txt file, with one difference:
The first phenotype row is optional and depends on whether the user wants to 
obtain an error measure in addition to phenotype predictions. In case the user
doesn't want an error measure, the file should start with a SNP defining row.

It is important to note that train and test files must contain the same genes
and SNPs (amount and types).

To run InPhenotype we enter:
python InPhenotype.py "Train file path" "Test file path" "Include error" "#Genes" "#SNPs" "Individual threshold" "Improvement threshold" "#Trees" "Output path" "#Representing SNPs" "#Final forest genes"

where the arguments are:
Train file path - 		full path for training file
Test file path - 		full path for testing file
Include error - 		a binary value which determines whether to add an error
						measure to predictions - 0/1 values (1 to add)
#Genes - 				number of genes in input files
#SNPs - 				number of SNPs in input files
Individual threshold - 	a positive natural number, smaller than half of the
						samples amount, which determines the individuals 
						threshold described in the manuscript
Improvement threshold - a number between 0 and 1, which determines the percentage
						of minimal improvement between tree branches as described
						in the manuscript. For example - 0.05 demands a minimal
						improvement of 5%. 
#Trees - 				number of total trees created by the InPhenotype model
Output path - 			path for output files
#Representing SNPs - 	initial number of SNPs used for building a single tree.
						Important in case the number of SNPs is very large, which
						might lengthen the running time. In other cases should be
						similar to total number of SNPs
Use leaves weight		a binary value which determines whether to consider leaf scores
						for prediction as a weighting scheme - 0/1 values (1 to include)
#Final forest genes - 	an optional parameter for number of genes to include in
						the final model. As mentioned in the manuscript, genes
						are sorted according to their scores. This parameter
						determines which gene-based trees to include in the final
						prediction model. Default is 1.
					  
Running example:

python InPhenotype.py "Data/Train.txt" "Data/Test.txt" "1" "50" "50" "7" "0.05" "1000" "Data/" "50" "1" "3"

Note that InPhenotype.py here is the full path for the InPhenotype script file. 

Requirements
===========

InPhenotype runs on python 2.7.6 version and requires the following modules:
scipy v1.1.0
numpy v1.15.2



If you have other questions on how to use InPhenotype, please feel free to send 
me an email to: harel@mail.tau.ac.il 
