This folder contains a python script, MLC2_test_set.py (and its corresponding Jupiter notebook) to take in a .csv file of protein coding genes and predict if a given sample has Tuberculosis meningitis (TBM, class 1) or Other Neurological diseases (OND, class 0). Please ensure that the input .csv file is in the correct format (go to folder 1.format_input). The file used for training uses ENSG IDs for gene IDs. The output of this python script has two columns, with column 1 = sample ID and column 2, final prediction (1 = TBM or 0 = OND). Each sample is predicted 100 times with a random samples set = 80% of training_validation.csv

Demo
Usage: python MLC2_test_set.py input_formatted_test.csv final_predictions.csv

The output of this python script has two columns, with column 1 = sample ID and column 2, final prediction (1 = TBM or 0 = OND). The expected run time for demo (149 samples) on normal desktop computer is 7m 25.572s.


Instructions for use

python MLC2_test_set.py data_files/combined_test.csv data_files/test_preds_MLC2.csv


