# ML-Project----Predicting-solubility-
A personal project on using sklearn library and RDkit to build a linear model to predict the solubility of molecules. 

SMILES and solubility data were obtained from John S. Delaney, Journal of Chemical Information and Computer Sciences 2004 44 (3), 1000-1005. DOI: 10.1021/ci034243x

This repository mainly consists of 4 folders and 1 txt file. To run the codes succesfully, you may clone the repository. Enter the following into the terminal to clone. 

```
git clone https://github.com/sherwin97/ML-Project----Predicting-solubility-.git
```

1. data_csv_files folder contains csv files of the orginal data from Delaney. SMILES and solubilities data were put into separated CSV files. Features CSV file contain the molecular descriptions obtained from RDkit and is generated using prepare_data.py found in python_scripts folder. 

2. model folder contains the trained_linear_model. Alternatively, you may save your own trained model into this folder. 

3. notebooks folder contains jupyter notebooks used for testing codes.

4. python_scripts folder contains 3 python scripts. 
    1. prepare_data.py takes in the SMILES data and return a csv file containing the molecular descriptions of molecules. Five main descriptions are obtained: LogP, Molecular Weight, Number of Rotatable Bonds, Total Polar Surface Area and Aromatic Proportion. This python scripts take in two arguments. 
      1. path_smiles: path of csv containing SMILES data
      2. path_features: path to save features data
    You may enter the following into the terminal. 
    ```
    $ python ./python_scripts/prepare_data.py --path_smiles ./data_csv_files/smiles.csv --path_features ./data_csv_files/features.csv
    ```
    2. train_linear_model.py takes in the features file, created with the first python script, and solubility file and splits them to training and testing set. Then, the training data will be used to train the model. A trained model will then be returned and saved using the pickle library. In addition, the r2 score will also be return to show the accuracy of model (r2 = 0.75 for Delaney SMILES data). 
       1. path_features: enter the path of features csv
       2. path_solubilities: enter the path of solubilities csv
       3. test_size: enter a float value to indicate size of test sample. 0.2 = 20% taken as testing data
       4. random_state: enter a numerical value to initialise the random state of splitting the data to ensure reproducibility of model
       5. path_to_trained_model: enter path to save the trained model
    You may enter the following into the terminal. 

    ```
    $ python ./python_scripts/train_linear_model.py --path_features ./data_csv_files/features.csv --path_solubilities ./data_csv_files/solubility.csv --test_size 0.2 --random_state 123 --path_to_trained_model ./model/trained_linear_model.pkl
    ```

    3. prediction.py takes in the features file and the saved model to predict the solubilities of molecules. A csv of the predicted solubility is returned. The script takes in 3 arguments. 
         1. path_features: enter the path of features csv
         2. path_predictions: enter the path to save the predicted solubility values as csv file 
         3. path_to_trained_model: enter the path to load the saved model. joblib library is used hence .pkl extension should be used to save and load model successfully. 
    You may enter the following into the terminal.

    ```
    $ python ./python_scripts/prediction.py --path_features ./data_csv_files/features.csv --path_predictions ./data_csv_files/predictions.csv --path_to_trained_model ./model/trained_linear_model.pkl
    ```

5. requirements.txt contains the version informations for libraries used for the codes to function. 