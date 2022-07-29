# ML-Project----Predicting-solubility-
A personal project to predict solubility of given molecules using their molecular descriptors. 

Comprises of three python scripts. Scripts can be ran in command line using the follow instructions. 
1. prepare_data.py allow users to obtain the molecules desciptors and prepare files to train the model. Arguments required includes paths to features file, target file, training files and test files. 

python prepare_data.py --X_input "path" --y_input "path" --X_train "path" --X_test "path", --y_train "path", --y_test "path"

2. train_linear_model.py allow users to construct a linear model to predict the solubility of a molecule. Arguments required includes paths to load training files and to save model.

python train_linear_model.py  --X_train "path" , --y_train "path", --Model "path"

3. prediction.py allow users to predict given the data. Arguments required inlcude path to load features file, path to save predicted values and path to load the saved model. 

python prepare_data.py --X "path" --y "path" --saved_model "path" 


Data obtained from:
John S. Delaney
Journal of Chemical Information and Computer Sciences 2004 44 (3), 1000-1005
DOI: 10.1021/ci034243x

Codes adapted from:
Data Professor, https://www.youtube.com/watch?v=VXFFHHoE1wk&list=PLtqF5YXg7GLlQJUv9XJ3RWdd5VYGwBHrP&index=10
