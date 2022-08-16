import argparse
import pandas as pd
import joblib
from sklearn import linear_model
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score

    

def train(path_features, path_solubilities, test_size, random_state, path_to_trained_model):
    '''
    train the model using X_train and y_train. return a model 
    '''
    features = pd.read_csv(path_features)
    solubilities = pd.read_csv(path_solubilities)
    test_size = float(test_size)
    random_state = int(random_state)

    X_train, X_test, y_train, y_test = train_test_split(features, solubilities, test_size = test_size, random_state = random_state)
    
    model = linear_model.LinearRegression()
    model.fit(X_train, y_train)

    y_pred_train = model.predict(X_train)

    print("Coefficient:", model.coef_)
    print("Intercept:", model.intercept_)
    print("Mean squared error:", mean_squared_error(y_train, y_pred_train))
    print("r2 value:", r2_score(y_train, y_pred_train))

    joblib.dump(model, path_to_trained_model)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--path_features", help="Enter the file path to load molecular descriptors")
    parser.add_argument("--path_solubilities", help="Enter the file path to load solubilities data")
    parser.add_argument("--test_size", help="Enter the numerical value from 0 to 1.0 to indicate the size of the test sample")
    parser.add_argument("--random_state", help="Enter a numerical value for the random seed")
    parser.add_argument("--path_to_trained_model", help='Enter the file path to save the linear model')
    args = parser.parse_args()

    train(args.path_features, args.path_solubilities, args.test_size, args.random_state, args.path_to_trained_model)
