import argparse
import pandas as pd
import joblib
from sklearn import linear_model
from sklearn.metrics import mean_squared_error, r2_score

from linear_model import X_test


def train(X_train, y_train, saved_model):
    '''
    train the model using X_train and y_train. return a model 
    '''
    X_train = pd.read_csv(X_train)
    y_train = pd.read_csv(y_train)
    model = linear_model.LinearRegression()
    model.fit(X_train, y_train)

    y_pred_train = model.predict(X_train)

    print("Coefficient:", model.coef_)
    print("Intercept:", model.intercept_)
    print("Mean squared error:", mean_squared_error(y_train, y_pred_train))
    print("r2 value:", r2_score(y_train, y_pred_train))

    joblib.dump(model, saved_model)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--X_train", help="Enter the file path to load X_train")
    parser.add_argument("--y_train", help="Enter the file path to load y_train")
    parser.add_argument("--Model", help="Enter the file path to save model")
    args = parser.parse_args()
    train(args.X_train, args.y_train, args.Model)