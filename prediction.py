import argparse
import joblib
import pandas as pd
from sklearn.metrics import mean_squared_error, r2_score

def prediction(x, y, saved_model):
    model = joblib.load(saved_model)
    X = pd.read_csv(x)
    y_pred = model.predict(X)
    print(y_pred)
    y_pred = pd.DataFrame(y_pred)
    y_pred.to_csv(y, header=False, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--X", help="Enter the file path to load features file")
    parser.add_argument("--y", help="Enter the file path to save the predictions")
    parser.add_argument("--saved_model", help="Enter the file path to load saved_model")
    args = parser.parse_args()

    prediction(args.X, args.y, args.saved_model)