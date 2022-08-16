import argparse
import joblib
import pandas as pd

from sklearn.metrics import r2_score


def prediction(path_features, path_prediction, path_to_trained_model):
    #check if necessary columns are present in the dataframe
    col_names = ['MolLogP', 'MolWt','NumRotableBonds', 'TPSA', 'AromaticProportion']
    features_ori = pd.read_csv(path_features)
    #check if features_df consist of necessary informations, else throw error message 
    for name in col_names:
        if name not in features_ori.columns:
            print(f'{name} not present')
        else:
            features_subset = features_ori.loc[:, col_names]
            model = joblib.load(path_to_trained_model)        
            y_pred = model.predict(features_subset)
            y_pred = pd.DataFrame(y_pred)
            y_pred.to_csv(path_prediction, header=False, index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("--path_features", help="Enter the file path to load features file")
    parser.add_argument("--path_predictions", help="Enter the file path to save the predictions")
    parser.add_argument("--path_to_trained_model", help="Enter the file path to load saved_model")
    args = parser.parse_args()

    prediction(args.path_features, args.path_predictions, args.path_to_trained_model)


