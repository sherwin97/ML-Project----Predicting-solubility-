from linear_model import model
from linear_model import X_test, y_test

from sklearn.metrics import mean_squared_error, r2_score

y_pred_test = model.predict(X_test)
print("Coefficient:", model.coef_)
print("Intercept:", model.intercept_)
print("Mean squared error:", mean_squared_error(y_test, y_pred_test))
print("r2 value:", r2_score(y_test, y_pred_test))