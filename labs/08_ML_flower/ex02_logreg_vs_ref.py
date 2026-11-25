from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import classification_report, accuracy_score
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder, StandardScaler

HANDLE = "dtimeea"

DATA_CSV = Path(
    f"data/work/{HANDLE}/lab08/expression_matrix_{HANDLE}.csv"
)

TEST_SIZE = 0.2
RANDOM_STATE = 42
N_ESTIMATORS = 200
MAX_ITER_LOGREG = 1000

OUT_DIR = Path(f"labs/08_ML_flower/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_REPORT_TXT = OUT_DIR / f"rf_vs_logreg_report_{HANDLE}.txt"


def ensure_exists(path):
    if not path.is_file():
        raise FileNotFoundError(path)


def load_dataset(path):
    ensure_exists(path)

    df = pd.read_csv(path)

    X = df.iloc[:, :-1].apply(pd.to_numeric)
    y = df.iloc[:, -1]

    return X, y


def encode_labels(y):
    le = LabelEncoder()
    y_enc = le.fit_transform(y)

    return y_enc, le


def train_models(X_train, y_train):

    rf = RandomForestClassifier(
        n_estimators=N_ESTIMATORS,
        random_state=RANDOM_STATE,
        n_jobs=-1
    )

    rf.fit(X_train, y_train)

    scaler = StandardScaler()
    X_train_scaled = scaler.fit_transform(X_train)

    logreg = LogisticRegression(
        max_iter=MAX_ITER_LOGREG,
        random_state=RANDOM_STATE
    )

    logreg.fit(X_train_scaled, y_train)

    return rf, logreg, scaler


def compare_models(rf, logreg, scaler,
                   X_test, y_test, label_encoder, out_txt):

    y_pred_rf = rf.predict(X_test)

    X_test_scaled = scaler.transform(X_test)
    y_pred_lr = logreg.predict(X_test_scaled)

    classes = label_encoder.classes_
    labels = np.arange(len(classes))

    report_rf = classification_report(
        y_test,
        y_pred_rf,
        labels=labels,
        target_names=classes,
        zero_division=0
    )

    report_lr = classification_report(
        y_test,
        y_pred_lr,
        labels=labels,
        target_names=classes,
        zero_division=0
    )

    acc_rf = accuracy_score(y_test, y_pred_rf)
    acc_lr = accuracy_score(y_test, y_pred_lr)

    report = (
        "=== Random Forest ===\n"
        + report_rf
        + "\n=== Logistic Regression ===\n"
        + report_lr
        + "\n=== Model comparison ===\n"
        + f"Random Forest accuracy: {acc_rf:.4f}\n"
        + f"Logistic Regression accuracy: {acc_lr:.4f}\n"
    )

    print(report)

    out_txt.write_text(report, encoding="utf-8")

    print("Raport salvat:", out_txt)


def main():

    X, y = load_dataset(DATA_CSV)
    y_enc, le = encode_labels(y)

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y_enc,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y_enc
    )

    print("Dataset:", X.shape)
    print("Train:", X_train.shape)
    print("Test:", X_test.shape)

    rf, logreg, scaler = train_models(
        X_train,
        y_train
    )

    compare_models(
        rf,
        logreg,
        scaler,
        X_test,
        y_test,
        le,
        OUT_REPORT_TXT
    )


if __name__ == "__main__":
    main()
