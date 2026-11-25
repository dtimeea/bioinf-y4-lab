from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

from sklearn.cluster import KMeans
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import classification_report, confusion_matrix
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import LabelEncoder

HANDLE = "dtimeea"

DATA_CSV = Path(
    f"data/work/{HANDLE}/lab08/expression_matrix_{HANDLE}.csv"
)

TEST_SIZE = 0.2
RANDOM_STATE = 42
N_ESTIMATORS = 200

OUT_DIR = Path(f"labs/08_ML_flower/submissions/{HANDLE}")
OUT_DIR.mkdir(parents=True, exist_ok=True)

OUT_CONFUSION = OUT_DIR / f"confusion_rf_{HANDLE}.png"
OUT_REPORT = OUT_DIR / f"classification_report_{HANDLE}.txt"
OUT_FEATIMP = OUT_DIR / f"feature_importance_{HANDLE}.csv"
OUT_CLUSTER_CROSSTAB = OUT_DIR / f"cluster_crosstab_{HANDLE}.csv"


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


def train_random_forest(X_train, y_train,
                        n_estimators, random_state):

    model = RandomForestClassifier(
        n_estimators=n_estimators,
        random_state=random_state,
        n_jobs=-1
    )

    model.fit(X_train, y_train)

    return model


def evaluate_model(model, X_test, y_test,
                   label_encoder, out_png, out_txt):

    y_pred = model.predict(X_test)
    classes = label_encoder.classes_
    labels = np.arange(len(classes))

    report = classification_report(
        y_test,
        y_pred,
        labels=labels,
        target_names=classes,
        zero_division=0
    )

    print("\n=== Random Forest ===")
    print(report)

    out_txt.write_text(report, encoding="utf-8")

    cm = confusion_matrix(
        y_test,
        y_pred,
        labels=labels
    )

    plt.figure(figsize=(7, 5))

    sns.heatmap(
        cm,
        annot=True,
        fmt="d",
        cmap="Blues",
        xticklabels=classes,
        yticklabels=classes
    )

    plt.xlabel("Predicted label")
    plt.ylabel("Actual label")
    plt.title("Random Forest - ER status")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()


def compute_feature_importance(model, feature_names, out_csv):

    importance = pd.DataFrame({
        "Feature": feature_names,
        "Importance": model.feature_importances_
    })

    importance = importance.sort_values(
        "Importance",
        ascending=False
    )

    importance.to_csv(out_csv, index=False)

    print("\nTop 10 gene:")
    print(importance.head(10).to_string(index=False))


def run_kmeans_and_crosstab(X, y, label_encoder,
                            n_clusters, out_csv):

    kmeans = KMeans(
        n_clusters=n_clusters,
        random_state=RANDOM_STATE,
        n_init="auto"
    )

    clusters = kmeans.fit_predict(X)

    df = pd.DataFrame({
        "Label": label_encoder.inverse_transform(y),
        "Cluster": clusters
    })

    ctab = pd.crosstab(df["Label"], df["Cluster"])

    print("\n=== KMeans crosstab ===")
    print(ctab)

    ctab.to_csv(out_csv)


def main():

    X, y = load_dataset(DATA_CSV)
    y_enc, le = encode_labels(y)

    print("Dataset:", X.shape)
    print("Clase:", list(le.classes_))

    X_train, X_test, y_train, y_test = train_test_split(
        X,
        y_enc,
        test_size=TEST_SIZE,
        random_state=RANDOM_STATE,
        stratify=y_enc
    )

    print("Train:", X_train.shape)
    print("Test:", X_test.shape)

    model = train_random_forest(
        X_train,
        y_train,
        N_ESTIMATORS,
        RANDOM_STATE
    )

    evaluate_model(
        model,
        X_test,
        y_test,
        le,
        OUT_CONFUSION,
        OUT_REPORT
    )

    compute_feature_importance(
        model,
        X.columns,
        OUT_FEATIMP
    )

    run_kmeans_and_crosstab(
        X,
        y_enc,
        le,
        len(le.classes_),
        OUT_CLUSTER_CROSSTAB
    )

    print("\nRezultate salvate in:", OUT_DIR)


if __name__ == "__main__":
    main()
