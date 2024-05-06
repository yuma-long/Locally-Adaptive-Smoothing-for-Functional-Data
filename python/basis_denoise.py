import numpy as np
import pandas as pd
from sklearn.model_selection import KFold, train_test_split, GridSearchCV
import tv_denoiser

def main():
    # trainval
    coefs_trainval_df = pd.read_csv("../basis_coefs_trainval.csv", index_col=0, header=0)

    # 5-fold CV で最適なハイパーパラメータを求める
    x = coefs_trainval_df.iloc[0].to_numpy()
    Y = coefs_trainval_df[1:].to_numpy() # 各行が一つの基底関数の係数
    
    lam_candidates = np.logspace(-2, 0, 7)
    kf = KFold(n_splits=5, shuffle=True, random_state=0)
    lam_error_list = []
    for lam in lam_candidates:
        error_list = []
        for train_index, test_index in kf.split(Y.T):
            x_train, x_test = x[train_index], x[test_index]
            Y_train, Y_test = Y[:, train_index], Y[:, test_index]
            
            argsort_train = x_train.argsort()
            x_train_sorted = x_train[argsort_train]
            
            Y_train_sorted = Y_train[:, argsort_train]

            denoised = np.array([tv_denoiser.run_flsa(row, lam) for row in Y_train_sorted])
            
            prediction_list = []
            
            for x_ in x_test:
                idx = np.searchsorted(x_train_sorted, x_)
                x_prediction = []
                for i in range(len(denoised)):
                    if idx == 0:
                        x_prediction.append(denoised[i, 0])
                    elif idx == len(x_train_sorted):
                        x_prediction.append(denoised[i, -1])
                    else:
                        x_prediction.append(((x_train_sorted[idx] - x_) * denoised[i, idx-1] + (x_ - x_train_sorted[idx-1]) * denoised[i, idx])/ (x_train_sorted[idx] - x_train_sorted[idx-1]))
                prediction_list.append(x_prediction)
            
            prediction_matrix = np.array(prediction_list).T
            
            squared_error = np.sum((Y_test - prediction_matrix) ** 2)
            error_list.append(squared_error)
        print(f"Average squared error for lambda = {lam}: {np.mean(error_list)}")
        lam_error_list.append((lam, np.mean(error_list)))
        
    best_lam = min(lam_error_list, key=lambda x: x[1])[0]
    print(f"Best lambda: {best_lam}")
    
    x_argsort = x.argsort()
    x_sorted = x[x_argsort]
    Y_sorted = Y[:, x_argsort]
    denoised = np.array([tv_denoiser.run_flsa(row, best_lam) for row in Y_sorted])

    # test
    coefs_test_df = pd.read_csv("../basis_coefs_test.csv", index_col=0, header=0)
    
    x = coefs_test_df.iloc[0].to_numpy()
    Y = coefs_test_df[1:].to_numpy()
    
    prediction_list = []
    for x_ in x:
        idx = np.searchsorted(x_sorted, x_)
        x_prediction = []
        for i in range(len(denoised)):
            if idx == 0:
                x_prediction.append(denoised[i, 0])
            elif idx == len(x_sorted):
                x_prediction.append(denoised[i, -1])
            else:
                x_prediction.append(((x_sorted[idx] - x_) * denoised[i, idx-1] + (x_ - x_sorted[idx-1]) * denoised[i, idx])/ (x_sorted[idx] - x_sorted[idx-1]))
        prediction_list.append(x_prediction)
    
    prediction_matrix = np.array(prediction_list).T
    prediction_df = pd.DataFrame(prediction_matrix, index=coefs_test_df.index[1:], columns=coefs_test_df.columns)
    prediction_df.to_csv("../basis_coefs_denoised_test.csv")


if __name__ == "__main__":
    main()
