# program to calculate similarity for machine learning implemetation
import numpy as np
import pandas as pd

data_dict = {'protein1': [1, 2, 3, 4, 5], 'protein2': [6, 5, 4, 3, 2], 'protein3': [6, 5, 4, 3, 2],
             'protein4': [[2, 3, 4, 5, 2]], 'protein5': [6, 5, 4, 3, 2]}


class similarty_measure(object):
    def __init__(self, prot1_key_freq, prot2_key_freq):  # accepts two parameters
        self.x = prot1_key_freq
        self.y = prot2_key_freq

    def corelation_cofficient_similarity(X, Y):
        mean_x = np.mean(X)
        mean_y = np.mean(Y)

        if len(X) == len(Y):
            Cov_x_y = 0
            Var_x_squared = 0
            Var_y_squared = 0
            for i in range(len(X)):
                Cov_x_y += ((X[i] - mean_x) * (Y[i] - mean_y))
                Var_x_squared += (X[i] - mean_x) ** 2
                Var_y_squared += (Y[i] - mean_y) ** 2

            corelation_cofficient = Cov_x_y / (np.sqrt(Var_x_squared * Var_y_squared))
            return corelation_cofficient

    def jaccard_similary(self):
        x_intersetion_y = 0
        x_union_y = 0
        if len(self.x) == len(self.y):
            for i in range(len(self.x)):
                x_intersetion_y += min(self.x[i], self.y[i])
                x_union_y += max(self.x[i], self.y[i])
            Similairy = x_intersetion_y / x_union_y
            return Similairy


def Similarity_matrix(*args):
    key_freq = args[0]
    prot_names = key_freq.keys()
    dict_values = list(key_freq.values())

    for i in range(len(prot_names)):
        for j in range(i + 1, len(prot_names)):
            print(dict_values[i], dict_values[j])
            print(similarty_measure.corelation_cofficient_similarity(dict_values[i], dict_values[j]))


# Similarity_matrix(data_dict)

# similairy=similarty_measure(data_dict)
# print(similairy.corelation_cofficient_similarity())
similarity = similarty_measure([5, 3, 2, 1], [3, 2, 5, 0])
print(similarity.jaccard_similary())
# print(similairy.similarity_matrix(protein1,protein2))






