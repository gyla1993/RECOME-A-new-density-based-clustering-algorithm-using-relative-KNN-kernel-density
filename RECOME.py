#!/usr/bin/python
# -*- coding=utf-8 -*-

import os
import numpy as np
import math

EC = 2.7182818284


def dist2(a, b, dim):
    distance = 0
    for i in range(dim):
        distance += (a[i] - b[i])*(a[i] - b[i])
    return distance


def calKNN(n, mxK, dim, features):
    td = np.zeros(n, dtype=np.float)
    index = np.zeros(n, dtype=np.int)
    ind_KNN = np.zeros((n, mxK+1), dtype=np.int)
    for i in range(n):
        print(i)
        cnt = 0
        for j in range(n):
            td[j] = dist2(features[i], features[j], dim)
            if i == j:
                td[j] = float('inf')
                continue
            index[cnt] = j
            cnt += 1
        min_index = np.where(td == np.min(td))[0][0]
        for j in range(mxK):
            ind_KNN[i][j+1] = min_index
            td[min_index] = float('inf')
            min_index = np.where(td == np.min(td))[0][0]
    return ind_KNN


def RECOME(n, dim, K, alpha, features, ind_KNN):
    sigma = 0
    rho = np.zeros(n, dtype=np.float)
    for i in range(n):
        tn = 0
        sigma += math.sqrt(dist2(features[i], features[ind_KNN[i][K]], dim))
    sigma /= n
    # calculate NKD
    for i in range(n):
        rho[i] = 0
        for j in range(1, K):
            rho[i] += math.pow(EC, -math.sqrt(dist2(features[i], features[ind_KNN[i][j]], dim)) / sigma)
    # calculate RNKD
    rhoStar = np.zeros(n, dtype=np.float)
    for i in range(n):
        mxr = rho[i]
        for j in range(1, K):
            if rho[ind_KNN[i][j]] > mxr:
                mxr = rho[ind_KNN[i][j]]
        rhoStar[i] = rho[i] / mxr
    # find Higher Density Nearest-neighbor (HDN)
    pre = np.zeros(n, dtype=np.int)
    connect_matrix = np.zeros((n, n), dtype=np.int)
    for i in range(n):
        if rhoStar[i] == 1:
            pre[i] = i
            connect_matrix[i, pre[i]] = 1
            continue
        md = 1e+30
        for j in range(1, K):
            tj = ind_KNN[i][j]
            if rho[tj] > rho[i] and md > dist2(features[i], features[tj], dim):
                md = dist2(features[i], features[tj], dim)
                pre[i] = tj
        connect_matrix[i, pre[i]] = 1
    return connect_matrix


def start_recursion(point_list, start_point, connect_matrix, predict_results, sample_index, atom_cluster_list):
    point_list.append(start_point)
    sample_index[start_point] = 1
    connect_content = connect_matrix[start_point, :]
    target_point = np.where(connect_content == 1)[0][0]
    if target_point == start_point:
        predict_results[point_list] = atom_cluster_list.index(target_point)
        return True
    else:
        start_recursion(point_list, target_point, connect_matrix, predict_results, sample_index, atom_cluster_list)


def generate_clusters(n, connect_matrix):
    atom_center_list = []
    for i in range(n):
        if connect_matrix[i, i] == 1:
            atom_center_list.append(i)
    predict_results = np.zeros(n, dtype=np.int)
    sample_index = np.zeros(n, dtype=np.int)
    for i in range(n):
        if sample_index[i] == 0:
            point_list = []
            start_recursion(point_list, i, connect_matrix, predict_results, sample_index, atom_center_list)
        else:
            continue
    return predict_results


if __name__ == '__main__':
    test_file = "data/in31.txt"
    file = open(test_file, "r")
    content = file.readlines()
    first_line = content[0].replace("\n", "")
    dim = int(first_line.split(" ")[0])
    clusters = int(first_line.split(" ")[1])
    n = len(content) - 1
    features = np.zeros((n, dim), dtype=np.float)
    labels = np.zeros(n, dtype=np.int)
    for i in range(n):
        line = content[i + 1].replace("\n", "")
        line = line.split(" ")
        features[i, 0] = float(line[0])
        features[i, 1] = float(line[1])
        labels[i] = int(line[2])

    mxK = int(math.sqrt(n)) + 1
    ind_KNN = calKNN(n, mxK, dim, features)
    #
    K = int(math.sqrt(n))
    alpha = 0.95
    #
    connect_matrix = RECOME(n, dim, K, alpha, features, ind_KNN)
    # connect_matrix = np.load("./connect_matrix.npy")
    cluster_results = generate_clusters(n, connect_matrix)
    print("Over.")
