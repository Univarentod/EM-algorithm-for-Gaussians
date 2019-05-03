import random
from math import sqrt
from math import exp
from math import pow
from math import log
import copy
from numpy import multiply
from numpy import subtract
from numpy import add
from numpy import dot
from numpy import all
from numpy.linalg import inv
from numpy.linalg import det
from numpy.linalg import eigvals

def from_vector_to_matrix(a, b):
    if len(a) != len(b):
        return
    k = len(a)
    matrix = list()
    for i in range(k):
        matrix.append(list())
        for j in range(k):
            matrix[i].append((a[i]-b[i])*(a[j]-b[j]))
    return matrix

def get_unduplicated_random_number(required_num, range_left, range_right):
    random_numbers = list()
    for i in range(required_num):
        add_fail = True
        while add_fail:
            r = random.randint(range_left, range_right)
            if r not in random_numbers:
                random_numbers.append(r)
                add_fail = False
    return random_numbers

#input_data would be a list of lists
def EM(
input_data,
k = 2
):
    if k < 1:
        return
    if k == 1:
        return input_data
    if type(input_data) != list:
        return
    if len(input_data) < k:
        return
    pow_n_of_2pi = pow(3.14159265*2, len(input_data[0]))                            #Initialization begins
    clusters_information = list()
    clusters = list()
    sum_of_weights = list()
    origin = list()
    for i in range(len(input_data[0])):
        origin.append(0)
    sum_for_new_mean = list()
    sum_for_new_cov = list()
    det_of_matrix = list()
    inv_of_matrix = list()
    random_numbers = get_unduplicated_random_number(k, 0, len(input_data)-1)
    modified_points_flags = list()
    modified_points = list()
    for i in range(k):
        modified_points_flags.append(list())
        for j in range(len(input_data[0])):
            modified_points_flags[i].append(False)
        modified_points.append(list())
        for m in range(len(input_data[0])):
            modified_points[i].append(list())
            for n in range(len(input_data[0])):
                modified_points[i][m].append(0)
        clusters.append(list())
        clusters_information.append(list())
        clusters_information[i].append(1/k)
        clusters_information[i].append(copy.deepcopy(input_data[random_numbers[i]]))
        clusters_information[i].append(from_vector_to_matrix(origin, origin))
        clusters_information[i].append(list())
        for j in range(len(input_data)):
            clusters_information[i][3].append(0)
        sum_of_weights.append(0)
        sum_for_new_mean.append(copy.deepcopy(origin))
        sum_for_new_cov.append(from_vector_to_matrix(origin, origin))
        det_of_matrix.append(0)
        inv_of_matrix.append(from_vector_to_matrix(origin, origin))
    for i in range(k):
        while not all(modified_points_flags[i]):
            random_numbers = get_unduplicated_random_number(1, 0, len(input_data)-1)
            clusters_information[i][2] = from_vector_to_matrix(input_data[random_numbers[0]], clusters_information[i][1])
            eigen_values = eigvals(clusters_information[i][2]).tolist()
            for j in range(len(input_data[0])):
                if not modified_points_flags[i][j]:
                    if eigen_values[j] > 0.001:
                        modified_points[i][j][j] = eigen_values[j]
                        modified_points_flags[i][j] = True
        clusters_information[i][2] = modified_points[i]
    del random_numbers
    del modified_points_flags                                                             #Initialization ends
    change = True
    step = 0
    while change and step < 100000000:
        change = False
        step += 1
        for i in range(k):
            sum_for_new_mean[i] = copy.deepcopy(origin)
            sum_for_new_cov[i] = from_vector_to_matrix(origin, origin)
            det_of_matrix[i] = sqrt(pow_n_of_2pi*det(clusters_information[i][2]))
            inv_of_matrix[i] = inv(clusters_information[i][2]).tolist()
            sum_of_weights[i] = 0
        for j in range(len(input_data)):                                                       #E step
            sum_of_ownership = 0
            first_value = log(clusters_information[0][0]/det_of_matrix[0])
            first_value += -0.5*dot(subtract(input_data[j], clusters_information[0][1]), dot(inv_of_matrix[0], subtract(input_data[j], clusters_information[0][1])))
            for i in range(k):
                clusters_information[i][3][j] = log(clusters_information[i][0]/det_of_matrix[i])
                clusters_information[i][3][j] += -0.5*dot(subtract(input_data[j], clusters_information[i][1]), dot(inv_of_matrix[i], subtract(input_data[j], clusters_information[i][1])))
                clusters_information[i][3][j] -= first_value
                clusters_information[i][3][j] = exp(clusters_information[i][3][j])
                sum_of_ownership += clusters_information[i][3][j]
            for i in range(k):
                clusters_information[i][3][j] /= sum_of_ownership
                sum_of_weights[i] += clusters_information[i][3][j]
                sum_for_new_mean[i] = add(sum_for_new_mean[i], multiply(clusters_information[i][3][j], input_data[j])).tolist()
                sum_for_new_cov[i] = add(sum_for_new_cov[i], multiply(clusters_information[i][3][j], from_vector_to_matrix(input_data[j], clusters_information[i][1]))).tolist()
        for i in range(k):                                                                   #M step
            temp_weight = sum_of_weights[i]/len(input_data)
            if (clusters_information[i][0] - temp_weight) > 0.000001:
                change = True
            clusters_information[i][0] = temp_weight                                        #computing the weights of each distribution
            clusters_information[i][1] = multiply(1/sum_of_weights[i], sum_for_new_mean[i]).tolist()
            clusters_information[i][2] = multiply(1/sum_of_weights[i], sum_for_new_cov[i]).tolist()
    for i in range(len(input_data)):
        max_probability = clusters_information[0][3][i]
        max_index = 0
        for j in range(1, k):
            if clusters_information[j][3][i] > max_probability:
                max_probability = clusters_information[j][3][i]
                max_index = j
        clusters[max_index].append(input_data[i])
    return [clusters, [elements[0:3] for elements in clusters_information]]
