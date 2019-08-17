## Contains some simple linear algebra functions for matrices over Z/2Z

def identity_matrix(n):
    '''returns the identity matrix of size n'''
    result = []
    for i in range(n):
        new_row = []
        for j in range(n):
            if j == i:
                new_row.append(1)
            else:
                new_row.append(0)
        result.append(new_row)
    return result

def list_sum_mod2(l1, l2):
    '''adds two lists (of equal length) term by term, modulo 2
    assumes both lists are populated by numbers'''
    return [(l1[i]+l2[i])%2 for i in range(len(l1))]

def print_matrix(M):
    for row in M:
        print row

def matrix_multiply_mod2(M1, M2):
    '''M1 and M2 are n x m and m x k matrices, with coefficients in Z/2Z
    returns the product'''

    if not len(M1[0]) == len(M2):
        print 'error, matrices have wrong dimension for multiplication'
        return None

    n = len(M1)
    m = len(M1[0])
    k = len(M2[0])

    result = []
    for i in range(n):
        new_row = []
        for j in range(k):
            new_entry = 0
            for l in range(m):
                new_entry += M1[i][l]*M2[l][j]
            new_row.append(new_entry%2)
        result.append(new_row)
    return result

def product_of_elementary_matrices_mod2(n, list_of_indices):
    '''returns the product of a list of n x n elementary matrices
    The input is a list of tuples (i,j), each representing the elementary matrix A_i,j,
    which has 1s on the diagonal and in the entry (i,j), and 0s elsewhere

    Input does not allow for transposition elementary matrices, but T_i,j can be replaced
    with A_i,j A_j,i A_i,j'''

    result = identity_matrix(n)
    for (i,j) in reversed(list_of_indices):
        result[i] = list_sum_mod2(result[i], result[j])  # add the jth row to the ith row
    return result


def LPU_decomp_mod2(M):
    '''M is square invertible matrix over Z/2Z
    returns a triple (L, P, U) where M = LPU and L is lower triangular,
    P is a permutation matrix, and U is upper triangular'''

    n = len(M)
    original_matrix = [row[:] for row in M]

    #base case
    if n == 1:
        return ([[1]],[[1]],[[1]])

    if not 1 in [row[0] for row in M]:
        print 'error: matrix not invertible'
    pivot_i = 0
    while not M[pivot_i][0] == 1:
        pivot_i += 1
    #pivot_i is the index i of the first nonzero entry in the 0th column

    L_list = []
    for i2 in range(pivot_i+1,n):
        if M[i2][0] == 1:
            L_list.append(i2)   #This corresponds to the (lower triangular) elementary matrix A_i2,pivot_i
            M[i2] = list_sum_mod2(M[pivot_i],M[i2])    #adding (pivot_i)th row to the (i2) row, corresponding to left multiplication by A_i2,pivot_i

    ## at this point now we have that M = (L_k x ... x L_1 x original matrix), where L_m is the elementary matrix
    ## corresponding to the pair (pivot_i,i2) pair, where i2 in the mth entry of L_list

    U_list = []
    for j in range(1,n):
        if M[pivot_i][j] == 1:
            U_list.append(j)    #This corresponds to the (upper triangular) elementary matrix A_0,j
            for k in range(n):
                M[k][j] = (M[k][j] + M[k][0])%2     #adding the 0th column to the jth column, corresponding to right mult by A_0,j

    ## at this point we have that M = (L_k x ... x L_1 x original matrix x U_1 x ... U_l),
    ## where U_m is the elementary matrix A_0,j, with j the mth entry in U_list
    ## M now has a 1 in the entry (0,pivot_i), and no other 1s in the 0th column or the (pivot_i)th row

    ## We have [original matrix] = (L_1 x ... x L_k x M x U_l x ... x U_1)

    # Compute L = L_1 x ... x L_k:
    L = product_of_elementary_matrices_mod2(n,[(i,pivot_i) for i in L_list])

    # Compute U = U_l x ... x U_1
    U = product_of_elementary_matrices_mod2(n, [(0, j) for j in reversed(U_list)])

    ## now (original matrix) = L x M x U

    # print 'L = '
    # print_matrix(L)
    # print
    #
    # print 'M = '
    # print_matrix(M)
    # print
    #
    # print 'U = '
    # print_matrix(U)
    # print
    #
    # check = matrix_multiply_mod2(L,M)
    # check = matrix_multiply_mod2(check, U)
    # print 'check: ', check == original_matrix

    M_prime = []
    for i in range(n):
        if not i == pivot_i:
            new_row = M[i][:]
            new_row.pop(0)
            M_prime.append(new_row)
    ## M_prime is M with the (pivot_i)th row and 0th column removed

    # print 'M_prime = '
    # print_matrix(M_prime)
    # print

    (L_prime, P_prime, U_prime) = LPU_decomp_mod2(M_prime)

    # print 'L_prime = '
    # print_matrix(M_prime)
    # print
    #
    # print 'P_prime = '
    # print_matrix(M_prime)
    # print
    #
    # print 'U_prime = '
    # print_matrix(M_prime)
    # print
    #
    # check = matrix_multiply_mod2(L_prime,P_prime)
    # check = matrix_multiply_mod2(check, U_prime)
    # print 'check2: ', check == M_prime
    #
    # print 'pivot_i = ', pivot_i

    #insert a new 0th row and 0th column to U_prime, with a 1 as the intersecting entry
    U_prime = [[0]+row for row in U_prime]
    new_row = [0 for item in U_prime[0]]
    new_row[0] = 1
    U_prime = [new_row]+U_prime

    #insert a new (pivot_i)th row and 0th column to P_prime, with a 1 as the intersecting entry
    P_prime = [[0] + row for row in P_prime]
    new_row = [0 for item in P_prime[0]]
    new_row[0] = 1
    P_prime = P_prime[:pivot_i] + [new_row] + P_prime[pivot_i:]

    # insert a new (pivot_i)th row and (pivot_i)th column to L_prime, with a 1 as the intersecting entry
    L_prime = [row[:pivot_i] + [0] + row[pivot_i:] for row in L_prime]
    new_row = [0 for item in L_prime[0]]
    new_row[pivot_i] = 1
    L_prime = L_prime[:pivot_i] + [new_row] + L_prime[pivot_i:]

    ## Now M = L_prime P_prime U_prime
    ## Thus (original matrix) = L x L_prime x P_prime x U_prime x U

    # print 'L_prime = '
    # print_matrix(M_prime)
    # print
    #
    # print 'P_prime = '
    # print_matrix(M_prime)
    # print
    #
    # print 'U_prime = '
    # print_matrix(M_prime)
    # print
    #
    # check = matrix_multiply_mod2(L_prime,P_prime)
    # check = matrix_multiply_mod2(check, U_prime)
    # print 'check3: ', check == M


    L_result = matrix_multiply_mod2(L, L_prime)
    U_result = matrix_multiply_mod2(U_prime, U)

    ### check:
    #if not matrix_multiply_mod2(L_result,matrix_multiply_mod2(M,U_result)) == original_matrix:
    #    print 'error in LPU decomposition'


    # print 20*'<'
    return (L_result, P_prime, U_result)

# M = product_of_elementary_matrices_mod2(9, [(5,2),(5,6),(3,4),(2,7),(2,5),(2,6),(2,7),(8,1),(8,5),(8,6),(8,7),(1,8),(1,0),(0,1),(0,5),(0,6),(0,7),(0,6)])


def decompose_lower_triangular(L):
    '''L is an invertible n x n lower triangular matrix, over Z/2Z
    It can be realized as a product L_1 x ... x L_k, where each L_m is a type A_i,j elementary matrix
    with i > j and with (i - j) decreasing
    In terms of arrows, this corresponds to a sequence of upward arrows from strand i to strand j, with shorter arrows to the right
    returns list of tuples (i,j) representing such a decomposition'''

    result = []

    n = len(L)
    for d in range(1,n):    # d = i-j
        for j in range(0,n-d):
            i = j + d
            if L[i][j] == 1:
                result.append((i,j))
                #modify L by multiplying on the right by A_i,j. That is, add ith column to jth column
                for r in range(n):
                    L[r][j] = (L[r][i] + L[r][j])%2
    if not L == identity_matrix(n):
        print 'error in decompose_lower_triangular'

    return result[::-1]

def decompose_upper_triangular(U):
    '''U is an invertible n x n upper triangular matrix, over Z/2Z
    It can be realized as a product U_1 x ... x U_k, where each U_m is a type A_i,j elementary matrix
    with i < j and with (j - i) increasing
    In terms of arrows, this corresponds to a sequence of downward arrows from strand i to strand j, with shorter arrows to the left
    returns list of tuples (i,j) representing such a decomposition'''

    result = []

    n = len(U)
    for d in range(1,n):    # d = j-i
        for i in range(0,n-d):
            j = i + d
            if U[i][j] == 1:
                result.append((i,j))
                #modify U by multiplying on the left by A_i,j. That is, add jth row to ith row
                U[i] = list_sum_mod2(U[i],U[j])
    if not U == identity_matrix(n):
        print 'error in decompose_upper_triangular'

    return result



