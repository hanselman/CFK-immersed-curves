from BifilteredComplex import *

def curves_from_CFKs(filename):
    '''reads in a file immersed_curves/filename.txt which is an output from (the modified
    version of) Szabo's program computing UV=0 knot Floer complexes. This file contains some number of complexes
    associated to knots. This function computes a collection of immersed curves from each complex. The output
    is a new text file input_complex/filename containing the immersed multicurves'''
    f_output = open('immersed_curves/'+filename+'.txt', 'w')
    f_input = open('input_complex/'+filename+'.txt')

    line = f_input.readline()
    reading_complex = False
    knot_counter = 0
    is_alternating = False
    while line:
        if line[0:4] == 'Knot': #this indicates we have reached the start of a new complex in the file
            reading_complex = True
            counter = 0
        if reading_complex:
            if counter == 1 and line[0:4] == 'Rank':   #If K is alternating the output of Szabo's program has this form
                is_alternating = True
                rank_by_grading = {}
                tau = None
            if is_alternating:
                if line[0:5] == 'Total':
                    reading_complex = False
                    delta_grading = tau
    
                    max_grading = max(rank_by_grading.keys())
                    if not max(rank_by_grading.keys()) == -min(rank_by_grading.keys()):
                        print 'error, gradings not symmetric for knot ', knot_counter
                        print rank_by_grading
                    generators_by_grading = [rank_by_grading[max_grading - i] for i in range(2 * max_grading + 1)]
                    # the ith entry in the list generators_by_grading is the number of generators with grading max_grading-i
                    # if knot_counter == 5493:
                    #     print generators_by_grading
    
                    # print 'thin knot 11n' + str(knot_counter), 'tau = ', tau
    
                    if tau <= 0:
                        non_trivial_component = range(tau, -tau + 1)
                    else:
                        non_trivial_component = range(tau, -tau - 1, -1)
    
                    condensed_curve = str(non_trivial_component)
                    # condensed_curve += ', ' + str(delta_grading)
    
                    for i in non_trivial_component:
                        generators_by_grading[max_grading - i] -= 1
    
                    trivial_components = []
                    while sum(generators_by_grading) > 0:
                        i = 0
                        while generators_by_grading[i] == 0:
                            i += 1
                        if not generators_by_grading[i] >= 1:
                            print 'error, top grading not occupied, knot ', knot_counter
                        if not generators_by_grading[i + 1] >= 2:
                            print 'error, top grading - 1 not occupied twice, knot ', knot_counter
                        if i + 2 >= len(generators_by_grading):
                            print 'error with knot ', knot_counter
                            generators_by_grading = [0 for i in range(2 * max_grading + 1)]
                        else:
                            if not generators_by_grading[i + 2] >= 1:
                                print 'error, top grading - 2 not occupied, knot ', knot_counter
    
                            generators_by_grading[i] -= 1
                            generators_by_grading[i + 1] -= 2
                            generators_by_grading[i + 2] -= 1
    
                            new_trivial_component = (max_grading - i - 1, delta_grading)
                            trivial_components.append(new_trivial_component)
    
                    # condensed_curve += ', ' + str(new_trivial_component)
    
                    trivial_components.sort()
                    trivial_components.reverse()
    
                    new_trivial_components = []
                    if len(trivial_components) > 0:
                        (i, j) = trivial_components[0]
                        new_trivial_components.append((i, j, 1))
                        for (i, j) in trivial_components[1:]:
                            if (new_trivial_components[-1][0], new_trivial_components[-1][1]) == (i, j):
                                new_trivial_components[-1] = (i, j, new_trivial_components[-1][2] + 1)
                            else:
                                new_trivial_components.append((i, j, 1))
    
                    for comp in new_trivial_components:
                        condensed_curve += ', ' + str(comp)
    
                    f_output.write(condensed_curve + '\n')
    
    
                else:
                    rank = ''
                    next_letter = line[0]
                    line = line[1:]
                    while not next_letter == ' ':
                        rank += next_letter
                        next_letter = line[0]
                        line = line[1:]
                    rank = int(rank)
                    bigrading = string_to_tuple(line)
                    (i, j) = (int(bigrading[0]), int(bigrading[1]))
                    if tau == None:
                        tau = i - j
                    rank_by_grading[i] = rank

            else:
                if counter == 5:
                    A = string_to_list(line)
                    A = [int(s) for s in A]
                elif counter == 8:
                    M = string_to_list(line)
                    M = [int(s) for s in M]
                elif counter == 11:
                    print knot_counter + 1
                    D = string_to_list(line)
                    D = [string_to_tuple(s) for s in D]
                    D = [(int(a), int(b), int(c)) for (a, b, c) in D]
                elif line[0:3] == 'Tau':
                    tau = int(line[6:len(line)-1])

                    reading_complex = False
                    knot_counter += 1

                    diff = []
                    for i in range(len(A)):
                        diff.append([])

                    for (a, b, c) in D:
                        # a = init_gen, b = final_gen, c = coefficient (should always be 1 for mod 2)
                        if c == 1:
                            Upower = (M[b] + 1 - M[a]) / 2
                            Vpower = A[a] - A[b] + Upower
                            diff[a].append((b, Upower, Vpower))

                    CFK = BifilteredComplex(A, M, diff)

                    if CFK.is_thin():
                        #compute curve components for CFK
                        delta_grading = A[0]-M[0]

                        max_grading = max(A)
                        generators_by_grading = [0 for i in range(2*max_grading + 1)]

                        for item in A:
                            generators_by_grading[max_grading - item] += 1
                        # the ith entry in the list generators_by_grading is the number of generators with grading max_grading-i

                        # print 'thin knot 11n' + str(knot_counter), 'tau = ', tau

                        if tau <= 0:
                            non_trivial_component = range(tau, -tau+1)
                        else:
                            non_trivial_component = range(tau, -tau-1, -1)



                        condensed_curve = str(non_trivial_component)
                        #condensed_curve += ', '+str(delta_grading)


                        for i in non_trivial_component:
                            generators_by_grading[max_grading - i] -= 1

                        trivial_components = [0 for i in range(2*max_grading -1)]
                        #the ith entry of trivial_components is the number of figure 8s whose top is at grading max_grading-i
                        while sum(generators_by_grading) > 0:
                            i = 0
                            while generators_by_grading[i] == 0:
                                i += 1
                            if not generators_by_grading[i] >= 1:
                                print 'error, top grading not occupied'
                            if not generators_by_grading[i+1] >= 2:
                                print 'error, top grading - 1 not occupied twice'
                            if not generators_by_grading[i+2] >= 1:
                                print 'error, top grading - 2 not occupied'

                            generators_by_grading[i] -= 1
                            generators_by_grading[i+1] -= 2
                            generators_by_grading[i+2] -= 1

                            trivial_components[i] += 1

                        for i in range(len(trivial_components)):
                            if trivial_components[i] > 0:
                                condensed_curve += ', ('+ str(max_grading-i-1) + ',' + str(tau) + ')'
                            if trivial_components[i] > 1:
                                condensed_curve += '^' + str(trivial_components[i])

                        f_output.write(condensed_curve + '\n')
                        # print condensed_curve

                    else:
                        train_track = CFK.train_track()
                        print filename + '-' + str(knot_counter), ' # generators = ', len(train_track.alexander_gradings)

                        condensed_curve = train_track.condensed_curve()
                        non_trivial_component = train_track.nontrivial_component

                        f_output.write(condensed_curve + '\n')




                        # print condensed_curve

            counter += 1
        line = f_input.readline()
    f_input.close()
    f_output.close()



def thickness_from_table_file(filename):
    #f_output = open('immersed_curves/'+filename+'.txt', 'w')
    f_input = open('input_complex/'+filename+'.txt')

    line = f_input.readline()
    reading_complex = False
    reading_gradings = False
    knot_counter = 0
    thicknesses = []
    while line:
        if reading_complex:
            if line[0:5] == 'Total':
                reading_gradings = False
                thickness = max(occupied_delta_gradings) - min(occupied_delta_gradings)
                thicknesses.append(thickness)
                if thickness > 2:
                    print 'thickness ', thickness, '  ', '16n'+str(knot_counter)

            if reading_gradings:
                pair = line
                while pair[0] != ' ':
                    pair = pair[1:]
                pair = string_to_tuple(pair)
                pair = (int(pair[0]), int(pair[1]))
                occupied_delta_gradings.append( pair[0]-pair[1] )

            if line[0:4] == 'Seif':
                genus = int(line[15:])
            if line[0:4] == 'Epsi':
                epsilon = int(line[9:])
                reading_complex = False
                knot_counter += 1
                if genus == 2 and epsilon == 0:
                    print 'g=2, e=0', '  ', filename+str(knot_counter)

        if line[0:4] == 'Rank':
            reading_complex = True
            reading_gradings = True
            occupied_delta_gradings = []


        line = f_input.readline()

    print '# knots checked: ', knot_counter
    f_input.close()










curves_from_CFKs('12n')


# curves_from_CFKs('14n-1')
# curves_from_CFKs('14n-2')
#
# curves_from_CFKs('15n-0')
# curves_from_CFKs('15n-1')
# curves_from_CFKs('15n-2')
# curves_from_CFKs('15n-3')
# curves_from_CFKs('15n-4')

# curves_from_CFKs('15n-5')
# curves_from_CFKs('15n-6')
# curves_from_CFKs('15n-7')
# curves_from_CFKs('15n-8')

# curves_from_CFKs('15n-9')
# curves_from_CFKs('15n-10')
# curves_from_CFKs('15n-11')
# curves_from_CFKs('15n-12')
#
# curves_from_CFKs('15n-13')
# curves_from_CFKs('15n-14')
# curves_from_CFKs('15n-15')
# curves_from_CFKs('15n-16')

# curves_from_CFKs('16n-0')
# curves_from_CFKs('16n-1')
# curves_from_CFKs('16n-2')
# curves_from_CFKs('16n-3')
# curves_from_CFKs('16n-4')
# curves_from_CFKs('16n-5')

