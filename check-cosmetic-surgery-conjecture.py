from BifilteredComplex import *


#checks for cosmetic crossings
def check_conjecture_for_knots(filename):
    '''finds all potential counterexamples to the Cosmetic surgery conjecture on a given set of knots,
    where input_complex/filename.txt is a file containing the CFK complex for these knots '''

    f_output = open('possible_cosmetic_surgeries/possible'+filename+'.txt', 'w')
    f_output_short = open('possible_cosmetic_surgeries/possible' + filename + 'short.txt', 'w')
    f_input = open('input_complex/'+filename+'.txt')

    line = f_input.readline()
    reading_complex = False
    knot_counter = 0        # keeps track of the total number of knots
    possible_counter = 0    # keeps track of the # of knots with epsilon = 0 and genus = 2
    definite_counter = 0    # keeps track of the # of knots for which CFK can not rule out truly cosmetic surgeries
    epsilon_counter = 0     # keeps track of the # of knots with epsilon = 0
    tau_counter = 0         # keeps track of the # of knots with tau = 0
    max_thickness = 0       # keeps track of the maximum Heegaard Floer thickness of all knots checked
    while line:
        if line[0:4] == 'Knot':     #Signals that we have started a new complex in the file
            reading_complex = True
            counter = 0
            is_alternating = False
        if reading_complex:
            #print counter, is_alternating, line
            if counter == 1 and line[0:4] == 'Rank':   #If K is alternating the output of Szabo's program has this form
                is_alternating = True
            if is_alternating:
                # In this case it is enough to read off the rank in each grading. We only need to record the
                # first five gradings, since we will throw out any alternating knot with genus not equal to 2
                if counter == 2:
                    r1 = line
                elif counter == 3:
                    r2 = line
                elif counter == 4:
                    r3 = line
                elif counter == 5:
                    r4 = line
                elif counter == 6:
                    r5 = line
            else:   # If K is not alternating, the we read off the Alexander gradings, Maslov gradings, and differential
                if counter == 5:
                    A = string_to_list(line)
                    A = [int(s) for s in A]
                elif counter == 8:
                    M = string_to_list(line)
                    M = [int(s) for s in M]
                elif counter == 11:
                    D = string_to_list(line)
                    D = [string_to_tuple(s) for s in D]
                    D = [(int(a), int(b), int(c)) for (a, b, c) in D]
            if line[0:4] == 'Seif':
                genus = int(line[16:])      # the Seifert genus of K
            if line[0:3] == 'Tau':
                tau = int(line[6:])         # the Ozsvath-Szabo tau invariant of K
            if line[0:4] == 'Epsi':
                epsilon = int(line[10:])    # epsilon of K

                reading_complex = False     # after epsilon, we are done reading what we need from this complex
                knot_counter += 1

                if is_alternating:
                    if tau == 0:
                        tau_counter += 1
                    if epsilon == 0:
                        epsilon_counter += 1
                    if genus == 2 and epsilon == 0:
                        possible_counter += 1
                        print filename + '-' + str(knot_counter),

                        while '(' in r1:
                            r1 = r1[:-1]
                        r1 = int(r1)
                        while '(' in r2:
                            r2 = r2[:-1]
                        r2 = int(r2)
                        while '(' in r3:
                            r3 = r3[:-1]
                        r3 = int(r3)
                        while '(' in r4:
                            r4 = r4[:-1]
                        r4 = int(r4)
                        while '(' in r5:
                            r5 = r5[:-1]
                        r5 = int(r5)

                        if not (r1 == r5 and r2 == r4):
                            print '**** error: not symmetric *****'
                            print 'knot_counter = ', knot_counter
                            print '*******************************'
                        if r2 == 4*r1 and r3 == 6*r1+1:
                            k = r1
                            if k > 1:
                                condensed_curve = '[0], (1,0)^' + str(k) + ', (0,0)^' + str(2*k) + ', (-1,0)^' + str(k)
                            if k == 1:
                                condensed_curve = '[0], (1,0), (0,0)^2, (-1,0)'

                            print condensed_curve
                            f_output.write(filename + '-' + str(knot_counter) + ' ' + condensed_curve + '\n')
                            f_output_short.write(filename + '-' + str(knot_counter) + '\n')
                            definite_counter += 1


                else:
                    delta_gradings = [A[i] - M[i] for i in range(len(A))]
                    thickness = max(delta_gradings) - min(delta_gradings)
                    if thickness > max_thickness:
                        max_thickness = thickness

                    if thickness > 2:
                        print 'thickness is ', thickness, '   ', filename+str(knot_counter)

                    if tau == 0:
                        tau_counter += 1
                    if epsilon == 0:
                        epsilon_counter += 1
                    if genus == 2 and epsilon == 0:
                        possible_counter += 1
                        print filename + '-' + str(knot_counter),

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

                        if thickness == 0:
                            max_grading = max(A)
                            generators_by_grading = [0 for i in range(2*max_grading + 1)]

                            for item in A:
                                generators_by_grading[max_grading - item] += 1
                            # the ith entry in the list generators_by_grading is the number of generators with grading max_grading-i
                            k =generators_by_grading[0]
                            if generators_by_grading == [k, 4*k, 6*k+1, 4*k, k]:
                                if k > 1:
                                    condensed_curve = '[0], (1,0)^'+str(k)+', (0,0)^'+str(2*k)+', (-1,0)^'+str(k)
                                if k == 1:
                                    condensed_curve = '[0], (1,0), (0,0)^2, (-1,0)'
                                print condensed_curve
                                f_output.write(filename + '-' + str(knot_counter) + ' ' + condensed_curve + '\n')
                                f_output_short.write(filename + '-' + str(knot_counter) + '\n')
                                definite_counter += 1


                        else:
                            train_track = CFK.train_track()
                            print 'computing: ' + filename + '-' + str(knot_counter), ' # generators = ', len(train_track.alexander_gradings)

                            condensed_curve = train_track.condensed_curve()

                            if not train_track.other_components == []:
                                print
                                print
                                print '**** example with non figure 8 component ****'
                                print condensed_curve
                                print


                            used_delta_gradings = range(min(delta_gradings), max(delta_gradings)+1)
                            is_possible_example = True

                            for delta in used_delta_gradings:
                                delta_fig8s = []
                                for item in train_track.fig8_components:
                                    if item[1] == delta:
                                        delta_fig8s.append(item)
                                delta_fig8s = [item[0] for item in delta_fig8s]

                                if len(delta_fig8s) > 0:
                                    if (max(delta_fig8s) == 1 and min(delta_fig8s) == -1):
                                        n1 = 0
                                        n0 = 0
                                        for item in delta_fig8s:
                                            if item == 0:
                                                n0 += 1
                                            if item == 1:
                                                n1 += 1
                                        if not n0 == 2*n1:
                                            is_possible_example = False
                                    else:
                                        is_possible_example = False

                            if is_possible_example:
                                print condensed_curve
                                f_output.write(filename + '-' + str(knot_counter) + ' ' + condensed_curve + '\n')
                                f_output_short.write(filename + '-' + str(knot_counter) + '\n')
                                definite_counter += 1

            counter += 1
        line = f_input.readline()

    print
    print 'knots checked: ', knot_counter
    f_output.write('\n')
    f_output.write('knots checked: ' + str(knot_counter) + '\n')

    print '# with tau = 0: ', tau_counter
    f_output.write('# with tau = 0: ' + str(tau_counter) + '\n')
    print '# with e = 0: ', epsilon_counter
    f_output.write('# with e = 0: ' + str(epsilon_counter) + '\n')
    print '# with g=2,e=0: ', possible_counter
    f_output.write('# with g=2,e=0: ' + str(possible_counter) + '\n')
    print '# still unknown: ', definite_counter
    f_output.write('# still unknown: ' + str(definite_counter) + '\n')
    print 'max thickenss: ', max_thickness
    f_output.write('max_thickness: ' + str(max_thickness) + '\n')

    f_input.close()
    f_output.close()
    f_output_short.close()


#### Uncomment the relevant line below to find all potential pairs of cosmetic surgeries on the given set of knots
#### Note that there must be a corresponding input file with the same name in the input_complex directory

# check_conjecture_for_knots('3')
# check_conjecture_for_knots('4')
# check_conjecture_for_knots('5')
# check_conjecture_for_knots('6')
# check_conjecture_for_knots('7')
# check_conjecture_for_knots('8')
# check_conjecture_for_knots('9')
# check_conjecture_for_knots('10')
# check_conjecture_for_knots('11a')
# check_conjecture_for_knots('11n')
# check_conjecture_for_knots('12a')
check_conjecture_for_knots('12n')
# check_conjecture_for_knots('13a')
# check_conjecture_for_knots('13n')
# check_conjecture_for_knots('14a')
# check_conjecture_for_knots('14n')
# check_conjecture_for_knots('15a')
# check_conjecture_for_knots('15n-0')
# check_conjecture_for_knots('15n-1')
# check_conjecture_for_knots('15n-2')
# check_conjecture_for_knots('15n-3')
# check_conjecture_for_knots('15n-4')
# check_conjecture_for_knots('15n-5')
# check_conjecture_for_knots('15n-6')
# check_conjecture_for_knots('15n-7')
# check_conjecture_for_knots('15n-8')
# check_conjecture_for_knots('15n-9')
# check_conjecture_for_knots('15n-10')
# check_conjecture_for_knots('15n-11')
# check_conjecture_for_knots('15n-12')
# check_conjecture_for_knots('15n-13')
# check_conjecture_for_knots('15n-14')
# check_conjecture_for_knots('15n-15')
# check_conjecture_for_knots('15n-16')
# check_conjecture_for_knots('16a')
# check_conjecture_for_knots('16n-0')
# check_conjecture_for_knots('16n-1')
# check_conjecture_for_knots('16n-2')
# check_conjecture_for_knots('16n-3')
# check_conjecture_for_knots('16n-4')
# check_conjecture_for_knots('16n-5')
# check_conjecture_for_knots('16n-6')
# check_conjecture_for_knots('16n-7')
# check_conjecture_for_knots('16n-8')
# check_conjecture_for_knots('16n-9')
# check_conjecture_for_knots('16n-10')