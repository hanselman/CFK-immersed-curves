# Defines a class BifilteredComplex, which is used to encode the UV = 0 version of the bifiltered
# complex CFK(K), working over Z/2Z coefficients. The class contains a function, train_track, which
# takes an instance of BifliteredComplex and converts it to a KnotTrainTrack, the class used to represent
# a collection of immersed curves in the punctured cylinder together with a collection of "crossover arrows"


from KnotTrainTrack import *

def reduce_list_mod2(lst):
    '''returns a copy of the list with identical pairs of elements removed
    until each element appears 0 or 1 times'''
    result = []
    for item in lst:
        if item in result:
            result.remove(item)
        else:
            result.append(item)
    return result

class BifilteredComplex:

    def __init__(self, alexander_gradings, maslov_gradings, differential):
        '''alexander and maslov gradings are lists of integers (one for each generator)
        differential is list of lists. The ith list has entries of the the form
        (j, m, n) indicating a differential from the ith generator to the jth generator weighted by U^m V^n

        for now: assumes generators are sorted by decreasing alexander grading
        assumes that the complex is reduced
        assumes that the differential U-V weights are consistent (could check this)'''

        if len(alexander_gradings) != len(maslov_gradings) or len(maslov_gradings) != len(differential):
            print 'error: length mismatch between alexander grading, maslov grading, and differential'
        self.alexander_gradings = alexander_gradings
        self.maslov_gradings = maslov_gradings
        self.differential = differential
        self.num_generators = len(alexander_gradings)
        self.generators = range(self.num_generators)
        self.max_alex = max(alexander_gradings)
        self.min_alex = min(alexander_gradings)
        self.name = 'Test'
        self.fig8_components = []

        self.sort_by_alexander_grading()

        self.num_gen_by_grading = {}  # dictionary: value at an integer n is number of generators with alexander grading n
        #        self.generators = []    #list of tuples (a,n), a is alexander grading, n is order within alexander grading
        for i in range(self.num_generators):
            if self.alexander_gradings[i] in self.num_gen_by_grading.keys():
                self.num_gen_by_grading[self.alexander_gradings[i]] += 1
            else:
                self.num_gen_by_grading[self.alexander_gradings[i]] = 1


    ##            if i == 0:
    ###                self.generators.append( (self.alexander_gradings[i], 0) )
    ##                self.num_gen_by_grading[ self.alexander_gradings[i] ] = 1
    ##            elif self.generators[-1][0] == self.alexander_gradings[i]:
    ###                self.generators.append( (self.alexander_gradings[i], self.generators[-1][1]+1))
    ##                self.num_gen_by_grading[ self.alexander_gradings[i] ] += 1
    ##            else:
    ###                self.generators.append( (self.alexander_gradings[i], 0) )
    ##                self.num_gen_by_grading[ self.alexander_gradings[i] ] = 1
    ##

    ##    def sort_by_alexander_grading(self):
    ##        '''rearranges geneartors (indexed by integers) so that they are in order of decreasing alexander grading'''
    ##        generators = range(self.num_generators)
    ##        generators.sort(key = lambda x: -self.alexander_gradings[x])
    ##        new_A = [self.alexander_gradings[ generators[i] ] for i in range(self.num_generators)]
    ##        new_M = [self.maslov_gradings[ generators[i] ] for i in range(self.num_generators)]
    ##        new_diff = [self.differential[ generators[i] ] for i in range(self.num_generators)]
    ##        for i in range(len(new_diff)):
    ##            for j in range(len(new_diff[i])):
    ##                (g,m,n) = new_diff[i][j]
    ##                new_g = generators.index(g)
    ##                new_diff[i][j] = (new_g, m, n)
    ##        self.alexander_gradings = new_A
    ##        self.maslov_gradings = new_M
    ##        self.differential = new_diff

    def find_fig8_components(self):
        num_generators = self.num_generators

    def sort_by_alexander_grading(self):
        '''rearranges geneartors (indexed by integers) so that they are in order of decreasing alexander grading'''
        self.generators.sort(key=lambda x: -self.alexander_gradings[x])


    def d_squared(self):
        '''returns True if d^2 = 0, False otherwise, prints d^2'''
        result = []
        for i in range(len(self.differential)):
            this_result = []
            for (g, m, n) in self.differential[i]:
                for (g2, m2, n2) in self.differential[g]:
                    new_term = (g2, m + m2, n + n2)
                    if new_term in this_result:
                        this_result.remove(new_term)
                    else:
                        this_result.append(new_term)

            result.append(this_result)

        is_all_zero = True
        for i in range(len(result)):
            is_zero = True
            for (g, m, n) in result[i]:
                if is_zero:
                    print i, '--> ',
                if not is_zero:
                    print ' + ',
                is_zero = False
                is_all_zero = False
                print 'U^', m, ' V^', n, ' g_', g,
            if not is_zero:
                print
        return is_all_zero


    def is_thin(self):
        x = [self.alexander_gradings[i] - self.maslov_gradings[i] for i in range(self.num_generators)]
        if max(x) == min(x):
            return True
        else:
            return False


    def find_vertical_basis(self):
        vertical_diff = []
        for i in self.generators:  # We assume that the list self.generators is sorted in order of decreasing alexander grading
            new_diff = []
            for j in range(len(self.differential[i])):
                if self.differential[i][j][1] == 0:
                    new_diff.append(self.generators.index(self.differential[i][j][0]))
            vertical_diff.append(new_diff)
        ### vertical_diff is a list of lists of integers.
        ### Each entry corresponds to the rows of a matrix over Z mod 2
        ### the integers tell you the nonzero entries of the row
        ### The indexing of the generators for this matrix refers to the position in the ordered list self.generators, not necessarily to the integer representing the generator
        ### since the complex is reduced, the matrix is strictly block upper triangular

        basis_changes = []
        col = self.num_gen_by_grading[self.max_alex]
        row = col
        while not col >= self.num_generators:
            ##            if col%100 == 0:
            ##                print col
            if not col in vertical_diff[row]:
                row -= 1
                if row < 0:
                    col += 1
                    row = col
            else:  # we have found the lowest row with a nonzero entry in col
                for r in range(0,
                               row):  # for all rows above, we add this row if necessary to clear out any other nonzero entries in col
                    if col in vertical_diff[r]:
                        vertical_diff[r] = reduce_list_mod2(vertical_diff[r] + vertical_diff[row])
                        basis_changes.append((r, row))
                #                        print 'left basis change ', (r, row)
                # note: no column operations are needed, since row < col so these operation would be on columns to the left of col, which are either zero or already dealt with

                for c in vertical_diff[row]:
                    if c > col:  # for all columns c to the right of col which have a nonzero entry in row
                        # add row c to row col
                        vertical_diff[col] = reduce_list_mod2(vertical_diff[col] + vertical_diff[c])
                        basis_changes.append((col, c))
                #                        print 'left basis change ', (col, c)
                # add column col (which only has one nonzero entry, in row) to columns c above
                vertical_diff[row] = [col]
                # now the row row and column col have a single (shared) nonzero entry
                col += 1
                row = col

        vert_matching = []
        for i in range(self.num_generators):
            if not vertical_diff[i] == []:
                if len(vertical_diff[i]) > 1:
                    print 'error: in vertically simplified basis, more than one arrow out of generator ', str(i)
                vert_matching.append((i, vertical_diff[i][0]))
        vert_matched_generators = []
        for (i, j) in vert_matching:
            vert_matched_generators.append(i)
            vert_matched_generators.append(j)
        if not len(vert_matched_generators) == self.num_generators - 1:
            print 'error: incorrect number of matched generators in vertically simplified basis'
        vert_unmatched = range(self.num_generators)
        for i in vert_matched_generators:
            vert_unmatched.remove(i)
        if not len(vert_unmatched) == 1:
            print 'error: more than 1 unmatched generator in vertical basis'
        vert_unmatched_generator = vert_unmatched[0]

        # the matchings here have the form (i,j) where i and j represent the position of the generator in the list self.generators
        # basis changes have the same form, where i and j corresponds to a change of basis
        # [ith gen] --> [ith gen] + [jth gen], again with indices meaning position in the list self.generators
        # the list of basis changes is given in the order done to get from the original basis to a vertically simplified one
        # fror the corresponding arrows in a train track, this order moves outward from the middle.
        ##        print
        ##        print 'vert basis changes = ', basis_changes
        ##        print
        return (vert_matching, vert_unmatched_generator, basis_changes)


    def find_horizontal_basis(self):
        horizontal_diff = []
        for i in self.generators:  # We assume that the list self.generators is sorted in order of decreasing alexander grading
            new_diff = []
            for j in range(len(self.differential[i])):
                if self.differential[i][j][2] == 0:
                    new_diff.append(self.generators.index(self.differential[i][j][0]))
            horizontal_diff.append(new_diff)
        ### horizontal_diff is a list of lists of integers.
        ### Each entry corresponds to the rows of a matrix over Z mod 2
        ### the integers tell you the nonzero entries of the row
        ### The indexing of the generators for this matrix refers to the position in the ordered list self.generators, not necessarily to the integer representing the generator
        ### since the complex is reduced, the matrix is strictly block lower triangular

        basis_changes = []
        col = self.num_generators - 1 - self.num_gen_by_grading[self.min_alex]
        # this is the rightmost column excluding those of the lowest alexander grading
        row = col
        while not col < 0:
            if not col in horizontal_diff[row]:
                row += 1
                if row >= self.num_generators:
                    col -= 1
                    row = col
            else:  # we have found the highest row with a nonzero entry in col
                for r in range(row + 1,
                               self.num_generators):  # for all rows below, we add this row if necessary to clear out any other nonzero entries in col
                    if col in horizontal_diff[r]:
                        horizontal_diff[r] = reduce_list_mod2(horizontal_diff[r] + horizontal_diff[row])
                        basis_changes.append((r, row))
                ##                        print 'right basis change ', (r, row)
                # note: no column operations are needed, since row > col so these operation would be on columns to the right of col, which are either zero or already dealt with

                for c in horizontal_diff[row]:
                    if c < col:  # for all columns c to the left of col which have a nonzero entry in row
                        # add row c to row col
                        horizontal_diff[col] = reduce_list_mod2(horizontal_diff[col] + horizontal_diff[c])
                        basis_changes.append((col, c))
                ##                        print 'right basis change ', (col, c)
                # add column col (which only has one nonzero entry, in row) to columns c above
                horizontal_diff[row] = [col]
                # now the row row and column col have a single (shared) nonzero entry
                col -= 1
                row = col

        horiz_matching = []
        for i in range(self.num_generators):
            if not horizontal_diff[i] == []:
                if len(horizontal_diff[i]) > 1:
                    print 'error: in horizontally simplified basis, more than one arrow out of generator ', str(i)
                horiz_matching.append((i, horizontal_diff[i][0]))
        horiz_matched_generators = []
        for (i, j) in horiz_matching:
            horiz_matched_generators.append(i)
            horiz_matched_generators.append(j)
        if not len(horiz_matched_generators) == self.num_generators - 1:
            print 'error: incorrect number of matched generators in horizontally simplified basis'
        horiz_unmatched = range(self.num_generators)
        for i in horiz_matched_generators:
            horiz_unmatched.remove(i)
        if not len(horiz_unmatched) == 1:
            print 'error: more than 1 unmatched generator in horizontal basis'
        horiz_unmatched_generator = horiz_unmatched[0]

        # the matchings here have the form (i,j) where i and j represent the position of the generator in the list self.generators
        # basis changes have the same form, where i and j corresponds to a change of basis
        # [ith gen] --> [ith gen] + [jth gen], again with indices meaning position in the list self.generators
        # the list of basis changes is given in the order done to get from the original basis to a vertically simplified one
        # fror the corresponding arrows in a train track, this order moves outward from the middle.
        ##        print
        ##        print 'horiz basis changes = ', basis_changes
        ##        print
        return (horiz_matching, horiz_unmatched_generator, basis_changes)


    def train_track(self):
        if verbose:
            print 'finding vertical basis .... ',
        (vert_matching, vert_unmatched, vert_basis_changes) = self.find_vertical_basis()
        if verbose:
            print 'vertical basis found'
            print 'finding vertical basis .... ',
        (horiz_matching, horiz_unmatched, horiz_basis_changes) = self.find_horizontal_basis()
        if verbose:
            print 'horizontal basis found'
        newA = [self.alexander_gradings[i] for i in self.generators]
        newM = [self.maslov_gradings[i] for i in self.generators]

        arrows = [vert_basis_changes[-i - 1] for i in range(len(vert_basis_changes))] + horiz_basis_changes

        return KnotTrainTrack(newA, newM, vert_matching, horiz_matching, arrows, self.name)





def string_to_list(s):
    '''parses a list that is saved as a string.
    Identifies the first '[' and last ']', and in between commas that are not within nested [,] or (,)
    returns a list whose entries are the strings between successive commas'''
    i = 0
    while not s[i] == '[':
        i += 1
    layer = 0
    i += 1
    result = []
    next_entry = ''
    while not (s[i] == ']' and layer == 0):
        if s[i] in ['(','[']:
            layer += 1
        if s[i] in [')',']'] and layer > 0:
            layer -= 1
        if s[i] == ',' and layer == 0:
            result.append(next_entry)
            next_entry = ''
        else:
            next_entry = next_entry + s[i]
        i += 1
    result.append(next_entry)
    return result

def string_to_tuple(s):
    '''parses a tuple that is saved as a string.
    identifies the initial assumes there is a single '(' and a single ')'
    Returns a tuple whose entries are the strings between successive commas'''
    #there can not be commas wihtin an entrie, i.e. no nested lists/tuples are allowed
    i = 0
    while not s[i] == '(':
        i += 1
    i += 1
    result = []
    next_entry = ''
    while not s[i] == ')':
        if s[i] == ',':
            result.append(next_entry)
            next_entry = ''
        else:
            next_entry = next_entry + s[i]
        i += 1
    result.append(next_entry)
    return tuple(result)


def CFK_from_file(filename):
    f = open('input_complex/' + filename + '.txt', 'r')
    lines = f.readlines()
    A = string_to_list(lines[9])
    A = [int(s) for s in A]
    M = string_to_list(lines[12])
    M = [int(s) for s in M]
    D = string_to_list(lines[15])
    D = [string_to_tuple(s) for s in D]
    D = [(int(a), int(b), int(c)) for (a, b, c) in D]

    f.close()

    diff = []
    for i in range(len(A)):
        diff.append([])

    for (a, b, c) in D:
        # a = init_gen, b = final_gen, c = coefficient (should always be 1 for mod 2)
        if c == 1:
            Upower = (M[b] + 1 - M[a]) / 2
            Vpower = A[a] - A[b] + Upower
            diff[a].append((b, Upower, Vpower))

    result = BifilteredComplex(A, M, diff)
    result.name = filename
    return result

def CFKs_from_file(filename, names, start_index = 0, end_index = 'max'):
    f = open('input_complex/' + filename + '.txt', 'r')
    lines = f.readlines()
    if end_index == 'max':
        end_index = len(lines)
    j = 0
    counter = 0
    result = []
    while j < len(lines) and counter < end_index:
        if lines[j][0:4] == 'Knot':
            if counter >= start_index:
                A = string_to_list(lines[j + 5])
                A = [int(s) for s in A]
                M = string_to_list(lines[j + 8])
                M = [int(s) for s in M]
                D = string_to_list(lines[j + 11])
                D = [string_to_tuple(s) for s in D]
                D = [(int(a), int(b), int(c)) for (a, b, c) in D]

                diff = []
                for i in range(len(A)):
                    diff.append([])

                for (a, b, c) in D:
                    # a = init_gen, b = final_gen, c = coefficient (should always be 1 for mod 2)
                    if c == 1:
                        Upower = (M[b] + 1 - M[a]) / 2
                        Vpower = A[a] - A[b] + Upower
                        diff[a].append((b, Upower, Vpower))

                next_result = BifilteredComplex(A, M, diff)
                next_result.name = names[counter]
                result.append(next_result)
            counter += 1
        j += 1

    f.close()
    return result



