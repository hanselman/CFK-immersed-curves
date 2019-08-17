# Defines the class KnotTrainTrack to encode collection of immersed curves in the infinitely punctured cylinder
# decorated with a collection of crossover arrows, such as the structure that arises from the bifiltered
# chain complex CFK.
# Contains functions implementing the arrow sliding algorithm of Hanselman-Rasmussen-Watson to simplify
# an instance of KnotTrainTrack to a collection of immersed curves
# (*does not currently support local systems on curves, but will print an error if local systems are needed... this has
# not yet happen, and it is an open question whether local systems are ever needed for knots in S^3)


from linear_algebra import *

verbose = False


class Arrow:
    def __init__(self, init_gen, final_gen, train_track):
        self.init_gen = init_gen
        self.final_gen = final_gen
        self.train_track = train_track
        self.alexander_grading = train_track.alexander_gradings[init_gen]

    def left_weight(self):
        return self.train_track.left_weight(self.init_gen, self.final_gen)

    def right_weight(self):
        return self.train_track.right_weight(self.init_gen, self.final_gen)


class KnotTrainTrack:
    def __init__(self, alexander_gradings, maslov_gradings, left_matching, right_matching, arrows, name='NoName'):
        '''alexander_gradings: list of 2n+1 integers (one for each generator). For now: assume this is decreasing
        maslov_gradings: list of 2n+1 integers (one for each generator)
        left/right_matching: list of n tuples of integers in range(2n-1)
        arrows: list of 2-tuples of integers in range(2n-1), with each tuple representing
                the initial and final generator of an arrow. Arrows are ordered from left to right

        *** for now: assumes there is exactly one unmatched generator in each matching (should build in a check for this)'''

        self.num_generators = len(alexander_gradings)
        self.alexander_gradings = alexander_gradings
        self.min_alex = min(alexander_gradings)
        self.max_alex = max(alexander_gradings)
        self.maslov_gradings = maslov_gradings
        self.left_matching = left_matching
        self.right_matching = right_matching
        left_matched = [pair[0] for pair in left_matching] + [pair[1] for pair in left_matching]
        for i in range(self.num_generators):
            if not i in left_matched:
                self.left_unmatched = i
        right_matched = [pair[0] for pair in right_matching] + [pair[1] for pair in right_matching]
        for i in range(self.num_generators):
            if not i in right_matched:
                self.right_unmatched = i

        self.left_match_dict = {}
        for (i, j) in left_matching:
            self.left_match_dict[i] = j
            self.left_match_dict[j] = i
        self.left_match_dict[self.left_unmatched] = None

        self.right_match_dict = {}
        for (i, j) in right_matching:
            self.right_match_dict[i] = j
            self.right_match_dict[j] = i
        self.right_match_dict[self.right_unmatched] = None

        self.name = name
        self.is_simplified = False

        self.arrows = []
        for ar in arrows:
            if alexander_gradings[ar[0]] == alexander_gradings[ar[1]]:  # we keep only arrows between same alexander grading
                self.arrows.append(Arrow(ar[0], ar[1], self))

        self.generators_by_grading = {}
        for alex in range(self.min_alex, self.max_alex + 1):
            self.generators_by_grading[alex] = []
        for i in range(self.num_generators):
            self.generators_by_grading[self.alexander_gradings[i]].append(i)

        self.arrows_by_grading = {}
        for alex in range(self.min_alex, self.max_alex + 1):
            self.arrows_by_grading[alex] = []
        for ar in self.arrows:
            self.arrows_by_grading[ar.alexander_grading].append(ar)

    def right_matching_switch(self, i, j):
        '''modifies self.right_matching and self.right_match_dict corresponding to adding a crossing at
        the right end of generator i and generator j'''
        i2 = self.right_match_dict[i]
        j2 = self.right_match_dict[j]
        self.right_match_dict[i] = j2
        self.right_match_dict[j] = i2
        if not i2 == None:
            self.right_match_dict[i2] = j
        if not j2 == None:
            self.right_match_dict[j2] = i
        if i2 == None:
            self.right_unmatched = j
        if j2 == None:
            self.right_unmatched = i
        new_matching = []
        for (a, b) in self.right_matching:
            if a == i:
                new_matching.append((j, b))
            elif b == i:
                new_matching.append((a, j))
            elif a == j:
                new_matching.append((i, b))
            elif b == j:
                new_matching.append((a, i))
            else:
                new_matching.append((a, b))
        self.right_matching = new_matching

    def left_matching_switch(self, i, j):
        '''modifies self.left_matching and self.left_match_dict corresponding to adding a crossing at
        the left end of generator i and generator j'''
        i2 = self.left_match_dict[i]
        j2 = self.left_match_dict[j]
        self.left_match_dict[i] = j2
        self.left_match_dict[j] = i2
        if not i2 == None:
            self.left_match_dict[i2] = j
        if not j2 == None:
            self.left_match_dict[j2] = i
        if i2 == None:
            self.left_unmatched = j
        if j2 == None:
            self.left_unmatched = i
        new_matching = []
        for (a, b) in self.left_matching:
            if a == i:
                new_matching.append((j, b))
            elif b == i:
                new_matching.append((a, j))
            elif a == j:
                new_matching.append((i, b))
            elif b == j:
                new_matching.append((a, i))
            else:
                new_matching.append((a, b))
        self.left_matching = new_matching

    def left_color(self, gen, depth):
        '''returns the depth-m color of the left end of the segment for generator gen
        returns a length m list of numbers (floats)
        a right turn where the next handle is n units away corresponds to 1/n
        a left turn wehre the next handle is n units away corresponds to -1/n
        going straight (the unmatched edge) corresponds to 0, and the rest of the list is 0
        e.g. a path left[2]-left[1]-right[1]-straight returns [-.5,-1,1,0]

        The relation < is such that (more left path) < (more right path)
        '''

        if depth == 0:
            return []
        elif self.left_match_dict[gen] == None:
            return [0]*depth
        else:
            next_gen = self.left_match_dict[gen]
            grading_difference = self.alexander_gradings[next_gen] - self.alexander_gradings[gen]
            return [1.0/grading_difference]+self.right_color(next_gen, depth-1)


    def right_color(self, gen, depth):
        '''returns the depth-m color of the right end of the segment for generator gen
        returns a length, a sequence of 0 (for left turn), 1 (for straight), or 2 (for right turn)
        e.g. a path left-left-right-stright returns '0021'
        The string stops after an instance of 1, otherwise it has length [depth]
        '''

        if depth == 0:
            return []
        elif self.right_match_dict[gen] == None:
            return [0]*depth
        else:
            next_gen = self.right_match_dict[gen]
            grading_difference = self.alexander_gradings[gen] - self.alexander_gradings[next_gen]
            return [1.0 / grading_difference] + self.left_color(next_gen, depth - 1)

    def right_color_old(self, gen, depth):
        '''returns the depth-m color of the right end of the segment for generator gen
        returns a length, a sequence of 0 (for left turn), 1 (for straight), or 2 (for right turn)
        e.g. a path left-left-right-stright returns '0021'
        The string stops after an instance of 1, otherwise it has length [depth]
        '''

        if depth == 0:
            return ''
        elif self.right_match_dict[gen] == None:
            return '1'
        else:
            next_gen = self.right_match_dict[gen]
            if self.alexander_gradings[next_gen] > self.alexander_gradings[gen]:
                return '0' + self.left_color(next_gen, depth - 1)
            else:
                return '2' + self.left_color(next_gen, depth - 1)

    def push_left(self, arrow_to_move):
        alex = arrow_to_move.alexander_grading
        index = self.arrows_by_grading[alex].index(arrow_to_move)
        # we assume that index is not 0
        arrow_to_pass = self.arrows_by_grading[alex][index - 1]
        (ap1, ap2) = (arrow_to_pass.init_gen, arrow_to_pass.final_gen)
        (am1, am2) = (arrow_to_move.init_gen, arrow_to_move.final_gen)
        if (am1, am2) == (ap1, ap2):  # if the two arrows are the same, they cancel each other
            self.arrows.remove(arrow_to_move)
            self.arrows_by_grading[alex].remove(arrow_to_move)
            self.arrows.remove(arrow_to_pass)
            self.arrows_by_grading[alex].remove(arrow_to_pass)
        elif (am1, am2) == (ap2,
                            ap1):  # if the two arrows have the same endpoints and opposite direction, we replace op_arrow - arrow with arrow - crossing
            self.arrows.remove(arrow_to_pass)
            self.arrows_by_grading[alex].remove(arrow_to_pass)  # now arrow_to_move is at position index-1
            for arrow in self.arrows_by_grading[alex][
                         index:]:  # we push the crossing to the right, which switches am1 and am2 for all arrows to the right of arrow_to_move
                if arrow.init_gen == am1:
                    arrow.init_gen = am2
                elif arrow.init_gen == am2:
                    arrow.init_gen = am1
                if arrow.final_gen == am1:
                    arrow.final_gen = am2
                elif arrow.final_gen == am2:
                    arrow.final_gen = am1
            self.right_matching_switch(am1, am2)
        elif am2 == ap1:  # if end of arrow_to_move is start of arrow_to_pass, we have to add a new concatenation arrow
            new_arrow = Arrow(am1, ap2, self)
            new_arrow_list = self.arrows_by_grading[alex][:index - 1] + [arrow_to_move, new_arrow, arrow_to_pass] + \
                             self.arrows_by_grading[alex][index + 1:]
            self.arrows_by_grading[alex] = new_arrow_list
            self.arrows.append(new_arrow)
        elif am1 == ap2:  # if start of arrow_to_move is end of arrow_to_pass, we have to add a new concatenation arrow
            new_arrow = Arrow(ap1, am2, self)
            new_arrow_list = self.arrows_by_grading[alex][:index - 1] + [arrow_to_move, new_arrow, arrow_to_pass] + \
                             self.arrows_by_grading[alex][index + 1:]
            self.arrows_by_grading[alex] = new_arrow_list
            self.arrows.append(new_arrow)
        else:  # otherwise, the arrows pass each other with no interaction
            new_arrow_list = self.arrows_by_grading[alex][:index - 1] + [arrow_to_move, arrow_to_pass] + \
                             self.arrows_by_grading[alex][index + 1:]
            self.arrows_by_grading[alex] = new_arrow_list

    def push_right(self, arrow_to_move):
        alex = arrow_to_move.alexander_grading
        index = self.arrows_by_grading[alex].index(arrow_to_move)
        # we assume that index is not len(arrows_by_grading[alex])-1
        arrow_to_pass = self.arrows_by_grading[alex][index + 1]
        (ap1, ap2) = (arrow_to_pass.init_gen, arrow_to_pass.final_gen)
        (am1, am2) = (arrow_to_move.init_gen, arrow_to_move.final_gen)
        if (am1, am2) == (ap1, ap2):  # if the two arrows are the same, they cancel each other
            self.arrows.remove(arrow_to_move)
            self.arrows_by_grading[alex].remove(arrow_to_move)
            self.arrows.remove(arrow_to_pass)
            self.arrows_by_grading[alex].remove(arrow_to_pass)
        elif (am1, am2) == (ap2,
                            ap1):  # if the two arrows have the same endpoints and opposite direction, we replace arrow - op_arrow with crossing - arrow
            self.arrows.remove(arrow_to_pass)
            self.arrows_by_grading[alex].remove(arrow_to_pass)
            for arrow in self.arrows_by_grading[alex][
                         :index]:  # we push the crossing to the left, which switches am1 and am2 for all arrows to the left of arrow_to_move
                if arrow.init_gen == am1:
                    arrow.init_gen = am2
                elif arrow.init_gen == am2:
                    arrow.init_gen = am1
                if arrow.final_gen == am1:
                    arrow.final_gen = am2
                elif arrow.final_gen == am2:
                    arrow.final_gen = am1
            self.left_matching_switch(arrow_to_move.init_gen, arrow_to_move.final_gen)
        elif am2 == ap1:  # if end of arrow_to_move is start of arrow_to_pass, we have to add a new concatenation arrow
            new_arrow = Arrow(am1, ap2, self)
            new_arrow_list = self.arrows_by_grading[alex][:index] + [arrow_to_pass, new_arrow, arrow_to_move] + \
                             self.arrows_by_grading[alex][index + 2:]
            self.arrows_by_grading[alex] = new_arrow_list
            self.arrows.append(new_arrow)
        elif am1 == ap2:  # if start of arrow_to_move is end of arrow_to_pass, we have to add a new concatenation arrow
            new_arrow = Arrow(ap1, am2, self)
            new_arrow_list = self.arrows_by_grading[alex][:index] + [arrow_to_pass, new_arrow, arrow_to_move] + \
                             self.arrows_by_grading[alex][index + 2:]
            self.arrows_by_grading[alex] = new_arrow_list
            self.arrows.append(new_arrow)
        else:  # otherwise, the arrows pass each other with no interaction
            new_arrow_list = self.arrows_by_grading[alex][:index] + [arrow_to_pass, arrow_to_move] + \
                             self.arrows_by_grading[alex][index + 2:]
            self.arrows_by_grading[alex] = new_arrow_list

    def depth(self):
        '''returns the absolute value of minimum weight over both sides of all arrows
        or returns None if there are no arrows or returns the string 'infty' if all arrows have infinite weight'''
        if len(self.arrows) == 0:
            return None
        weights = []
        for arrow in self.arrows:
            (i, j) = (arrow.left_weight(), arrow.right_weight())
            if not i == 'infty':
                weights.append(i)
            if not j == 'infty':
                weights.append(j)
        if len(weights) == 0:
            ##            print '*******************'
            ##            print self.name, ' has a nontrivial local system'
            ##            print '*******************'
            return 'infty'
        else:
            return min([abs(x) for x in weights])



    def sort_arrows_left(self, alexander_grading, first_arrow, last_arrow):
        '''sorts a subset of the arrows at a given alexander_grading, specifically the range
        from first_arrow to last_arrow, so that lower left_weight arrows are on the right
        for arrows with same left weight, they are ordered by init_gen and then by final_gen'''

        if not first_arrow <= last_arrow:
            return
            print 'error'
        num_arrows_on_left = first_arrow
        num_arrows_on_right = len(self.arrows_by_grading[alexander_grading]) - last_arrow - 1
        i = len(self.arrows_by_grading[alexander_grading]) - num_arrows_on_right - 2
        while not i < num_arrows_on_left:
            arrow = self.arrows_by_grading[alexander_grading][i]
            next_arrow = self.arrows_by_grading[alexander_grading][i + 1]

            if arrow.left_weight() < next_arrow.left_weight():
                self.push_right(arrow)
                i = min(i + 2, len(self.arrows_by_grading[alexander_grading]) - num_arrows_on_right - 2)
            elif arrow.left_weight() == next_arrow.left_weight() and arrow.init_gen > next_arrow.init_gen:
                self.push_right(arrow)
                i = min(i + 2, len(self.arrows_by_grading[alexander_grading]) - num_arrows_on_right - 2)
            elif arrow.left_weight() == next_arrow.left_weight() and arrow.init_gen == next_arrow.init_gen and arrow.final_gen >= next_arrow.final_gen:
                self.push_right(arrow)
                i = min(i + 2, len(self.arrows_by_grading[alexander_grading]) - num_arrows_on_right - 2)
            else:
                i -= 1

    def sort_arrows_right(self, alexander_grading, first_arrow, last_arrow):
        '''sorts a subset of the arrows at a given alexander_grading, specifically the range
        from first_arrow to last_arrow, so that lower right_weight arrows are on the left
        for arrows with same right weight, they are ordered by init_gen and then by final_gen'''

        if not first_arrow <= last_arrow:
            return
            print 'error'
        num_arrows_on_left = first_arrow
        num_arrows_on_right = len(self.arrows_by_grading[alexander_grading]) - last_arrow - 1
        i = num_arrows_on_left + 1
        while not i > len(self.arrows_by_grading[alexander_grading]) - num_arrows_on_right - 1:
            arrow = self.arrows_by_grading[alexander_grading][i]
            next_arrow = self.arrows_by_grading[alexander_grading][i - 1]

            right_weight = arrow.right_weight()
            if right_weight < 0:
                right_weight = -right_weight
            next_right_weight = next_arrow.right_weight()
            if next_right_weight < 0:
                next_right_weight = -next_right_weight

            if right_weight < next_right_weight:
                self.push_left(arrow)
                i = max(i - 2, num_arrows_on_left + 1)
            elif right_weight == next_right_weight and arrow.init_gen > next_arrow.init_gen:
                self.push_left(arrow)
                i = max(i - 2, num_arrows_on_left + 1)
            elif right_weight == next_right_weight and arrow.init_gen == next_arrow.init_gen and arrow.final_gen >= next_arrow.final_gen:
                self.push_left(arrow)
                i = max(i - 2, num_arrows_on_left + 1)
            else:
                i += 1

    def simplify(self):
        self.is_simplified = True
        m = self.depth()
        while not m in [None, 'infty']:
            m = self.depth()
            if verbose:
                print '# arrows = ', len(self.arrows)

            for alex in range(self.min_alex, self.max_alex + 1):
                self.simplify_handle(alex, m, range(self.min_alex, alex))
            m = self.depth()

        if verbose:
            print 'finding components...'

        self.find_components()

        if verbose:
            print 'done simplifying'
            print
            
    def simplify_handle(self, alexander_grading, min_weight, already_simplified_handles):
        '''For the 1-handle at height alexander_grading, performs arrow slides until there are arrows
        with left or right weight <= min_weight.
        already_simplified_handles is a list of heights of other handles which we assume are simplified to this level.
        in the process we do not at any arrows with weight < min_weight to any other handles and do not add any
        arrows with weight <= min_weight to already simplified handles.'''

        # confirm that already_simplified_handles is correct
        for old_alex in already_simplified_handles:
            for ar in self.arrows_by_grading[old_alex]:
                if not ar.left_weight() == 'infty':
                    if abs(ar.left_weight()) <= min_weight or abs(ar.right_weight()) <= min_weight:
                        print 'error: ', old_alex, ' not fully simplified'

        arrows = self.arrows_by_grading[alexander_grading]
        if len(arrows) == 0:
            return
        generators = self.generators_by_grading[alexander_grading]
        n = len(generators)


        M = identity_matrix(n)
        for ar in arrows[::-1]:   #add arrows one at a time from right to left
            i = generators.index(ar.init_gen)
            j = generators.index(ar.final_gen)
            M[i] = list_sum_mod2(M[i],M[j])   #add the jth row to the ith row

        ## M now represents the parallel strands and arrows, with respect to the ordering given by the list [generators]

        left_colors = {}
        for gen in generators:
            left_colors[gen] = self.left_color(gen, min_weight+1)
        right_colors = {}
        for gen in generators:
            right_colors[gen] = self.right_color(gen, min_weight + 1)

        left_generators = sorted(generators, key = lambda x: left_colors[x], reverse = True)
        right_generators = sorted(generators, key=lambda x: right_colors[x])

        M2 = []
        for i in range(len(generators)):
            new_row = []
            shifted_i = generators.index(left_generators[i])
            for j in range(len(generators)):
                shifted_j = generators.index(right_generators[j])
                new_row.append(M[shifted_i][shifted_j])
            M2.append(new_row)

        # print 'left_generators: ', left_generators
        # print 'right_generators: ', right_generators
        # print 'M2 = '
        # print_matrix(M2)

        #The i,j entry of M2 now represents the count of paths from left_generators[i] on left to right_generators[j] on right

        m = min_weight
        if verbose:
            print 'simplifying grading ', alexander_grading, ' at depth ', min_weight, '...'
            # print [(ar.init_gen, ar.final_gen) for ar in self.arrows_by_grading[alexander_grading]]


        (L,P,U) = LPU_decomp_mod2(M2)
        if verbose:
            print 'LPU decomposition found'
            
        P_index = [row.index(1) for row in P]
        # the permutation matrix P takes i on the left to P_index[i] on the right (ie the i, P_index[i] entry of P is 1)

        # print 'L='
        # print_matrix(L)
        # print 'P='
        # print_matrix(P)
        # print 'U='
        # print_matrix(U)
        # print
        # print 'P_index = ', P_index
        # print


        # print 'new left arrows: ',
        pre_left_arrows = decompose_lower_triangular(L)
        left_arrows = [(left_generators[ar[0]],left_generators[ar[1]]) for ar in pre_left_arrows]
        left_arrows = [Arrow(ar[0], ar[1], self) for ar in left_arrows]
        # print [(ar.init_gen, ar.final_gen) for ar in left_arrows]

        # print 'new right arrows: ',
        pre_right_arrows = decompose_upper_triangular(U)
        #slide arrows through permutation P
        pre_right_arrows2 = [(P_index.index(ar[0]),P_index.index(ar[1])) for ar in pre_right_arrows]
        right_arrows = [(left_generators[ar[0]], left_generators[ar[1]]) for ar in pre_right_arrows2]
        right_arrows = [Arrow(ar[0], ar[1], self) for ar in right_arrows]
        # print [(ar.init_gen, ar.final_gen) for ar in right_arrows]

        new_permutation = []
        for i in range(n):
            gen = generators[i]
            left_index = left_generators.index(gen)
            right_index = P_index[left_index]
            new_gen = right_generators[right_index]
            j = generators.index(new_gen)
            new_permutation.append(j)

        #new_permutation is the permutation list, where generator[i] on the left is connected to generator[new_permutation[i]] on the right

        for i1 in range(n):
            j1 = new_permutation[i1]
            if not j1 == i1:
                j2 = i1
                i2 = new_permutation.index(j2)

                self.right_matching_switch(generators[j1], generators[j2])
                # print 'switching right matching: ', generators[j1], generators[j2]
                new_permutation[i1] = j2
                new_permutation[i2] = j1

        ######

        ## remove old arrows, add new ones
        for ar in self.arrows_by_grading[alexander_grading]:
            self.arrows.remove(ar)
        for ar in left_arrows:
            self.arrows.append(ar)
        for ar in right_arrows:
            self.arrows.append(ar)
        self.arrows_by_grading[alexander_grading] = left_arrows + right_arrows

        num_left_arrows = len(left_arrows)
        num_right_arrows = len(right_arrows)
        ## At this point the arrows at the given alexander grading have been sorted so that
        ## all of the left arrows have left weight not equal to -m or -(m+1) and
        ## all of the right arrows have right weight not equal to -m or -(m+1)

        ## Check:
        for ar in self.arrows_by_grading[alexander_grading][:num_left_arrows]:
            if not (ar.left_weight() == 'infty'):
                if ar.left_weight() >= -(m+1) and ar.left_weight() < m:
                    print 'error with left_weights in arrow check'
        for ar in self.arrows_by_grading[alexander_grading][num_left_arrows:]:
            if not (ar.right_weight() == 'infty'):
                if ar.right_weight() >= -(m + 1) and ar.right_weight() < m:
                    print 'error with right_weights in arrow check'
                    print ar.right_weight()
                    print (ar.init_gen, ar.final_gen)
                    for g in range(self.min_alex, self.max_alex+1):
                        print g, self.generators_by_grading[g]
                    print 'left_match = ', self.left_matching
                    print 'right_match = ', self.right_matching
                    print num_left_arrows, num_right_arrows
                    print right_generators
                    print self.right_color(4,2), self.right_color(6,2)
                    print self.left_color(2, 1), self.left_color(3, 1)

                    print

        if verbose:
            print "new arrow configuration assembled"
    
        ## remove all instances of outer weight = +m
        for ar in right_arrows[::-1]:
            if ar.right_weight() == m:
                self.arrows.remove(ar)
                self.arrows_by_grading[alexander_grading].remove(ar)
                num_right_arrows -= 1
        for ar in left_arrows[::-1]:
            if ar.left_weight() == m:
                self.arrows.remove(ar)
                self.arrows_by_grading[alexander_grading].remove(ar)
                num_left_arrows -= 1

        if verbose:
            print "outer weight m arrows removed"

        #for any arrows on right with left_weight = +/- m, push right off of the handle
        i = len(self.arrows_by_grading[alexander_grading]) - 1
        while i >= num_left_arrows and len(self.arrows_by_grading[alexander_grading]) - num_left_arrows > 0:
            # print i, ' out of ',len(self.arrows_by_grading[alexander_grading]), num_left_arrows, num_right_arrows
            if not self.arrows_by_grading[alexander_grading][i].left_weight() in [m, -m]:
                i -= 1
            else:
                arrow_to_move = self.arrows_by_grading[alexander_grading][i]
                if i == len(self.arrows_by_grading[alexander_grading]) - 1:
                    (a1, a2) = (arrow_to_move.init_gen, arrow_to_move.final_gen)
                    (b1, b2) = (self.right_match_dict[a1], self.right_match_dict[a2])
                    new_arrow = Arrow(b1, b2, self)
                    new_alex = new_arrow.alexander_grading

                    self.arrows.remove(arrow_to_move)
                    self.arrows_by_grading[alexander_grading].remove(arrow_to_move)
                    if not (arrow_to_move.right_weight() == m + 1 and new_alex in already_simplified_handles):
                        # print 'pushing arrow right to new handle', new_arrow.init_gen, new_arrow.final_gen
                        self.arrows.append(new_arrow)
                        self.arrows_by_grading[new_alex].append(new_arrow)
                    i -= 1
                else:
                    # print 'pushing right'
                    self.push_right(arrow_to_move)
                    i = min(i+2, len(self.arrows_by_grading[alexander_grading]) - 1)

        num_right_arrows = len(self.arrows_by_grading[alexander_grading])-num_left_arrows
        # for any arrows on left with right_weight = +/- m, push left off of the handle
        i = 0
        while i < len(self.arrows_by_grading[alexander_grading]) - num_right_arrows and len(self.arrows_by_grading[alexander_grading]) - num_right_arrows > 0:
            if not self.arrows_by_grading[alexander_grading][i].right_weight() in [m, -m]:
                i += 1
            else:
                arrow_to_move = self.arrows_by_grading[alexander_grading][i]
                if i == 0:
                    (a1, a2) = (arrow_to_move.init_gen, arrow_to_move.final_gen)
                    (b1, b2) = (self.left_match_dict[a1], self.left_match_dict[a2])
                    new_arrow = Arrow(b1, b2, self)
                    new_alex = new_arrow.alexander_grading

                    self.arrows.remove(arrow_to_move)
                    self.arrows_by_grading[alexander_grading].remove(arrow_to_move)
                    if not (arrow_to_move.left_weight() == m + 1 and new_alex in already_simplified_handles):
                        # print 'pushing arrow left to new handle', new_arrow.init_gen, new_arrow.final_gen
                        self.arrows = [new_arrow] + self.arrows
                        self.arrows_by_grading[new_alex] = [new_arrow] + self.arrows_by_grading[new_alex]
                else:
                    # print 'pushing left'
                    self.push_left(arrow_to_move)
                    i = max(i - 2, 0)

        # check:
        for ar in self.arrows_by_grading[alexander_grading]:
            if not ar.left_weight() == 'infty':
                if min(abs(ar.left_weight()),abs(ar.right_weight())) <= m:
                    print 'error in simplify_handle'

        if verbose:
            print "handle simplified"



    def left_weight(self, i, j, stop_after='max'):
        '''determines the left weight of an arrow from gen i to gen j
        the absolute value of the weight is how many times the arrow must be slid along parallel arcs
        (starting moving leftward) before the curves the arrow connects diverges.
        The sign is + if the arrow moves left to right when they diverge, and - if it moves right to left.
        if the curves are parallel, the weight is the string 'infty'  '''

        if stop_after == 'max':
            stop_after = self.num_generators
        if stop_after == 0:
            return 'infty'

        current_alex = self.alexander_gradings[i]
        i2 = self.left_match_dict[i]
        j2 = self.left_match_dict[j]
        if not i2 == None:
            alex_i2 = self.alexander_gradings[i2]
        if not j2 == None:
            alex_j2 = self.alexander_gradings[j2]

        if i2 == None and alex_j2 > current_alex:  # in this case, the strand from j diverges right
            return +1
        elif i2 == None and alex_j2 < current_alex:  # here j diverges left
            return -1
        elif j2 == None and alex_i2 < current_alex:  # here i diverges left
            return +1
        elif j2 == None and alex_i2 > current_alex:  # i diverges right
            return -1
        elif alex_i2 < current_alex and alex_j2 > current_alex:
            return +1
        elif alex_i2 > current_alex and alex_j2 < current_alex:
            return -1
        elif alex_i2 > alex_j2:  # in this case both curves either bend up or bend down, and i is to the left
            return +1
        elif alex_i2 < alex_j2:
            return -1
        elif alex_i2 == alex_j2:
            recursive_step = self.right_weight(i2, j2, stop_after - 1)
            if recursive_step == 'infty':
                return 'infty'
            elif recursive_step > 0:
                return recursive_step + 1
            elif recursive_step < 0:
                return recursive_step - 1
        else:
            print 'error: unexpected case'

    def right_weight(self, i, j, stop_after='max'):
        '''determines the right weight of an arrow from gen i to gen j
        the absolute value of the weight is how many times the arrow must be slid along parallel arcs
        (starting moving rightward) before the curves the arrow connects diverges.
        The sign is + if the arrow moves left to right when they diverge, and - if it moves right to left.
        if the curves are parallel, the weight is the string 'infty'  '''

        if stop_after == 'max':
            stop_after = self.num_generators
        if stop_after == 0:
            return 'infty'

        current_alex = self.alexander_gradings[i]
        i2 = self.right_match_dict[i]
        j2 = self.right_match_dict[j]
        if not i2 == None:
            alex_i2 = self.alexander_gradings[i2]
        if not j2 == None:
            alex_j2 = self.alexander_gradings[j2]

        if i2 == None and alex_j2 < current_alex:  # in this case, the strand from j diverges right
            return +1
        elif i2 == None and alex_j2 > current_alex:  # here j diverges left
            return -1
        elif j2 == None and alex_i2 > current_alex:  # here i diverges left
            return +1
        elif j2 == None and alex_i2 < current_alex:  # i diverges right
            return -1
        elif alex_i2 < current_alex and alex_j2 > current_alex:
            return -1
        elif alex_i2 > current_alex and alex_j2 < current_alex:
            return +1
        elif alex_i2 > alex_j2:  # in this case both curves either bend up or bend down, and i is to the right
            return -1
        elif alex_i2 < alex_j2:
            return +1
        elif alex_i2 == alex_j2:
            recursive_step = self.left_weight(i2, j2, stop_after - 1)
            if recursive_step == 'infty':
                return 'infty'
            elif recursive_step > 0:
                return recursive_step + 1
            elif recursive_step < 0:
                return recursive_step - 1
        else:
            print 'error: unexpected case'

    def find_components(self):
        '''assumes there are only infty weight arrows left, computes the following attributes
        self.nontrivial_component: list of integers, the alexander gradings hit by the unique nontrivial curve
        self.fig8_components: list of 2-tuples of integers, each is the alex and delta grading of the middle of a fig8 component
        self.other_components: a list of any remaining components, each given as a list of integers representing generators'''
        unused_generators = range(self.num_generators)
        
        nontrivial_component = [self.left_unmatched]
        unused_generators.remove(self.left_unmatched)
        while not nontrivial_component[-1] == self.right_unmatched:
            if len(nontrivial_component) % 2 == 1:
                next_gen = self.right_match_dict[nontrivial_component[-1]]
                nontrivial_component.append(next_gen)
                unused_generators.remove(next_gen)
            else:
                next_gen = self.left_match_dict[nontrivial_component[-1]]
                nontrivial_component.append(next_gen)
                unused_generators.remove(next_gen)
        self.nontrivial_component = [self.alexander_gradings[item] for item in nontrivial_component]

        
        self.fig8_components = []
        other_components = []
        while len(unused_generators) > 0:
            next_component = [unused_generators.pop(0)]
            while not next_component[-1] == self.left_match_dict[next_component[0]]:
                if len(next_component) % 2 == 1:
                    next_gen = self.right_match_dict[next_component[-1]]
                    next_component.append(next_gen)
                    unused_generators.remove(next_gen)
                else:
                    next_gen = self.left_match_dict[next_component[-1]]
                    next_component.append(next_gen)
                    unused_generators.remove(next_gen)

            is_8 = True
            top = max([self.alexander_gradings[item] for item in next_component])
            count = [0,0,0]            
            for gen in next_component:
                a = self.alexander_gradings[gen]
                if a == top:
                    count[0] += 1
                elif a == top-1:
                    count[1] += 1
                elif a == top-2:
                    count[2] += 1
                else:
                    is_8 = False
            if not count[2] == count[0]:
                is_8 = False
            if not count[1] == 2*count[0]:
                is_8 = False
                
            if is_8:
                alex = top-1
                delta = self.alexander_gradings[next_component[0]] - self.maslov_gradings[next_component[0]]
                for i in range(count[0]):
                    self.fig8_components.append( (alex,delta) )
            else:
                other_components.append(next_component)
                print
                print
                print '******** non-8 closed component in ', self.name, '********'
                print
                print
        self.other_components = other_components

    def condensed_curve(self):
        if not self.is_simplified:
            self.simplify()
        self.find_components()
        self.fig8_components.sort()
        self.fig8_components.reverse()

        condensed_fig8_components = []

        for fig in self.fig8_components[:1]:
            condensed_fig8_components.append( (fig[0], fig[1], 1) )
        for fig in self.fig8_components[1:]:
            if fig == (condensed_fig8_components[-1][0], condensed_fig8_components[-1][1]):
                condensed_fig8_components[-1] = (fig[0], fig[1], condensed_fig8_components[-1][2] + 1)
            else:
                condensed_fig8_components.append( (fig[0], fig[1], 1) )
        
        
        result = str(self.nontrivial_component)

        for comp in self.other_components:
            print 'other component! ', self.other_components
            first_generator = comp[0]
            delta = self.alexander_gradings[first_generator] - self.maslov_gradings[first_generator]
            result += ', (' + str([self.alexander_gradings[item] for item in comp]) + ','+ str(delta) + ')'

        for fig in condensed_fig8_components:
            result += ', (' + str(fig[0]) +','+  str(fig[1]) + ')'
            if fig[2] > 1:
                result += '^'+str(fig[2])

        return result


#### The following section of code enables a graphical output of a KnotTrainTrack, in the form of an SVG file
    def first_crossing_left(self, i, j, max_depth='max'):
        '''returns the number (>=0) of other generator pairs passed leaving generators i and j to the left
        before the strands cross, or returns None if they do not cross'''

        if max_depth == 0:
            return None
        if max_depth == 'max':
            max_depth = self.num_generators
        if self.alexander_gradings[i] != self.alexander_gradings[j]:
            return None
        h_i = self.heights[i]
        h_j = self.heights[j]

        new_i = self.left_match_dict[i]
        if new_i == None:
            new_h_i = 2 * self.max_alex + 2  # stand in for 'infty'... larger than any other height
            new_i_alex = None
        else:
            new_h_i = self.heights[new_i]
            new_i_alex = self.alexander_gradings[new_i]
        new_j = self.left_match_dict[j]
        if new_j == None:
            new_h_j = 2 * self.max_alex + 2  # stand in for 'infty'... larger than any other height
            new_j_alex = None
        else:
            new_h_j = self.heights[new_j]
            new_j_alex = self.alexander_gradings[new_j]

        if (h_i - h_j) * (h_i - new_h_j) * (new_h_i - h_j) * (new_h_i - new_h_j) < 0:
            return 0  # here the two strands cross immediately
        elif new_i_alex != new_j_alex:
            return None
        else:
            recursive_result = self.first_crossing_right(new_i, new_j, max_depth - 1)
            if recursive_result == None:
                return None
            else:
                return recursive_result + 1

    def first_crossing_right(self, i, j, max_depth='max'):
        '''returns the number (>=0) of other generator pairs passed leaving generators i and j to the left
        before the strands cross, or returns None if they do not cross'''
        if max_depth == 0:
            return None
        if max_depth == 'max':
            max_depth = self.num_generators
        if self.alexander_gradings[i] != self.alexander_gradings[j]:
            return None
        h_i = self.heights[i]
        h_j = self.heights[j]

        new_i = self.right_match_dict[i]
        if new_i == None:
            new_h_i = 2 * self.max_alex + 2  # stand in for 'infty'... larger than any other height
            new_i_alex = None
        else:
            new_h_i = self.heights[new_i]
            new_i_alex = self.alexander_gradings[new_i]
        new_j = self.right_match_dict[j]
        if new_j == None:
            new_h_j = 2 * self.max_alex + 2  # stand in for 'infty'... larger than any other height
            new_j_alex = None
        else:
            new_h_j = self.heights[new_j]
            new_j_alex = self.alexander_gradings[new_j]

        if (h_i - h_j) * (h_i - new_h_j) * (new_h_i - h_j) * (new_h_i - new_h_j) < 0:
            return 0  # here the two strands cross immediately
        elif new_i_alex != new_j_alex:
            return None
        else:
            recursive_result = self.first_crossing_left(new_i, new_j, max_depth - 1)
            if recursive_result == None:
                return None
            else:
                return recursive_result + 1

    def is_in_bigon(self, i, j):
        return not (self.first_crossing_left(i, j) == None or self.first_crossing_right(i, j) == None)

    def set_heights(self):
        self.heights = [0 for item in range(self.num_generators)]

        for a in self.generators_by_grading.keys():
            num = len(self.generators_by_grading[a])
            shift = .8 / (num + 1)
            h = a + shift * (num - 1) / 2.0
            for g in self.generators_by_grading[a]:
                self.heights[g] = h
                h -= shift

        done = False
        pairs = []
        for alex in range(self.min_alex, self.max_alex + 1):
            for i in self.generators_by_grading[alex]:
                for j in self.generators_by_grading[alex]:
                    if j > i:
                        pairs.append((i, j))
        while not done:
            done = True
            for (i, j) in pairs:
                if self.is_in_bigon(i, j):
                    done = False
                    left_length = self.first_crossing_left(i, j)
                    right_length = self.first_crossing_right(i, j)
                    (self.heights[i], self.heights[j]) = (self.heights[j], self.heights[i])

                    (next_i, next_j) = (i, j)
                    counter = 0
                    for n in range(left_length):
                        if counter % 2 == 0:
                            next_i = self.left_match_dict[next_i]
                            next_j = self.left_match_dict[next_j]
                        else:
                            next_i = self.right_match_dict[next_i]
                            next_j = self.right_match_dict[next_j]
                        (self.heights[next_i], self.heights[next_j]) = (self.heights[next_j], self.heights[next_i])
                        counter += 1

                    (next_i, next_j) = (i, j)
                    counter = 0
                    for n in range(right_length):
                        if counter % 2 == 0:
                            next_i = self.right_match_dict[next_i]
                            next_j = self.right_match_dict[next_j]
                        else:
                            next_i = self.left_match_dict[next_i]
                            next_j = self.left_match_dict[next_j]
                        (self.heights[next_i], self.heights[next_j]) = (self.heights[next_j], self.heights[next_i])
                        counter += 1

    def make_svg(self, scale=500):
        gen_width = 1  # length of the horizontal segments representing generators
        show_weights = False
        show_maslov = True

        x_offset = scale * self.max_alex
        y_offset = scale * (self.max_alex + 1)
        dwg = svgwrite.Drawing('SVG_files/' + self.name + '.svg', profile='full',
                               size=(str(2 * x_offset) + "px", str(2 * y_offset) + "px"))

        self.set_heights()
        heights = self.heights

        # draw green dot for each puncture, which occur at half integer heights
        for i in range(self.min_alex, self.max_alex):
            dwg.add(dwg.circle((0 + x_offset, -(i + .5) * scale + y_offset), 3, stroke='green'))

        # draw a segment from
        for i in range(self.num_generators):
            dwg.add(dwg.line((-gen_width * scale / 2.0 + x_offset, -heights[i] * scale + y_offset),
                             (gen_width * scale / 2.0 + x_offset, -heights[i] * scale + y_offset), stroke='black'))
            if show_maslov:
                dwg.add(dwg.text(str(self.maslov_gradings[i]), insert=(x_offset, -heights[i] * scale + y_offset),
                                 stroke='red', font_size=.1 * scale))

        # draw left arcs
        for (i, j) in self.left_matching:
            h1 = heights[i]
            h2 = heights[j]
            draw_left_arc(dwg, h1, h2, gen_width / 2.0, scale, x_offset, y_offset)

        h = heights[self.left_unmatched]
        dwg.add(dwg.line((-gen_width * scale / 2.0 + x_offset, -h * scale + y_offset),
                         (-self.max_alex * scale + x_offset, -h * scale + y_offset), stroke='black'))

        # draw right arcs
        for (i, j) in self.right_matching:
            h1 = heights[i]
            h2 = heights[j]
            draw_right_arc(dwg, h1, h2, gen_width / 2.0, scale, x_offset, y_offset)

        h = heights[self.right_unmatched]
        dwg.add(dwg.line((gen_width * scale / 2.0 + x_offset, -h * scale + y_offset),
                         (self.max_alex * scale + x_offset, -h * scale + y_offset), stroke='black'))

        # draw arrows
        for a in self.arrows_by_grading.keys():
            num = len(self.arrows_by_grading[a])
            if num > 0:
                shift = .8 / num
                x_pos = -shift * (num - 1) / 2.0
                for arrow in self.arrows_by_grading[a]:
                    arrow.x_position = x_pos
                    x_pos += shift
                    draw_arrow(dwg, heights[arrow.init_gen], heights[arrow.final_gen], arrow.x_position, scale,
                               x_offset, y_offset)
                    if show_weights:  # display arrow weights
                        left_weight = self.left_weight(arrow.init_gen, arrow.final_gen, self.num_generators)
                        right_weight = self.right_weight(arrow.init_gen, arrow.final_gen, self.num_generators)

                        dwg.add(dwg.text(str(left_weight), insert=(
                        (arrow.x_position - .15) * scale + x_offset, (-heights[arrow.init_gen]) * scale + y_offset),
                                         stroke='red', font_size=.1 * scale))
                        dwg.add(dwg.text(str(right_weight), insert=(
                        (arrow.x_position + .05) * scale + x_offset, (-heights[arrow.init_gen]) * scale + y_offset),
                                         stroke='red', font_size=.1 * scale))

        dwg.save()

def make_svg_partial(generators, arrows, counter, hscale=60, vscale=40):
    show_weights = True
    h_dim = hscale * (len(arrows) + 2)
    v_dim = vscale * (len(generators) + 1)

    dwg = svgwrite.Drawing('SVG_files/debug' + str(counter) + '.svg', profile='full',
                           size=(str(h_dim) + "px", str(v_dim) + "px"))

    num_generators = len(generators)
    heights = [i * vscale for i in range(num_generators)]

    # draw horizontal segments
    for i in range(num_generators):
        dwg.add(dwg.line((0, heights[i]), (h_dim, heights[i]), stroke='black'))

    # draw arrows
    num_arrows = len(arrows)
    for i in range(num_arrows):
        arrow = arrows[i]
        x_pos = hscale * (i + 1)
        init_height = heights[generators.index(arrow.init_gen)]
        final_height = heights[generators.index(arrow.final_gen)]
        dwg.add(dwg.line((x_pos, init_height), (x_pos, final_height), stroke='red'))
        if init_height < final_height:
            dwg.add(dwg.line((x_pos - .2 * vscale, final_height - .2 * vscale), (x_pos, final_height), stroke='red'))
            dwg.add(dwg.line((x_pos + .2 * vscale, final_height - .2 * vscale), (x_pos, final_height), stroke='red'))
        if init_height > final_height:
            dwg.add(dwg.line((x_pos - .2 * vscale, final_height + .2 * vscale), (x_pos, final_height), stroke='red'))
            dwg.add(dwg.line((x_pos + .2 * vscale, final_height + .2 * vscale), (x_pos, final_height), stroke='red'))

        if show_weights:  # display arrow weights
            left_weight = arrow.left_weight()
            right_weight = arrow.right_weight()

            dwg.add(dwg.text(str(left_weight), insert=((x_pos - 15), (init_height + final_height) / 2), stroke='red',
                             font_size=12))
            dwg.add(dwg.text(str(right_weight), insert=((x_pos + 3), (init_height + final_height) / 2), stroke='red',
                             font_size=12))

    dwg.save()
##### end of code for SVG graphical output



def train_track_data_sort(alexander_gradings, maslov_gradings, left_matching, right_matching, arrows):
    '''rearranges data so that alexander_gradings is sorted in decreasing order'''
    generators = range(len(alexander_gradings))
    generators = sorted(generators, key = lambda x: alexander_gradings[x], reverse = True)
    new_alex = [alexander_gradings[g] for g in generators]
    new_maslov = [maslov_gradings[g] for g in generators]
    new_left_matching = []
    for (i,j) in left_matching:
        new_left_matching.append( (generators.index(i), generators.index(j)) )
    new_right_matching = []
    for (i, j) in right_matching:
        new_right_matching.append((generators.index(i), generators.index(j)))
    new_arrows = []
    for (i, j) in arrows:
        new_arrows.append((generators.index(i), generators.index(j)))
    return (new_alex, new_maslov, new_left_matching, new_right_matching, new_arrows)


def train_track_from_curve(curve, delta_shift):
    alexander_gradings = curve
    maslov_gradings = [curve[0] - delta_shift]
    left_matching = []
    right_matching = []
    arrows = []

    for i in range(len(curve)-1):
        if i%2 == 0:
            right_matching.append( (i, i+1) )
            deltaA = curve[i+1]-curve[i]
            deltaM = 2*deltaA - deltaA/abs(deltaA)
            maslov_gradings.append( maslov_gradings[i] + deltaM)
        else:
            left_matching.append( (i, i+1) )
            deltaA = curve[i + 1] - curve[i]
            deltaM = deltaA / abs(deltaA)
            maslov_gradings.append(maslov_gradings[i] + deltaM)

    (alex, maslov, left, right, arrows) = train_track_data_sort(alexander_gradings, maslov_gradings, left_matching, right_matching, arrows)
    return KnotTrainTrack(alex, maslov, left, right, arrows)

def tensor_product(track1, track2):
    '''assumes no arrows in input tracks'''
    new_generators = []

    new_alex = []
    new_maslov = []
    for g1 in range(track1.num_generators):
        for g2 in range(track2.num_generators):
            new_generators.append( (g1, g2) )
            new_alex.append( track1.alexander_gradings[g1]+track2.alexander_gradings[g2])
            new_maslov.append(track1.maslov_gradings[g1] + track2.maslov_gradings[g2])

    new_arrows = []

    new_right_matching = []
    for (i,j) in track1.right_matching:
        new_match_i = (i, track2.right_unmatched)
        new_match_j = (j, track2.right_unmatched)
        new_right_matching.append( (new_generators.index(new_match_i), new_generators.index(new_match_j)) )

    for (i, j) in track2.right_matching:
        if track2.alexander_gradings[i] < track2.alexander_gradings[j]:
            (i,j) = (j,i)
        new_match_i = (track1.right_unmatched, i)
        new_match_j = (track1.right_unmatched, j)
        new_right_matching.append((new_generators.index(new_match_i), new_generators.index(new_match_j)))
        for (k,l) in track1.right_matching:
            if track1.alexander_gradings[k] < track1.alexander_gradings[l]:
                (k, l) = (l, k)

                new_match_a = (k,i)
                if track1.alexander_gradings[k] + track2.alexander_gradings[j] > track1.alexander_gradings[l] + track2.alexander_gradings[i]:
                    new_match_b = (k,j)
                    new_match_c = (l,i)
                else:
                    new_match_b = (l,i)
                    new_match_c = (k,j)
                new_match_d = (l,j)

                new_right_matching.append((new_generators.index(new_match_a), new_generators.index(new_match_b)))
                new_right_matching.append((new_generators.index(new_match_c), new_generators.index(new_match_d)))

                if track1.alexander_gradings[k] + track2.alexander_gradings[j] == track1.alexander_gradings[l] + track2.alexander_gradings[i]:
                    new_arrows.append((new_generators.index(new_match_c), new_generators.index(new_match_b)))

    new_left_matching = []
    for (i, j) in track1.left_matching:
        new_match_i = (i, track2.left_unmatched)
        new_match_j = (j, track2.left_unmatched)
        new_left_matching.append((new_generators.index(new_match_i), new_generators.index(new_match_j)))

    for (i, j) in track2.left_matching:
        if track2.alexander_gradings[i] < track2.alexander_gradings[j]:
            (i, j) = (j, i)
        new_match_i = (track1.left_unmatched, i)
        new_match_j = (track1.left_unmatched, j)
        new_left_matching.append((new_generators.index(new_match_i), new_generators.index(new_match_j)))
        for (k, l) in track1.left_matching:
            if track1.alexander_gradings[k] < track1.alexander_gradings[l]:
                (k, l) = (l, k)

                new_match_a = (k, i)
                if track1.alexander_gradings[k] + track2.alexander_gradings[j] > track1.alexander_gradings[l] + track2.alexander_gradings[i]:
                    new_match_b = (k, j)
                    new_match_c = (l, i)
                else:
                    new_match_b = (l, i)
                    new_match_c = (k, j)
                new_match_d = (l, j)

                new_left_matching.append((new_generators.index(new_match_a), new_generators.index(new_match_b)))
                new_left_matching.append((new_generators.index(new_match_c), new_generators.index(new_match_d)))

                if track1.alexander_gradings[k] + track2.alexander_gradings[j] == track1.alexander_gradings[l] + track2.alexander_gradings[i]:
                    new_arrows = [(new_generators.index(new_match_b), new_generators.index(new_match_c))] + new_arrows


    (alex, maslov, left, right, arrows) = train_track_data_sort(new_alex, new_maslov, new_left_matching, new_right_matching, new_arrows)
    return KnotTrainTrack(alex, maslov, left, right, arrows)
