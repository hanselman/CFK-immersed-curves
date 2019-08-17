upper = 'abcdefghijklmnopqrstuvwxyz'
lower = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
letter_to_number = {}
for i in range(26):
    letter_to_number[lower[i]] = i+1
    letter_to_number[upper[i]] = -(i+1)
    
def str_to_DT(string):
    return [2*letter_to_number[char] for char in string]

def DTcrossing_to_PDcrossing(crossing):
    (a,b) = crossing
    is_negative = False
    if a < 0:
        is_negative = True
        a = abs(a)
    elif b < 0:
        is_negative = True
        b = abs(b)
    ## the edge labels at this crossing in PD notation are {a, a+1, b, b+1}

def planar_graph_step(vertices, vertex_count, vertex_or, regions, current_end):
    
    #reindex regions so that current_end is in regions[0]
    i = 0
    while not current_end in regions[i]:
        i += 1
    regions = regions[i:]+regions[:i]

    #if every vertex has been visited twice
    if vertex_count == 2*len(vertices):
        i = 0
        while not 1 in vertices[i]:
            i+=1
        next_vertex = vertices[i]

        if (next_vertex, 0) in regions[0]:
            return vertex_or
        else:
            return None
            

    #find the next vertex
    i = 0
    while not vertex_count + 1 in vertices[i]:
        i += 1
    next_vertex = vertices[i]

    #if the next vertex hasn't been seen yet, just insert it at current_end
    if min(next_vertex) == vertex_count+1: 
        i = regions[0].index( current_end )
        new_region0 = regions[0][:i] + [(next_vertex, 1), (next_vertex, 2), (next_vertex, 3)] + regions[0][i+1:]
        new_regions = [r for r in regions]
        new_regions[0] = new_region0
        next_end = (next_vertex, 2)
        
        return planar_graph_step(vertices, vertex_count+1, vertex_or + [0], new_regions, next_end)

    else:
        if (next_vertex, 1) in regions[0]:
            i = regions[0].index( current_end )

            region0 = regions[0][i:] + regions[0][:i]  #cyclically permute so that current end is first
            
            j = region0.index( (next_vertex, 1) )

            new_region0a = region0[1:j]
            new_region0b = region0[j+1:]
            new_regions = [new_region0a, new_region0b]
            new_regions += regions[1:]
            
            next_end = (next_vertex, 3)

            result = planar_graph_step(vertices, vertex_count+1, vertex_or + [1], new_regions, next_end)
            if not result == None:
                return result

        if (next_vertex, 3) in regions[0]:
            i = regions[0].index( current_end )

            region0 = regions[0][i:] + regions[0][:i]  #cyclically permute so that current end is first
            
            j = region0.index( (next_vertex, 3) )

            new_region0a = region0[1:j]
            new_region0b = region0[j+1:]
            new_regions = [new_region0a, new_region0b]
            new_regions += regions[1:]
            
            next_end = (next_vertex, 1)

            result = planar_graph_step(vertices, vertex_count+1, vertex_or + [3], new_regions, next_end)
            if not result == None:
                return result
        return None

def DT_to_PD(DTcode):
    '''input is list of even integers (with sign)'''

    n = len(DTcode)
    vertices = []
    for i in range(n):
        vertices.append( (2*i+1, abs(DTcode[i]) ) )
    DTsigns = {}
    for i in range(1,n+1):
        if 2*i in DTcode:
            DTsigns[2*i] = 2*i
        else:
            DTsigns[2*i] = -2*i

    # a quadrant of a vertex is a tuple (vertex, k) where k is 0,1,2, or 3.
    # These represent the 4 possible edges leaving the vertex. We can think of 0,1,2,3 as E,S,W,N

    # As we gradually build up the graph, we keep track of the regions of the complement of the graph in S^2
    #each region is labeled by a list of quadrants of vertices along the outside of the region,
    #in order given by the orientation of the boundary of that region
        
    starting_region = [(vertices[0], 0), (vertices[0],1), (vertices[0],2), (vertices[0],3)]
    regions = [starting_region]
    vertex_count = 1
    vertex_or = [0]
    current_end = (vertices[0],2)

    vertex_orientations = planar_graph_step(vertices, vertex_count, vertex_or, regions, current_end)
    ## vertex_orientations is now a length  integer (0,1,or 3)
    ## each time a vertex is passed, the integer is 0 the first time and 1 or 3 the second time

    #print vertex_orientations 
    PD_code = []
    for (a,b) in vertices:
        # a is odd, b is even
        # n is negative iff it is the overcrossing

        if DTsigns[b] > 0:
            if b > a:
                if vertex_orientations[b-1] == 1:
                    newPDvertex = [b,a,b+1,a+1]
                elif vertex_orientations[b-1] == 3:
                    newPDvertex = [b,a+1,b+1,a]
                else:
                    print 'error a', a, b
            elif b < a:
                if vertex_orientations[a-1] == 1:
                    newPDvertex = [b,a+1,b+1,a]
                elif vertex_orientations[a-1] == 3:
                    newPDvertex = [b,a,b+1,a+1]
                else:
                    print 'error b'
            else:
                print 'error c'
        elif DTsigns[b] < 0:
            if b > a:
                if vertex_orientations[b-1] == 1:
                    newPDvertex = [a,b+1,a+1,b]
                elif vertex_orientations[b-1] == 3:
                    newPDvertex = [a,b,a+1,b+1]
                else:
                    print 'error d'
            elif b < a:
                if vertex_orientations[a-1] == 1:
                    newPDvertex = [a,b,a+1,b+1]
                elif vertex_orientations[a-1] == 3:
                    newPDvertex = [a,b+1,a+1,b]
                else:
                    print 'error e'
            else:
                print 'error f'
        else:
            print 'error: ', DTsigns[b]
        for i in range(4):
            if newPDvertex[i] == 2*n+1:
                newPDvertex[i] = 1
        PD_code.append(newPDvertex)
                

    return PD_code




######## To run for other sets of knots, change filename and n (crossing number) here
###      filename points to a txt file with a DT code on each line, all knots must have the same crossing number n
filename = '12n'
n = 12
###



f = open(filename + 'DT.txt', 'r')
lines = f.readlines()
f.close()
DTcodes = [str_to_DT(item[0:n]) for item in lines]
PDcodes = [DT_to_PD(item) for item in DTcodes]

f = open(filename+'PD.txt', 'w')
for code in PDcodes:
    f.write('PD '+str(code) + '\n')
f.close()

