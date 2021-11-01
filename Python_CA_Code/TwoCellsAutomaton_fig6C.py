import numpy as np

import matplotlib.pyplot as plt

import time

import matplotlib.pyplot as plt

from matplotlib.colors import ListedColormap

from colorama import Fore

import os



#methods

# # Toroidal BC

# def copy_boundaries(array, range0):

    

#     #copy left and right

#     for j in range(range0, size + range0):

#         for i in range(0, range0):

#             array[i][j] = array[size + i][j]

#             array[size + range0 + i][j] = array[range0 + i][j]



#     #copy top and bottom

#     for k in range(0, size + 2*range0):

#         for i in range(0, range0):

#             array[k][i] = array[k][size + i]

#             array[k][snize + range0 + i] = array[k][range0 + i]



# concatenates the part of grid we care about into a new 3 x 3 grid

def extendArray(array, size):

    trunc_array = array[size: -size, size: -size]

    arrayx3 = np.concatenate((trunc_array, trunc_array, trunc_array),axis=1)

    array3x3 = np.concatenate((arrayx3, arrayx3, arrayx3),axis=0)

    return array3x3





#calculate maximum payoff

def calcMaxPayoff(mumax, range0, range1):

    x=0

    y=0

    #calc x=max payoff of type 1, y=max payoff of type 0

    if payoffType==1:

        x=1/km

        y=mumax/km



    if payoffType==2:

        x=((2*range1+1)*(2*range1+1)-1)/km

        y=mumax*((2*range0+1)*(2*range0+1)-1)/km



    if payoffType==3:

        x=((2*range1+1)*(2*range1+1)-1)/(km+((2*range1+1)*(2*range1+1)-1))

        y=mumax*((2*range0+1)*(2*range0+1)-1)/(km+((2*range0+1)*(2*range0+1)-1))



    return max(x,y)





def initialize_grid(initial_fraction):

    for i in range(size, 2*size):

        for j in range(size, 2*size):

            #Initialize with customized ratio of particle 0 to 1

            if np.random.rand() >= initial_fraction: 

                grid_init[i][j] = 0

            else: 

                grid_init[i][j] = 1

    

    return extendArray(grid_init, size)



# EDIT SVV: avoid using global variables by explictly passing grid As input

def compute_patch_size2(grid):

    # Count patch sizes horizontally, taking toroidal BC into account

    patch0ArrayHorizontal = np.zeros((size, size), dtype = int)

    patch1ArrayHorizontal = np.zeros((size, size), dtype = int)

    zeroCount = 0

    oneCount = 0

    for i in range(size, 2*size):

        for j in range(size, 2*size):

            if grid[i][j] == 0:

                zeroCount = zeroCount + 1

                patch1ArrayHorizontal[i-size][j-size]  = oneCount

                oneCount = 0

            elif grid[i][j] == 1:

                oneCount = oneCount + 1

                patch0ArrayHorizontal[i-size][j-size]  = zeroCount

                zeroCount = 0



    # Count patch sizes vertically, taking toroidal BC into account

    patch0ArrayVertical = np.zeros((size, size), dtype = int)

    patch1ArrayVertical = np.zeros((size, size), dtype = int)

    zeroCount = 0

    oneCount = 0

    for i in range(size, 2*size):

        for j in range(size, 2*size):

            if grid[j][i] == 0:

                zeroCount = zeroCount + 1

                patch1ArrayVertical[j-size][i-size] = oneCount

                oneCount = 0

            elif grid[j][i] == 1:

                oneCount = oneCount + 1

                patch0ArrayVertical[j-size][i-size] = zeroCount

                zeroCount = 0





    avg0_horizontal = patch0ArrayHorizontal[patch0ArrayHorizontal.nonzero()].mean()

    avg0_vertical = patch0ArrayVertical[patch0ArrayVertical.nonzero()].mean()

    avg1_horizontal = patch1ArrayHorizontal[patch1ArrayHorizontal.nonzero()].mean()

    avg1_vertical = patch1ArrayVertical[patch1ArrayVertical.nonzero()].mean()

    patchSize0 = avg0_horizontal * avg0_vertical

    patchSize1 = avg1_horizontal * avg1_vertical



    return patchSize0, patchSize1





def local_fracs(a0, a1, gridSum):

  

    local_frac0s = a0*(km/mumax)/(gridSum/(size*size))

    local_frac1s = a1*km/((size*size-gridSum)/(size*size))



    return local_frac0s, local_frac1s







def compute_payoffs(mumax, range0, range1):

  

    neighbors0 = (2*range0+1)*(2*range0+1)-1

    neighbors1 = (2*range1+1)*(2*range1+1)-1

    neighborsOnes = 0

    diff_neighb = 0

    

    for row in range(size, 2*size):

        for col in range(size, 2*size):

        #initialize count of ones to zero

            neighborsOnes = 0

            if grid[row][col] == 0:      

                #count neighbors that are ones

                for i in range(row - range0, row + range0 + 1):

                    for j in range(col - range0, col + range0 + 1):

                        neighborsOnes += grid[i][j]

                diff_neighb = neighborsOnes

                #calculate payoff

                if payoffType==1:

                    payoff[row][col] = mumax * diff_neighb / neighbors0 / km #average

                if payoffType==2:

                    payoff[row][col] = mumax * diff_neighb / km #sum

                if payoffType==3:

                    payoff[row][col] = mumax * diff_neighb / (neighbors0 + km) #monod



            if grid[row][col] == 1:      

                #count neighbors that are ones

                for i in range(row - range1, row + range1 + 1):

                    for j in range(col - range1, col + range1 + 1):

                        neighborsOnes += grid[i][j]

                #subtract element grid[row][col] = 1

                neighborsOnes -= grid[row][col]

                #number of 0-type neighbors of element (row,col)

                diff_neighb  = neighbors1 - neighborsOnes

                #calculate payoff

                if payoffType==1:

                    payoff[row][col] = 1 * diff_neighb / neighbors1 / km #average

                if payoffType==2:

                    payoff[row][col] = 1 * diff_neighb / km #sum

                if payoffType==3:

                    payoff[row][col] = 1 * diff_neighb / (neighbors1 + km) #monod



    return extendArray(payoff, size)





# only compute payoffs that have changed (important for program speed) #EDIT SVV add grid and payoff as inputs

def compute_changed_payoffs(grid, payoff, mumax, range0, range1, actionRow, actionCol, previous_val, new_val):

  

    if previous_val != new_val:



        neighbors0 = (2*range0+1)*(2*range0+1)-1

        neighbors1 = (2*range1+1)*(2*range1+1)-1

        neighborsOnes = 0

        diff_neighb = 0

        

        for row in range(actionRow - range0, actionRow + range0 + 1):

            for col in range(actionCol - range0, actionCol + range0 + 1):



                if grid[row][col]==0 and abs(row-actionRow)<=range0 and abs(col-actionCol)<=range0 and (row != actionRow or col != actionCol):

                    if payoffType==1:

                        payoff[row][col] = payoff[row][col] - (-1) ** new_val * mumax / neighbors0 / km #average

                    if payoffType==2:

                        payoff[row][col] = payoff[row][col] - (-1) ** new_val * mumax / km #sum

                    if payoffType==3:

                        payoff[row][col] = payoff[row][col] - (-1) ** new_val * mumax / (neighbors0 + km) #monod



                if grid[row][col]==1 and abs(row-actionRow)<=range1 and abs(col-actionCol)<=range1 and (row != actionRow or col != actionCol):

                    if payoffType==1:

                        payoff[row][col] = payoff[row][col] + (-1) ** new_val * 1 / neighbors1 / km #average

                    if payoffType==2:

                        payoff[row][col] = payoff[row][col] + (-1) ** new_val * 1 / km #sum

                    if payoffType==3:

                        payoff[row][col] = payoff[row][col] + (-1) ** new_val * 1 / (neighbors1 + km) #monod



                if grid[row][col]==0 and row == actionRow and col == actionCol:

                    for i in range(row - range0, row + range0 + 1):

                        for j in range(col - range0, col + range0 + 1):

                            neighborsOnes += grid[i][j]

                    #how many equal neighbors element(row col) has

                    diff_neighb = neighborsOnes

                    #calculate payoff

                    if payoffType==1:

                        payoff[row][col] = mumax * diff_neighb / neighbors0 / km #average

                    if payoffType==2:

                        payoff[row][col] = mumax * diff_neighb / km #sum

                    if payoffType==3:

                        payoff[row][col] = mumax * diff_neighb / (neighbors0 + km) #monod



                if grid[row][col]==1 and row == actionRow and col == actionCol:

                    for i in range(row - range1, row + range1 + 1):

                        for j in range(col - range1, col + range1 + 1):

                            neighborsOnes += grid[i][j]

                    #subtract element(row,col)

                    neighborsOnes -= 1

                    #how many equal neighbors element(row col) has

                    diff_neighb  = neighbors1 - neighborsOnes

                    #calculate payoff

                    if payoffType==1:

                        payoff[row][col] = 1 * diff_neighb / neighbors1 / km #average

                    if payoffType==2:

                        payoff[row][col] = 1 * diff_neighb / km #sum

                    if payoffType==3:

                        payoff[row][col] = 1 * diff_neighb / (neighbors1 + km) #monod



        #copy payoff to boundaries

        return extendArray(payoff, size)

    else:

        return payoff

#end of compute_changed_payoffs



# EDIT SVV: avoid using global variables by explictly passing grid and payoff as inputs

def truncate_arrays(grid, payoff):

    trunc_grid = grid[size: -size, size: -size]

    trunc_payoff = payoff[size: -size, size: -size]

    return trunc_grid, trunc_payoff

 

#avg payoffs of 0-type cells

def avg_payoffs0(trunc_grid, trunc_payoff):

    inv_trunc_grid = 1-trunc_grid

    return np.trace(np.inner(trunc_payoff, inv_trunc_grid)) / np.sum(inv_trunc_grid)



#avg payoffs of 1-type cells

def avg_payoffs1(trunc_grid, trunc_payoff):

    return np.trace(np.inner(trunc_payoff, trunc_grid)) / np.sum(trunc_grid)



def sum_grid(trunc_grid):

    return np.sum(trunc_grid)



#one particle reproduces and kills a random neighbour

def replace_and_kill(grid, payoff):



    while True:

        #this section finds who will reproduce

        aRand = int(np.random.rand() * (size * size))

        #find corresponding individual. Note that int() truncates the decimal, hence row < rangeMax + size

        row = int(aRand/size)+size

        col = int(aRand%size)+size



        #decide whether individual will reproduce or search another

        aRand2 = np.random.rand()

        prob_aus=payoff[row][col]/maxPayoff



        if(aRand2 < prob_aus):

            break



    #calculate range where to pick neighbors to be replaced

    range_replace = (1-replaceAdjacentNeighbor) * (range0*(1-grid[row][col]) + range1*grid[row][col]) + 1 * replaceAdjacentNeighbor

   

    # this loop prevents cell from reproducing and replacing itself

    while True:

        #select random number between 0 and totalNeighbors-1

        aRand3 = int(np.random.rand() * (2*range_replace+1) * (2*range_replace+1))

        #find corresponding individual

        relRow = int(aRand3/(2*range_replace+1))

        relCol = int(aRand3%(2*range_replace+1))

        if relRow != range_replace or relCol != range_replace:

            break



    actionRow = row-range_replace+relRow

    actionCol = col-range_replace+relCol

    previous_val = grid[actionRow][actionCol]

    #replace 

    grid[actionRow][actionCol]=grid[row][col]

    new_val = grid[actionRow][actionCol]



    grid = extendArray(grid, size)

    # EDIT SVV: output grid explicitly: not clear if update in place works with extendArray function

    return grid, actionRow, actionCol, previous_val, new_val





def plot_grid(trunc_grid, path=None, title=None, pathLength = None, lineColor='red'):

    fig, ax = plt.subplots()

    

    colors = 'white red'.split()

    cmap = ListedColormap(colors, name='colors', N=None)

    

    img = ax.matshow(trunc_grid, zorder=1, cmap=cmap)

    

    size = len(trunc_grid)

    labels = np.arange(0,size,5) # controls tick label spacing

    ax.spines['left'].set_position(('axes', 0))

    ax.spines['bottom'].set_position(('axes', 0))

    # ax.spines['right'].set_color('none')

    ax.invert_yaxis()

    # fig.colorbar(img)

    ax.xaxis.set_ticks_position('bottom')

    ax.yaxis.set_ticks_position('left')

    

    # Shift ticks to be at 0.5, 1.5, etc

    locs = np.arange(0,size,5)

    locs2 =np.arange(0,size)

    for axis in [ax.xaxis, ax.yaxis]:

        axis.set_ticks(locs2 + 0.5, minor=True) # controls gridline spacing via locs2

        axis.set(ticks=locs, ticklabels=labels) # controls tickmark spacing via locs

    # Turn on the grid for the minor ticks

    ax.grid(True, which='minor')

    

    if path!=None:

        x_vals = [x[0] for x in path]

        y_vals = [x[1] for x in path]

        ax.plot(x_vals, y_vals, zorder=2, Color=lineColor)

        start = (x_vals[0], y_vals[0])

        ax.annotate('start', xy = start,Color=lineColor)

        if pathLength!=None:

            end = (x_vals[pathLength-1], y_vals[pathLength-1])

            ax.annotate('end', xy = (10,10),Color=lineColor)

        

    if title!=None:

        ax.set_title(title)

        

    if pathLength!=None:

        txt = ('Path Length = '+str(pathLength))

        plt.figtext(0.5, 0.01, txt, wrap=True, horizontalalignment='center', fontsize=16)

        

    name = int(t + 100*r + s)

    fig.savefig('F6C{}.png'.format(name))

    return



# this method for printing payoff in terminal works, but not very useful        

def print_payoff(range0):

    sizeTot = size + 2*range0

    for i in range(sizeTot):

        for j in range(sizeTot):

            if grid[i][j] == 1:

                print(Fore.RED + str(payoff[i][j]), end=' ')

            else:

                print(Fore.BLUE + str(payoff[i][j]), end=' ')

        print()




def saveArrays(fracA_in_time):
    file = open("Fig6C.m", "a")

    for i in range(num_initial_fractions):
        file.write("fracA_in_time({},:) = [".format(i+1))
        for j in range(int(timeSteps/save_freq)):
            file.write(" {:.5f}".format(fracA_in_time[i][j]))
        if j < int(timeSteps/save_freq) - 1:
            file.write(";")
        else:
            file.write("];\n")


    file.write("sizeGrid = {};\n".format(size))
    file.write("range0 = {};\n".format(range0))
    file.write("range1 = {};\n".format(range1))
    file.write("mumax = {}; ... mumax=mu0/mu1\n".format(mumax))
    file.write("payoffType = {};\n".format(payoffType))
    file.write("replaceAdjacentNeighbor = {};\n".format(int(replaceAdjacentNeighbor)))
    file.write("km = {};\n".format(km))
    file.write("initial_fractions = {};\n".format(initial_fractions))
    file.write("timeSteps = {};\n".format(timeSteps))
    file.write("runs = {};\n".format(runs))
    file.write("runtime = {};\n".format(time.time() - t0))

    
    file.close()

# end of methods

############################################# begin main "method" ##########################################



t0 = time.time()



# initialize variables

size = 100

km = 49

replaceAdjacentNeighbor = True

timeSteps = 100000

runs = 1

PRINT = False

SAVE = True

payoffType = 1

num_initial_fractions = 22

save_freq = 10




grid_init = np.zeros((3*size, 3*size), dtype = int)

payoff = np.zeros((3*size, 3*size))

initial_fractions = np.linspace(1/24, 23/24, num=num_initial_fractions, endpoint=True)

# initial_fractions = np.linspace(3/4, 1, num=num_initial_fractions, endpoint=False)


# remove file with previous name, if it exists

try:

    os.remove(r"C:\Users\kofri\CellularAutomaton1\Fig6C.m")

except OSError:

    pass



# preallocate memory for payoff and fraction in time

# sums = np.zeros((num_initial_fractions, tail))

# payoffs0 = np.zeros((num_initial_fractions, tail))

# payoffs1 = np.zeros((num_initial_fractions, tail))

# patch0 = np.zeros((num_initial_fractions, tail))

# patch1 = np.zeros((num_initial_fractions, tail))

# local_fraction0s = np.zeros((num_initial_fractions,tail))

# local_fraction1s = np.zeros((num_initial_fractions,tail))

fracA_in_time = np.zeros((num_initial_fractions, int(timeSteps/save_freq)))



## BEGIN COMPUTING AUTOMATA

for s in range(num_initial_fractions):

    initial_fraction = initial_fractions[s]

    mumax = 3.6

    range0 = 5

    range1 = 1

    maxPayoff = calcMaxPayoff(mumax, range0, range1)


    # initialize grid and calculate its initial payoff

    grid = initialize_grid(initial_fraction)

    # payoff = compute_payoffs(mumax, range0, range1)

    payoff = compute_payoffs(mumax, range0, range1)



    for t in range(timeSteps):



        # kill one particle and replace

        (grid, actionRow, actionCol, previous_val, new_val) = replace_and_kill(grid, payoff)

        # replace_and_kill(grid, payoff)

        # compute new payoffs



        #EDIT SVV add payoff as explicit input

        payoff = compute_changed_payoffs(grid, payoff, mumax, range0, range1, actionRow, actionCol, previous_val, new_val)

        # payoff = compute_payoffs(mumax, range0, range1)

        

        if PRINT == True and t == timeSteps - 1:

            trunc_grid, trunc_payoff = truncate_arrays(grid, payoff)

            plot_grid(trunc_grid)



        if SAVE == True and t % save_freq == 0:

            trunc_grid, trunc_payoff = truncate_arrays(grid, payoff)

            gridSum = sum_grid(trunc_grid)



            # patchSize0, patchSize1 = compute_patch_size2(grid)

            # local_fraction0s[s][j], local_fraction1s[s][j] = local_fracs(a0, a1, gridSum)

            fracA_in_time[s][int(t/save_freq)] = gridSum/size/size

            # sums[s][j] = gridSum

            # payoffs0[s][j] = a0

            # payoffs1[s][j] = a1

            # patch0[s][j] = patchSize0

            # patch1[s][j] = patchSize1





if SAVE == 1:

    saveArrays(fracA_in_time)

            

print("runtime = ", time.time() - t0, "s")



############################################# end main "method" #############################################