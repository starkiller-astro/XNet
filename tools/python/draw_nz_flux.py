def draw_nz_mf(datafile, k, num_species, x_limit='None', colors):

    '''
        1. Scatter plot color coded by mass fractions at a specified timestep.
        2. Draw vectors for the reaction flux on the NZ plane from reaction target
           to product colored by relative flux intensity.
        Inputs: datafile = must be a ts file
                k = desired time for thermodynamic and mass fraction data
                num_species = number of plotted species desired
                colors = turns colorbar on or off (mass fraction coloring)
        Outputs: scatterplot of number of neutrons vs number of protons, color coded for mass fraction
                 thermodynamic data at specific time

    '''
    import read_ts_file as rtf
    import numpy as numpy
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.ticker as plticker
    import matplotlib.colors
    import draw_nz_background as back

#-------------------------------------------------------------------------------
#   MASS FRACTION SCATTERPLOT

    #Create plot space.
    fig = plt.figure(1)

    #Read ts file and use variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile)

    #Plot proton number vs neutron number and color code for mass fraction (log).
    for counter in numpy.arange(0, num_species):
        nn = aa[counter] - zz[counter]
        nn = int(round(nn))
        if colors == True:
            s = plt.scatter(nn, zz[counter], s = 25, c = xmf[k, counter], cmap = 'YlOrBr', norm = matplotlib.colors.LogNorm(vmin=xmf[k].min(), vmax=xmf[k].max()), marker = 's')
        else:
            s = plt.scatter(nn, zz[counter], s = 5, c = 'w')

    #Make nn back to an array
    nn = aa - zz

    #Plot colorbar as legend.
    if colors == True:
        cb = plt.colorbar(s)
        cb.set_label('Log of Mass Fraction')

    #Print thermodynamic data.
    x_loc = 5
    y_loc = 37
    words = ('Time = %s sec') % time[k]
    plt.text(x_loc, y_loc, words)
    words = ('Temperature = %s GK') % temperature[k]
    plt.text(x_loc, y_loc - 1.5, words)
    words = (r'$\rho =$ %s g/$\mathregular{cm^3}$') % density[k]
    plt.text(x_loc, y_loc - 3, words)
    words = ('Timestep = %s') % k
    plt.text(x_loc, y_loc - 4.5, words)

    #Format x and y axes.
    plt.ylim(-1, 40)
    plt.ylabel('Z (Proton Number)')
    if x_limit == 'None':
            x_limit = 40
    plt.xlim(-2, x_limit)
    plt.xlabel('N (Neutron Number)')

    #Set axis ticks
    special_elements = [2, 8, 15, 28]
    loc = plticker.FixedLocator(special_elements)
    plt.xticks(special_elements, special_elements)
    plt.yticks(special_elements, special_elements)

    #Label elements and mass numbers
    z_max = max(zz)
    z_min = min(zz)
    n_max = max(nn)
    n_min = min(nn)
    leftn  = numpy.empty([0,])
    rightn = numpy.empty([0,])
    alabel = numpy.empty([0,])

    #Loop over all elements, including neutrons (z=0)
    z_max = int(z_max) + 1   #33
    for i in range(z_max):   #[0,1,2,...,31,32]
        zint = numpy.arange(z_max)   #[0 1 2 ... 31 32]
        zint[i] = i
        iz = numpy.where(zz == i)[0]

        #Find the lightest and heaviest isotope
        if not iz.size==0:
        # returns True if the condition is considered False
        # if proton no. array is not empty
            nn = nn.astype(int)
            leftn = numpy.append(leftn, min(nn[iz]))
            leftn = leftn.astype(int)
            rightn = numpy.append(rightn, max(nn[iz]))
            rightn = rightn.astype(int)
            alabel = numpy.append(alabel,rightn[i] + zint[i])
            alabel = alabel.astype(int)
            # max neutron no. + i in range(z_max)

    #Load element names
    element = rtf.build_element_symbol()
    element = numpy.asarray(element,dtype=object)

    for i in range(0,113):
        if i in zint:
            #Label elements to the left
            plt.text(leftn[i]-0.75, zint[i], element[i], fontsize=9, horizontalalignment='left')

            #Label mass numbers to the right
            plt.text(rightn[i]+1, zint[i]-0.5, alabel[i], fontsize=9, horizontalalignment='right', rotation=-45)

#-------------------------------------------------------------------------------
#   FLUX VECTORS

    #Identify starting and ending points for all fluxes
    zt = zz[flx_end[:,0]-1]
    zp = zz[flx_end[:,1]-1]
    nt = nn[flx_end[:,0]-1]
    np = nn[flx_end[:,1]-1]

    #Choose cycle(timestep) to plot
    flx_write = flx[k,:]

    #Choose maximum flux for calculating bin ranges
    flx_max = max(abs(flx_write))

    #Calculate vector origin and lengths
    aflx = abs(flx_write)
    zo = zt # vector origin
    no = nt # vector origin
    zv = zp-zt # vector length (z-component)
    nv = np-nt # vector length (n-component)

    #Reverse arrows for negative flux
    ineg = numpy.where(flx_write < 0.0)[0]
    zo[ineg] = zp[ineg]
    no[ineg] = np[ineg]
    zv[ineg] = -zv[ineg]
    nv[ineg] = -nv[ineg]

    #set number of colors/levels
    color_bins = numpy.array([1e-1, 1e-2, 1e-3, 1e-4, 1e-6])

    #set number of arrays
    colors = numpy.array(['k','r','g','b','c'])
    ncolors = len(color_bins)
    bin_upper = numpy.array([1])
    bin_upper = numpy.append(bin_upper,color_bins[0:4])
    bin_lower = color_bins

    #Build Legend array
    legend_array = numpy.array(range(ncolors), dtype=object)

    #Reverse ncolors so that for loop runs through color_bins backwards
    ncolors = range(ncolors)
    rev_ncolors = ncolors.reverse()
    ncolors = numpy.asarray(ncolors)
    # ultimately, heavier flux vectors will be plotted on top because
    # they are calculated later than lighter fluxes

    #Loop over colors
    for i in ncolors:

        #Plot Vectors
        flx_top = flx_max * bin_upper[i] #scalar multiplication
        flx_bot = flx_max * bin_lower[i]
        ii = numpy.where((aflx <= flx_top) & (aflx > flx_bot))[0]
        scale = 0.0 #prevents autoscaling
        if i == 0:
            lwidth= 0.002
        elif i == 1:
            lwidth= 0.0015
        else:
            lwidth= 0.001

        #Assemble handles of legend
        bin_lower[i] = str(bin_lower[i])
        label = [' > ','{:.1e}'.format(bin_lower[i]),' of max'] #prints bin_lower[i] in scientific notation
        label = ''.join(map(str,label)) #joins string elements of a sequence(i.e. label is a list)
        legend_array[i] = label

        #Draw vectors
        h = fig.add_subplot(111)
        h.quiver(no[ii],zo[ii],nv[ii],zv[ii],width=lwidth,color=colors[i],angles='xy',scale_units='xy',scale=1,label=legend_array[i])
        plt.grid()
        plt.draw()


    #Plot legend
    plt.legend(loc=4,title='Reaction Fluxes')
    # automatically pulls in label from quiver plot
    # and label=legend_array[i]

#-------------------------------------------------------------------------------
#   FINALIZE PLOT

    #Add grid lines.
    plt.grid(True)

    #Add title.
    plt.title(datafile)

    #Show plot.
    plt.show()
