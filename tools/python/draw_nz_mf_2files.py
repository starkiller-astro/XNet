def draw_nz_mf(datafile, datafile2, k, num_species, colors = False):

    '''
        Colorbar is currently not working. 
        Inputs: datafile = must be a ts file
                datafile2 = second datafile
                k = desired timestep for thermodynamic and mass fraction data
                num_species = number of plotted species desired
                colors = turn colorbar on or off, defaults to False
        Outputs: scatterplot of number of neutrons vs number of protons, color coded for mass fraction
                 thermodynamic data at specific time
    '''
    import read_ts_file as rtf
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    import matplotlib.ticker as plticker
    import matplotlib.colors

    fig = plt.figure(1)

    #Read ts file, use variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile)

    #For each species, determine neutron number by subtracting proton number from mass number, and plot neutron number by proton number. Apply a colormap, and normalize the data on a log scale.
    ax1 = plt.subplot2grid((1, 2), (0, 0))
    for counter in np.arange(0, num_species):
        nn = aa[counter] - zz[counter]
        nn = int(round(nn))
        if colors:
            s = ax1.scatter(nn, zz[counter], s = 60, c = xmf[k, counter], cmap = 'YlGnBu_r', norm = matplotlib.colors.LogNorm(vmin=xmf[k].min(), vmax=xmf[k].max()), marker = 's')
        else:
            s = ax1.scatter(nn, zz[counter], s = 60, marker = 's')

    #Configure and label axes, placing ticks at helium, oxygen, calcium, and nickel.
    plt.ylim(-1, 40)
    plt.ylabel('Z (Proton Number)')
    plt.xlim(-2, 40)
    plt.xlabel('N (Neutron Number)')

    special_elements = [2, 8, 20, 28]
    loc = plticker.FixedLocator(special_elements)
    plt.xticks(special_elements, special_elements)
    plt.yticks(special_elements, special_elements)

    #Add labels for each element.
    element = rtf.build_element_symbol()
    for item in range(0, 113):
        if item in zz:
            plt.text(-1.5, item - .25, element[item], fontsize = 7)

    #Set up grid structure for adding second plot.
    plt.grid(True)

    #Read second ts file, use variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile2)

    #Repeat plotting process from above, with data from the second ts file.
    ax2 = plt.subplot2grid((1, 2), (0, 1)) #positions the second plot on the grid
    for counter in np.arange(0, num_species):
        nn = aa[counter] - zz[counter]
        nn = int(round(nn))
        if colors:
            s = ax2.scatter(nn, zz[counter], s = 60, c = xmf[k, counter], cmap = 'YlGnBu_r', norm = matplotlib.colors.LogNorm(vmin=xmf[k].min(), vmax=xmf[k].max()), marker = 's')
        else:
            s = ax2.scatter(nn, zz[counter], s = 60, marker = 's')

    #Set space between the two plots.
    fig.subplots_adjust(wspace = 0.1)

    #Print thermodynamic data.
    words = ('Time = %s sec') % time[k]
    x_loc = 5
    y_loc = 37
    plt.text(x_loc, y_loc, words)
    words = ('Temperature = %s GK') % temperature[k]
    plt.text(x_loc, y_loc - 1.5, words)
    words = (r'$\rho =$ %s g/$\mathregular{cm^3}$') % density[k]
    plt.text(x_loc, y_loc - 3.75, words)

    #Configure and label axes, placing ticks at helium, oxygen, calcium, and nickel.
    plt.ylim(-1, 40)
    plt.xlim(-2, 40)

    special_elements = [2, 8, 20, 28]
    loc = plticker.FixedLocator(special_elements)
    plt.xticks(special_elements, special_elements)
    plt.yticks(special_elements, special_elements)

    #Add labels for each element.
    element = rtf.build_element_symbol()
    for item in range(0, 113):
        if item in zz:
            plt.text(-1.5, item - .25, element[item], fontsize = 7)

    plt.grid(True)

    #Show colorbar as legend.
    if colors:
        cb = plt.colorbar(s)
        cb.set_label('Log of Mass Fraction')

    plt.show()
