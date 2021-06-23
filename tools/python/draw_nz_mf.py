def draw_nz_mf(datafile, k, num_species, colors = False):

    '''
        Colorbar is currently not working.
        Inputs: datafile = must be a ts file, formatted as 'datafile'
                k = desired timestep for thermodynamic and mass fraction data
                num_species = number of plotted species desired
                colors = turns colorbar on or off, defaults to False
        Outputs: scatterplot of number of neutrons vs number of protons, color coded for mass fraction
                 thermodynamic data at specific time
    '''
    import read_ts_file as rtf #module to read ts file
    import numpy as np #module to provide calculation support
    import matplotlib.pyplot as plt #plotting module
    import matplotlib.cm as cm #colormap module
    import matplotlib.ticker as plticker #module to make and assign tick marks on axes
    import matplotlib.colors #colors module

    #Create plot space.
    fig = plt.figure(1)

    #Read ts file and use variables.
    zz, aa, xmf, time, temperature, density, timestep, edot, flx_end, flx = rtf.read_ts_file(datafile) #Each of these variables is returned by the function read_ts_file in the module rtf, so you just assign them to the same names here

    #Plot proton number vs neutron number and color code for mass fraction (log).
    for counter in np.arange(0, num_species): #generates a counter that goes from 0 to the number of species
        nn = aa[counter] - zz[counter] #Determine number of neutrons for each species. aa is atomic weight, and zz is proton number.
        nn = int(round(nn)) #Round the number of neutrons, and cast it to an integer so that it can be used for plotting
        #If colors has been set as true, plot a scatterplot of nn vs zz.
        #the first two terms are the x and y of your plot
        #nn = neutron number
        #zz[counter] = proton number, counted through for every species
        #s = a size determiner, not all that important
        #marker = lets you choose the point shape, 's' means square. It should default to round.
        #c = sets color. Here, I have it do a gradient based on mass fraction.
        #cmap = lets you choose colormap (there's a list of them on the matplotlib website)
        #norm = normalizes the data to fit between 0 and 1, here I use a log scale with largest and smallest mass fractions at timestep as parameters
        if colors == True:
            s = plt.scatter(nn, zz[counter], s = 50, c = xmf[k, counter], cmap = 'YlOrBr', norm = matplotlib.colors.LogNorm(vmin=xmf[k].min(), vmax=xmf[k].max()), marker = 's')
        else:
            s = plt.scatter(nn, zz[counter], s = 50, marker = 's')

    if colors == True:
        cb = plt.colorbar(s) #Plot colorbar as legend.
        cb.set_label('Log of Mass Fraction')

    #Print thermodynamic data.
    words = ('Time = %s sec') % time[k] #Generates the time at the specified timestep
    x_loc = 5 #Arbitrary x and y coordinates
    y_loc = 37
    plt.text(x_loc, y_loc, words) #Writes the time at specified x and y.
    words = ('Temperature = %s GK') % temperature[k] #Generates the temperature at the specified timestep
    plt.text(x_loc, y_loc - 1.5, words) #Writes temperature just below time.
    words = (r'$\rho =$ %s g/$\mathregular{cm^3}$') % density[k] #Generates the density at the specified timestep
    plt.text(x_loc, y_loc - 3.75, words) #Writes density just before temperature
    
    #Format x and y axes.
    plt.ylim(-1, 40) #Y axis will range from -1 to 40, and so will the x axis. They're labeled accordingly.
    plt.ylabel('Z (Proton Number)')
    plt.xlim(-2, 40)
    plt.xlabel('N (Neutron Number)')

    #Set axis ticks.
    special_elements = [2, 8, 20, 28] #List of values where I want ticks
    loc = plticker.FixedLocator(special_elements) #Says that I want fixed, specified locations for the ticks
    plt.xticks(special_elements, special_elements) #Labels axes with the numerical values of the ticks at those positions (prints each item in the list at the value of that item)
    plt.yticks(special_elements, special_elements)

    #Add labels for each element.
    element = rtf.build_element_symbol() #build_element_symbol returns a list of each element name, which you're using here
    for item in range(0, 118): #counts through all elements, and if their proton numbers are in the list of proton numbers plotted, their name appears along the y axis.
        if item in zz:
            plt.text(-1.5, item - .25, element[item], fontsize = 7)

    #Add grid lines,
    plt.grid(True)

    #Add title, just the name of the datafile input.
    plt.title(datafile)

    #Show plot, which makes the plot appear.
    plt.show()
