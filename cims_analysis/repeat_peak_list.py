start_time = time.time()
directory = 
f = 

i = 0
while i < 3:
    
    #initialise peak list object
    pl = Peaklist("I", ["CH2","O"]) 

    # Load data into dataframe
    if i > 0:
        fname = f+"_kmd%i.txt" % i
    else:
        fname = f+".txt"

    print "loading %s" % fname

    pl.Load_peak_list(directory, fname, sep="\t") 

    # Calculates mass defect and kendrick mass defects for 
    # passed kendrick bases. Matches kendrick mass defects 
    # and gets suggested formulas and error on assignment.
    pl.Run() 

    i += 1

    # Generate new peaklist based on suggestions
    pl.New_peaklist().to_csv(
        directory + f + "_kmd%i.txt" % i, 
        sep="\t", index=False)


print "--- Time taken: %s s ---" % (round(time.time() - start_time, 3))
