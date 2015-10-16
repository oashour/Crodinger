# Import modules
import os
import time
import subprocess


def main():

    devnull = open('/dev/null', 'w')

    # Run GNU Make
    print("Running GNU Make\n")
    subprocess.call("make", stdout=devnull)

    # Allocate dt for multiple runs
    w = [2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 9.0]
    #for i in range(num):
    #    w.append(start*(i+1))
    
    wS = []
    for i in range(len(w)):
        wS.append(format(w[i], '.4f'))

    # Check if results file exists
    fileName = "results_O2X.txt"
    if os.path.exists(fileName):
        os.remove(fileName)
    
    # Run code multiple times
    for i in wS:
        print("Running code with w={}".format(i))
        subprocess.call(["./nlse", '0.0001', '2', 'X', i], stdout=devnull)
        print("Running GNU Octave")
        subprocess.call(["octave", "main.m"], stdout=devnull)
        print("")

    print("Done looping!\n")

if __name__ == "__main__":
    main()

