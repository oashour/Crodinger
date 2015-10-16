# Import modules
import os
import time
import subprocess

devnull = open('/dev/null', 'w')

def run(order, type1, A):
    
    # Allocate dt for multiple runs
    start = 0.0001;
    num = 100
    dt = []
    for i in range(num):
        dt.append(start*(i+1))
    
    dtS = []
    for i in range(num):
        dtS.append(format(dt[i], '.4f'))

    # Check if results file exists
    fileName = "results_{}{}.txt".format(order, type1)
    if os.path.exists(fileName):
        os.remove(fileName)
    
    # Run code multiple times
    for i in dtS:
        print("Running code with dt={}, order={}, type ={}".format(i, order, type1))
        subprocess.call(["./nlse", i, str(order), str(type1), str(A)], stdout=devnull)
        print("Running GNU Octave")
        subprocess.call(["octave", "main.m"], stdout=devnull)
        print("")

    print("Done looping!\n")


def main():

    # Run GNU Make
    print("Running GNU Make\n")
    subprocess.call("make", stdout=devnull)

    A=1e-1
    print("Second order.")
    run(2, 'X', A)
    print("Fourth order.")
    run(4, 'S', A)

    # # Higher orders
    # orders = [4, 6, 8]
    # types = ['S', 'M'] # Only M and S. Remove the one you do not need
# 
    # for type1 in types:
        # for order in orders: 
            # print("Running order={}, type={}\n\n".format(order, type1))
            # run(order, type1, A)
    # print("Done!")

if __name__ == "__main__":
    main()

