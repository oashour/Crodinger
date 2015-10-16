# Import modules
import os
import time
import subprocess
import math

def main():
    
    os.chdir("..")
    dt = "0.0001"
    devnull = open('/dev/null', 'w')

    w = [0.01]

    q = []
    i = 0
    k = 0
    while k < 0.45:
        k = 0.05+i*0.0125;
        q.append(k);
        i = i+1;
    q.append(0.46)
    # Run GNU Make
    print("Running GNU Make\n")
    subprocess.call("make")

    # Allocate dt for multiple runs
    #q = [0.4296875, 0.468, 0.486, 0.4895]
    #q = []
    tm = [[0 for x in range(len(q))] for x in range(len(w))] 

    for j in range(len(w)):
        for i in range(len(q)):
            lamb = 2*math.sqrt(2*q[i])*math.sqrt(1-2*q[i])
            omega = 2*math.sqrt(1-2*q[i])
            t0 = (-math.log(w[j]/omega) + math.log(lamb))/lamb
            tm[j][i] = int(t0)+2

    for j in range(len(w)):
        print(tm[j])

    # Check if results file exists
#    fileName = "output/results_q.txt"
#    if os.path.exists(fileName):
#        os.remove(fileName)
    
    #Run code multiple times
    for j in range(len(w)):
        for i in range(len(q)):
            print("Running code with w={}, q = {}".format(w[j], q[i]))
            string = "echo \"{} {} {} {} {} {}\"".format(dt, tm[j][i], 
                                           '256', w[j], q[i], '2', 'X')

            p1 = subprocess.Popen(string, stdout=subprocess.PIPE, shell=True)
            p2 = subprocess.Popen("./nlse", stdin=p1.stdout, stdout=subprocess.PIPE)
            p1.stdout.close()  # Allow p1 to receive a SIGPIPE if p2 exits.
            output = p2.communicate()[0]
 
            print("Running GNU Octave")
            subprocess.call(["octave", "main.m"], stdout=devnull)
            print("")

    print("Done looping!\n")

if __name__ == "__main__":
    main()

