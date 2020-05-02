from QSimulator import *

def penny_flip4():

    #Alice prepares the initial state
    u=Q0
    d=Q1
    A=d
    state=A

    #Bob plays his move: Hadamard operator
    B=Hadamard(state)
    state=B

    #Now again Alice plays
    r=random.uniform(0, 1)
    A=PauliX(state) if r<0.5 else I2.dot(state)
    state=A

    #Last move by Bob: again Hadamard operator
    B=Hadamard(state)
    state=B
    plot_measure(measure(state))

penny_flip4()
