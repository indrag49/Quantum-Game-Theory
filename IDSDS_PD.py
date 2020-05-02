from QSimulator import *

alpha, beta, gamma, delta=3, 1, 0, 5
Alice=np.array(([alpha, gamma], [delta, beta]))
Bob=np.array(([alpha, delta], [gamma, beta]))

def IDSDS(P1, P2):
    check1, check2=0, 0
    while True:
        for i in range(len(P1)):
            for j in range(len(P1)):
                if list(P1[i, :]<P1[j, :]).count(True)==len(P1[0]):
                    P1=np.delete(P1, (i), axis=0)
                    P2=np.delete(P2, (i), axis=0)
                    check1=1
                    break
            if check1==1: break

        if check1==0:
            for i in range(len(P2[0])):
                for j in range(len(P2[0])):
                    if list(P2[:, i]<P2[:, j]).count(True)==len(P2):
                        P2=np.delete(P2, (i), axis=1)
                        P1=np.delete(P1, (i), axis=1)
                        check2=1
                        break
                if check2==1:break
        print(P1)
        print(P2)
        if check1==0 and check2==0: break
        check1, check2=0, 0
    return (P1, P2)

IDSDS(Alice, Bob)
