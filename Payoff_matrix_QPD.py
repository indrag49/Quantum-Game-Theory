from QPD import *

moves={1:I2, 2:PauliX(I2), 3:Hadamard(I2), 4:PauliZ(I2)}

w, x, y, z=3, 1, 0, 5
Alice=np.zeros([4, 4])
Bob=np.zeros([4, 4])
for i in range(1, 5):
    for j in range(1, 5):
        X=qpd(moves[i], moves[j], w, x, y, z)
        Alice[i-1, j-1], Bob[i-1, j-1]=X[0], X[1]

print(Alice, Bob)
