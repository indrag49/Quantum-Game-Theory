import numpy as np
import pandas as pd
import random
import itertools
import pylab


## Initialise the 1-Qubit system: {|0>, |1>}
Q0=np.array([1., 0.])
Q1=np.array([0., 1.])

## Initialise the 2-Qubit system: {|00>, |01>, |10>, |11>}
Q00=np.concatenate(np.outer(Q0, Q0))
Q01=np.concatenate(np.outer(Q0, Q1))
Q10=np.concatenate(np.outer(Q1, Q0))
Q11=np.concatenate(np.outer(Q1, Q1))

Dict={'0':Q0, '1':Q1, '00':Q00, '01':Q01, '10':Q10, '11':Q11}

X=[]
for n in range(3, 9):
    l=list(itertools.product([0, 1], repeat=n))
    X+=[l, ]
x=list(itertools.chain(*X))
L=[]
for i in x:
    L+=["".join([str(j) for j in i]), ]

for i in L:
    l=len(i)
    s='Q'+i
    s1=i[:-1]
    s2=i[-1]
    r=np.concatenate(np.outer(Dict[s1], Dict[s2]))
    globals()[s]=r
    Dict[s[1:]]=r

## Diagonal basis
Q_plus=np.array([1, 1])/np.sqrt(2)
Q_minus=np.array([1, -1])/np.sqrt(2)

## Circular basis
Q_clock=np.array([1, 1j])/np.sqrt(2)
Q_anticlock=np.array([1, -1j])/np.sqrt(2)

## Identity matrices
I2 = np.identity(2)
I4 = np.identity(2**2)
I8 = np.identity(2**3)
I16 = np.identity(2**4)
I32 = np.identity(2**5)

# Define Gate operations
def PauliX(n): return np.array(([0., 1.], [1., 0.])).dot(n)
def PauliY(n): return np.array(([0., -1.0j], [1.0j, 0.])).dot(n)
def PauliZ(n): return np.array(([1., 0.], [0., -1.])).dot(n)

def Dagger(n): return np.conj(n).T

def Hadamard(n): return np.array(([1., 1.], [1., -1.])).dot(n)/np.sqrt(2)
def Hadamard4(n):
        """ This represents a 4X4 Hadamard Matrix, n is a 2-Qubit system """
        h=Hadamard(I2)
        x=np.concatenate([h, h], axis=1)
        y=np.concatenate([h, -h], axis=1)
        h4=np.concatenate([x, y])
        return h4.dot(n)
def Hadamard8(n):
        """ This represents a 8X8 Hadamard Matrix, n is a 3-Qubit system """
        h4=Hadamard4(I4)
        x=np.concatenate([h4, h4], axis=1)
        y=np.concatenate([h4, -h4], axis=1)
        h8=np.concatenate([x, y])
        return h8.dot(n)
def Hadamard16(n):
        """ This represents a 16X16 Hadamard Matrix, n is a 4-Qubit system """
        h8=Hadamard8(I8)
        x=np.concatenate([h8, h8], axis=1)
        y=np.concatenate([h8, -h8], axis=1)
        h16=np.concatenate([x, y])
        return h16.dot(n)
def Hadamard32(n):
        """ This represents a 32X32 Hadamard Matrix, n is a 5-Qubit system """
        h16=Hadamard16(I16)
        x=np.concatenate([h16, h16], axis=1)
        y=np.concatenate([h16, -h16], axis=1)
        h32=np.concatenate([x, y])
        return h32.dot(n)

def CNOT(n):
        """CNOT gate on 2-Qubit system with control qubit = 0 and target qubit = 1"""
        x=np.copy(I4)
        t=np.copy(x[2,])
        x[2,]=x[3,]
        x[3,]=t
        return x.dot(n)

def CNOT2_10(n):
        """CNOT gate on 2-Qubit system with control qubit = 1 and target qubit = 0"""
        H=Hadamard(I2)
        x=CNOT(I4)
        y=np.kron(H, H)
        return y.dot(x).dot(y).dot(n)
        
def CNOT3_01(n):
        """CNOT gate on 3-Qubit system with control qubit = 0 and target qubit = 1"""
        return (np.kron(CNOT(I4), I2)).dot(n)

def CNOT3_02(n):
        """CNOT gate on 3-Qubit system with control qubit = 0 and target qubit = 2"""
        return np.array(([1., 0., 0., 0., 0., 0., 0., 0.], [0., 1., 0., 0., 0., 0., 0., 0.], [0., 0., 1., 0., 0., 0., 0., 0.], [0., 0., 0., 1., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 1., 0., 0.], [0., 0., 0., 0., 1., 0., 0., 0.], [0., 0., 0., 0., 0., 0., 0., 1.], [0., 0., 0., 0., 0., 0., 1., 0.])).dot(n)

def CNOT3_10(n):
        """CNOT gate on 3-Qubit system with control qubit = 1 and target qubit = 0"""
        return np.kron(CNOT2_10(I4), I2).dot(n)

def CNOT3_12(n):
        """CNOT gate on 3-Qubit system with control qubit = 1 and target qubit = 2"""
        return np.kron(I2, CNOT(I4)).dot(n)

def CNOT3_20(n):
        """CNOT gate on 3-Qubit system with control qubit = 2 and target qubit = 0"""
        return np.array(([1., 0., 0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 1., 0., 0.], [0., 0., 1., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0., 0., 1.], [0., 0., 0., 0., 1., 0., 0., 0.], [0., 1., 0., 0., 0., 0., 0., 0.], [0., 0., 0., 0., 0., 0., 1., 0.], [0., 0., 0., 1., 0., 0., 0., 0.])).dot(n)

def CNOT3_21(n):
        """CNOT gate on 3-Qubit system with control qubit = 2 and target qubit = 1"""
        return np.kron(I2, CNOT2_10(I4)).dot(n)

def CNOT4_12(n): return np.kron(CNOT3_12(I8), I2)
def CNOT5_12(n): return np.kron(CNOT4_12(I16), I2)


def cPauliY(n): return np.array(([1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 0., -1.0j], [0., 0., 1.0j, 0.])).dot(n)
def cPauliZ(n): return np.array(([1., 0., 0., 0.], [0., 1., 0., 0.], [0., 0., 1., 0.], [0., 0., 0., -1.])).dot(n)

def Rotate(n, t): return np.array(([np.cos(t), np.sin(t)], [-np.sin(t), np.cos(t)])).dot(n)

def Phase(n): return np.array(([1., 0.], [0., 1.0j])).dot(n)
def PhaseDagger(n): return np.array(([1., 0.], [0., -1.0j])).dot(n)
def T(n): return np.array(([1., 0.], [0., np.exp(1.0j*np.pi/4)])).dot(n)
def TDagger(n): return np.array(([1., 0.], [0., np.exp(-1.0j*np.pi/4)])).dot(n)

def R(n, Lambda): return np.array(([1., 0.], [0., exp(1j*Lambda)])).dot(n)

def SWAP(n):
        """n is a 4X4 matrix"""
        x=np.copy(I4)
        t=np.copy(x[1, ])
        x[1,]=x[2,]
        x[2,]=t
        return x.dot(n)
def Toffoli(n):
        """n must be a 8X8 matrix"""
        x=np.copy(I8)
        t=np.copy(x[6,])
        x[6,]=x[7,]
        x[7,]=t
        return x.dot(n)
def Fredkin(n):
        """n must be a 8X8 matrix"""
        x=np.copy(I8)
        t=np.copy(x[5,])
        x[5,]=x[6,]
        x[6,]=t
        return x.dot(n)
def Ising(n, phi):
        pi=np.pi
        e=np.exp
        return np.array(([1., 0., 0., e(1j*(phi-pi/2))], [0., 1., -1.0j, 0.], [0., -1.0j, 1., 0.], [e(1j*(-phi-pi/2)), 0., 0., 1.])).dot(n)

## Preparation of the Bell States: {|Beta_00>, |Beta_01>, |Beta_10>, |Beta_11>}
def bell(Qubit1, Qubit2):
        """Qubit1 and Qubit2 must be 0 or 1"""
        h=Hadamard(Q0) if Qubit1==0 else Hadamard(Q1)        
        x=np.concatenate(np.outer(h, Q1 if Qubit2==1 else Q0))
        return CNOT(x)
    
def Walsh(n): return Hadamard(n)

def Walsh4(n):
    h=Hadamard(I2)
    w=np.kron(h, h)
    return w.dot(n)

def Walsh8(n):
    h=Hadamard(I2)
    w=np.kron(np.kron(h, h), h)
    return w.dot(n)

def Walsh16(n):
    h=Hadamard(I2)
    w=np.kron(np.kron(np.kron(h, h), h), h)
    return w.dot(n)

def Walsh32(n):
    h=Hadamard(I2)
    w=np.kron(np.kron(np.kron(np.kron(h, h), h), h), h)
    return w.dot(n)

def QFT(y):
    """ For a given state |y> the quantum Fourier transform is calculated with QFT(y) """
    Y=bin(y)[2:]
    n=len(Y)
    return sum([np.exp(2*np.pi*1j*x*y/(2**n))*Dict['0'*(n-len(bin(x)[2:]))+bin(x)[2:]] for x in range(2**n)])/np.sqrt(2**n)

Beta_00=bell(0, 0)
Beta_01=bell(0, 1)
Beta_10=bell(1, 0)
Beta_11=bell(1, 1)

## An alternate way to prepare the bell states
def B(Qubit):
        H=Hadamard(I2)
        b=CNOT(np.kron(H, I2))
        return b.dot(Qubit)
b00=B(Q00)
b01=B(Q01)
b10=B(Q10)
b11=B(Q11)

def measure(n):
    l=len(n)
    n=pd.DataFrame(n)
    values=[]
    for i in range(l): values+=[abs(n.iloc[i][0])**2, ]
    p=pd.DataFrame(values).T
    if l==2: p.columns=["|0>", "|1>"]
    elif l==2**2: p.columns=["|00>", "|01>", "|10>", "|11>"]
    elif l==2**3: p.columns=["|000>","|001>","|010>","|011>","|100>","|101>","|110>","|111>"]
    elif l==2**4: p.columns=["|0000>","|0001>","|0010>","|0011>","|0100>","|0101>","|0110>","|0111>","|1000>","|1001>","|1010>","|1011>","|1100>","|1101>","|1110>","|1111>"]
    elif l==2**5: p.columns=["|00000>","|00001>","|00010>","|00011>","|00100>","|00101>","|00110>","|00111>","|01000>","|01001>","|01010>","|01011>","|01100>","|01101>","|01110>","|01111>","|10000>","|10001>","|10010>","|10011>","|10100>","|10101>","|10110>","|10111>","|11000>","|11001>","|11010>","|11011>","|11100>","|11101>","|11110>","|11111>"]
    elif l==2**6: p.columns=['|000000>', '|000001>', '|000010>', '|000011>', '|000100>', '|000101>', '|000110>', '|000111>', '|001000>', '|001001>', '|001010>', '|001011>', '|001100>', '|001101>', '|001110>', '|001111>', '|010000>', '|010001>', '|010010>', '|010011>', '|010100>', '|010101>', '|010110>', '|010111>', '|011000>', '|011001>', '|011010>', '|011011>', '|011100>', '|011101>', '|011110>', '|011111>', '|100000>', '|100001>', '|100010>', '|100011>', '|100100>', '|100101>', '|100110>', '|100111>', '|101000>', '|101001>', '|101010>', '|101011>', '|101100>', '|101101>', '|101110>', '|101111>', '|110000>', '|110001>', '|110010>', '|110011>', '|110100>', '|110101>', '|110110>', '|110111>', '|111000>', '|111001>', '|111010>', '|111011>', '|111100>', '|111101>', '|111110>', '|111111>']
    elif l==2**7: p.columns=['|0000000>', '|0000001>', '|0000010>', '|0000011>', '|0000100>', '|0000101>', '|0000110>', '|0000111>', '|0001000>', '|0001001>', '|0001010>', '|0001011>', '|0001100>', '|0001101>', '|0001110>', '|0001111>', '|0010000>', '|0010001>', '|0010010>', '|0010011>', '|0010100>', '|0010101>', '|0010110>', '|0010111>', '|0011000>', '|0011001>', '|0011010>', '|0011011>', '|0011100>', '|0011101>', '|0011110>', '|0011111>', '|0100000>', '|0100001>', '|0100010>', '|0100011>', '|0100100>', '|0100101>', '|0100110>', '|0100111>', '|0101000>', '|0101001>', '|0101010>', '|0101011>', '|0101100>', '|0101101>', '|0101110>', '|0101111>', '|0110000>', '|0110001>', '|0110010>', '|0110011>', '|0110100>', '|0110101>', '|0110110>', '|0110111>', '|0111000>', '|0111001>', '|0111010>', '|0111011>', '|0111100>', '|0111101>', '|0111110>', '|0111111>', '|1000000>', '|1000001>', '|1000010>', '|1000011>', '|1000100>', '|1000101>', '|1000110>', '|1000111>', '|1001000>', '|1001001>', '|1001010>', '|1001011>', '|1001100>', '|1001101>', '|1001110>', '|1001111>', '|1010000>', '|1010001>', '|1010010>', '|1010011>', '|1010100>', '|1010101>', '|1010110>', '|1010111>', '|1011000>', '|1011001>', '|1011010>', '|1011011>', '|1011100>', '|1011101>', '|1011110>', '|1011111>', '|1100000>', '|1100001>', '|1100010>', '|1100011>', '|1100100>', '|1100101>', '|1100110>', '|1100111>', '|1101000>', '|1101001>', '|1101010>', '|1101011>', '|1101100>', '|1101101>', '|1101110>', '|1101111>', '|1110000>', '|1110001>', '|1110010>', '|1110011>', '|1110100>', '|1110101>', '|1110110>', '|1110111>', '|1111000>', '|1111001>', '|1111010>', '|1111011>', '|1111100>', '|1111101>', '|1111110>', '|1111111>']
    else: p.columns=['|00000000>', '|00000001>', '|00000010>', '|00000011>', '|00000100>', '|00000101>', '|00000110>', '|00000111>', '|00001000>', '|00001001>', '|00001010>', '|00001011>', '|00001100>', '|00001101>', '|00001110>', '|00001111>', '|00010000>', '|00010001>', '|00010010>', '|00010011>', '|00010100>', '|00010101>', '|00010110>', '|00010111>', '|00011000>', '|00011001>', '|00011010>', '|00011011>', '|00011100>', '|00011101>', '|00011110>', '|00011111>', '|00100000>', '|00100001>', '|00100010>', '|00100011>', '|00100100>', '|00100101>', '|00100110>', '|00100111>', '|00101000>', '|00101001>', '|00101010>', '|00101011>', '|00101100>', '|00101101>', '|00101110>', '|00101111>', '|00110000>', '|00110001>', '|00110010>', '|00110011>', '|00110100>', '|00110101>', '|00110110>', '|00110111>', '|00111000>', '|00111001>', '|00111010>', '|00111011>', '|00111100>', '|00111101>', '|00111110>', '|00111111>', '|01000000>', '|01000001>', '|01000010>', '|01000011>', '|01000100>', '|01000101>', '|01000110>', '|01000111>', '|01001000>', '|01001001>', '|01001010>', '|01001011>', '|01001100>', '|01001101>', '|01001110>', '|01001111>', '|01010000>', '|01010001>', '|01010010>', '|01010011>', '|01010100>', '|01010101>', '|01010110>', '|01010111>', '|01011000>', '|01011001>', '|01011010>', '|01011011>', '|01011100>', '|01011101>', '|01011110>', '|01011111>', '|01100000>', '|01100001>', '|01100010>', '|01100011>', '|01100100>', '|01100101>', '|01100110>', '|01100111>', '|01101000>', '|01101001>', '|01101010>', '|01101011>', '|01101100>', '|01101101>', '|01101110>', '|01101111>', '|01110000>', '|01110001>', '|01110010>', '|01110011>', '|01110100>', '|01110101>', '|01110110>', '|01110111>', '|01111000>', '|01111001>', '|01111010>', '|01111011>', '|01111100>', '|01111101>', '|01111110>', '|01111111>', '|10000000>', '|10000001>', '|10000010>', '|10000011>', '|10000100>', '|10000101>', '|10000110>', '|10000111>', '|10001000>', '|10001001>', '|10001010>', '|10001011>', '|10001100>', '|10001101>', '|10001110>', '|10001111>', '|10010000>', '|10010001>', '|10010010>', '|10010011>', '|10010100>', '|10010101>', '|10010110>', '|10010111>', '|10011000>', '|10011001>', '|10011010>', '|10011011>', '|10011100>', '|10011101>', '|10011110>', '|10011111>', '|10100000>', '|10100001>', '|10100010>', '|10100011>', '|10100100>', '|10100101>', '|10100110>', '|10100111>', '|10101000>', '|10101001>', '|10101010>', '|10101011>', '|10101100>', '|10101101>', '|10101110>', '|10101111>', '|10110000>', '|10110001>', '|10110010>', '|10110011>', '|10110100>', '|10110101>', '|10110110>', '|10110111>', '|10111000>', '|10111001>', '|10111010>', '|10111011>', '|10111100>', '|10111101>', '|10111110>', '|10111111>', '|11000000>', '|11000001>', '|11000010>', '|11000011>', '|11000100>', '|11000101>', '|11000110>', '|11000111>', '|11001000>', '|11001001>', '|11001010>', '|11001011>', '|11001100>', '|11001101>', '|11001110>', '|11001111>', '|11010000>', '|11010001>', '|11010010>', '|11010011>', '|11010100>', '|11010101>', '|11010110>', '|11010111>', '|11011000>', '|11011001>', '|11011010>', '|11011011>', '|11011100>', '|11011101>', '|11011110>', '|11011111>', '|11100000>', '|11100001>', '|11100010>', '|11100011>', '|11100100>', '|11100101>', '|11100110>', '|11100111>', '|11101000>', '|11101001>', '|11101010>', '|11101011>', '|11101100>', '|11101101>', '|11101110>', '|11101111>', '|11110000>', '|11110001>', '|11110010>', '|11110011>', '|11110100>', '|11110101>', '|11110110>', '|11110111>', '|11111000>', '|11111001>', '|11111010>', '|11111011>', '|11111100>', '|11111101>', '|11111110>', '|11111111>']
    return p

def plot_measure(n):
        n.plot.bar()
        pylab.xlabel("Qubits")
        pylab.ylabel("Probabilitis")
        pylab.title("Probability Distribution")
        pylab.show()
