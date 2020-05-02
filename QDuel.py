from QSimulator import *
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d
import math
from tkinter import *


pi=np.pi
e=np.exp
cos=np.cos
acos=np.arccos
sin=np.sin
sqrt=np.sqrt

def entryfield(parent,label,labelwidth=12,**packopts):
    f = Frame(parent)
    f.pack(**packopts)
    l = Label(f,text=label,width=labelwidth)
    l.pack(side=LEFT,anchor=W)
    value = DoubleVar(f)
    e = Entry(f,textvariable=value)
    e.pack(side=RIGHT,fill=X,expand=True)
    return lambda: value.get()

def entryfieldInt(parent,label,labelwidth=12,**packopts):
    f = Frame(parent)
    f.pack(**packopts)
    l = Label(f,text=label,width=labelwidth)
    l.pack(side=LEFT,anchor=W)
    value = IntVar(f)
    e = Entry(f,textvariable=value)
    e.pack(side=RIGHT,fill=X,expand=True)
    return lambda: value.get()

def gqduels(Psi, n, a, b, alpha1, alpha2, beta1, beta2):

    Ab=(e(-1j*alpha1)*cos(acos(sqrt(a)))*Q11[:, np.newaxis]+1j*e(1j*beta1)*sin(acos(sqrt(a)))*Q10[:, np.newaxis])*Q11+(e(1j*alpha1)*cos(acos(sqrt(a)))*Q10[:, np.newaxis]+1j*e(-1j*beta1)*sin(acos(sqrt(a)))*Q11[:, np.newaxis])*Q10+Q00[:, np.newaxis]*Q00+Q01[:, np.newaxis]*Q01
    Ba=(e(-1j*alpha2)*cos(acos(sqrt(b)))*Q11[:, np.newaxis]+1j*e(1j*beta2)*sin(acos(sqrt(b)))*Q01[:, np.newaxis])*Q11+(e(1j*alpha2)*cos(acos(sqrt(b)))*Q01[:, np.newaxis]+1j*e(-1j*beta2)*sin(acos(sqrt(b)))*Q11[:, np.newaxis])*Q01+Q00[:, np.newaxis]*Q00+Q10[:, np.newaxis]*Q10
    A=np.linalg.matrix_power((Ba.dot(Ab)), n)
        
    Psif=A.dot(Psi)
        
    #return (1+abs(np.conj(Q10).dot(Psif))**2-abs(np.conj(Q01).dot(Psif))**2)/2.
    return abs(np.conj(Q10).dot(Psif))**2 + 0.5*abs(np.conj(Q11).dot(Psif))**2

def qduels2(Psi, a, b, alpha1, alpha2, beta1, beta2):
    #The improvement plots do not depend on the beta values
    Ab=(e(-1j*alpha1)*cos(acos(sqrt(a)))*Q11[:, np.newaxis]+1j*e(1j*beta1)*sin(acos(sqrt(a)))*Q10[:, np.newaxis])*Q11+(e(1j*alpha1)*cos(acos(sqrt(a)))*Q10[:, np.newaxis]+1j*e(-1j*beta1)*sin(acos(sqrt(a)))*Q11[:, np.newaxis])*Q10+Q00[:, np.newaxis]*Q00+Q01[:, np.newaxis]*Q01
    Ba=(e(-1j*alpha2)*cos(acos(sqrt(b)))*Q11[:, np.newaxis]+1j*e(1j*beta2)*sin(acos(sqrt(b)))*Q01[:, np.newaxis])*Q11+(e(1j*alpha2)*cos(acos(sqrt(b)))*Q01[:, np.newaxis]+1j*e(-1j*beta2)*sin(acos(sqrt(b)))*Q11[:, np.newaxis])*Q01+Q00[:, np.newaxis]*Q00+Q10[:, np.newaxis]*Q10
    
    C=Ba.dot(Ba.dot(Ab))
    Psi2=C.dot(Psi)
    
    return abs(np.conj(Q10).dot(Psi2))**2 + 0.5*abs(np.conj(Q11).dot(Psi2))**2

def qduels3(Psi, a, b, alpha1, alpha2, beta1, beta2):
    #The improvement plots do not depend on the beta values
    Ab=(e(-1j*alpha1)*cos(acos(sqrt(a)))*Q11[:, np.newaxis]+1j*e(1j*beta1)*sin(acos(sqrt(a)))*Q10[:, np.newaxis])*Q11+(e(1j*alpha1)*cos(acos(sqrt(a)))*Q10[:, np.newaxis]+1j*e(-1j*beta1)*sin(acos(sqrt(a)))*Q11[:, np.newaxis])*Q10+Q00[:, np.newaxis]*Q00+Q01[:, np.newaxis]*Q01
    Ba=(e(-1j*alpha2)*cos(acos(sqrt(b)))*Q11[:, np.newaxis]+1j*e(1j*beta2)*sin(acos(sqrt(b)))*Q01[:, np.newaxis])*Q11+(e(1j*alpha2)*cos(acos(sqrt(b)))*Q01[:, np.newaxis]+1j*e(-1j*beta2)*sin(acos(sqrt(b)))*Q11[:, np.newaxis])*Q01+Q00[:, np.newaxis]*Q00+Q10[:, np.newaxis]*Q10
    
    C=Ab.dot(Ba.dot(Ab))
    Psi2=C.dot(Psi)
    
    return abs(np.conj(Q10).dot(Psi2))**2 + 0.5*abs(np.conj(Q11).dot(Psi2))**2

def twoPerson():
    Root1=Tk()

    def initialState():
        Root2=Tk()

        def q01():
            Root3=Tk()
            Psi=Q01
            
            def plot1():
                Root4=Tk()

                def simulate(n, a, b, beta1, beta2):

                    x=np.linspace(-pi, pi, 20)
                    y=np.linspace(-pi, pi, 20)

                    X, Y=np.meshgrid(x, y)

                    Z1=[]
                    Z2=[]
                    z1=[]
                    for i in x:
                        for j in y:
                            z1+=[gqduels(Psi, n, a, b, i, j, beta1, beta2), ]
                        Z1+=[z1, ]
                        Z2+=[1-np.array(z1), ]
                        z1=[]
                    Z1=np.array(Z1)
                    Z2=np.array(Z2)
                    X, Y=np.meshgrid(x, y)
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 2, 1, projection='3d')
                    ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Alice}}$")
                    ax = fig.add_subplot(1, 2, 2, projection='3d')
                    ax.plot_surface(X, Y, Z2, rstride=1, cmap='Reds', cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Bob}}$")
                    plt.show()

                def plot1_gui(parent):
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    n=entryfieldInt(parent, "n:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot1_gui(Root4)
                Root4.mainloop()

            def plot2():
                Root4=Tk()

                def simulate(n, a, b, alpha1, alpha2, beta1, beta2):

                    

                    X=np.arange(1, n+1)
                    Y1=[gqduels(Psi, i, a, b, alpha1, alpha2, beta1, beta2) for i in X]
                    Y2=1-np.array(Y1)
                    
                    pylab.subplot(1, 2, 1)
                    pylab.plot(X, Y1)
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Alice}}$")

                    pylab.subplot(1, 2, 2)
                    pylab.plot(X, Y2, 'r-')
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Bob}}$")
                    
                    pylab.show()

                def plot2_gui(parent):
                    n=entryfieldInt(parent, "n(no. of rounds):", side=TOP, fill=X)
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), alpha1(), alpha2(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)
                plot2_gui(Root4)
                Root4.mainloop()

            def plot3():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels2(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{Alice Diff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot3_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot3_gui(Root4)
                Root4.mainloop()

            def plot4():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels3(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{Bob Diff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot4_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot4_gui(Root4)
                Root4.mainloop()
                    

            def q01_gui(parent):
                B1=Button(parent, text='Plot I', bg="blue", command=plot1)
                B1.pack(side=TOP, fill=X)
                B2=Button(parent, text='Plot II', bg="red", command=plot2)
                B2.pack(side=TOP, fill=X)
                B3=Button(parent, text='Plot III', bg="green", command=plot3)
                B3.pack(side=TOP, fill=X)
                B4=Button(parent, text='Plot IV', bg="orange", command=plot4)
                B4.pack(side=TOP, fill=X)

            q01_gui(Root3)
            Root3.mainloop()
            
        def q10():
            Root3=Tk()
            Psi=Q10
            
            def plot1():
                Root4=Tk()

                def simulate(n, a, b, beta1, beta2):

                    x=np.linspace(-pi, pi, 20)
                    y=np.linspace(-pi, pi, 20)

                    X, Y=np.meshgrid(x, y)

                    Z1=[]
                    Z2=[]
                    z1=[]
                    for i in x:
                        for j in y:
                            z1+=[gqduels(Psi, n, a, b, i, j, beta1, beta2), ]
                        Z1+=[z1, ]
                        Z2+=[1-np.array(z1), ]
                        z1=[]
                    Z1=np.array(Z1)
                    Z2=np.array(Z2)
                    X, Y=np.meshgrid(x, y)
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 2, 1, projection='3d')
                    ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Alice}}$")
                    ax = fig.add_subplot(1, 2, 2, projection='3d')
                    ax.plot_surface(X, Y, Z2, rstride=1, cmap='Reds', cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Bob}}$")
                    plt.show()

                def plot1_gui(parent):
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    n=entryfieldInt(parent, "n:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot1_gui(Root4)
                Root4.mainloop()

            def plot2():
                Root4=Tk()

                def simulate(n, a, b, alpha1, alpha2, beta1, beta2):

                    

                    X=np.arange(1, n+1)
                    Y1=[gqduels(Psi, i, a, b, alpha1, alpha2, beta1, beta2) for i in X]
                    Y2=1-np.array(Y1)
                    
                    pylab.subplot(1,2, 1)
                    pylab.plot(X, Y1)
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Alice}}$")

                    pylab.subplot(1, 2, 2)
                    pylab.plot(X, Y2, 'r-')
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Bob}}$")
                    
                    pylab.show()

                def plot2_gui(parent):
                    n=entryfieldInt(parent, "n(no. of rounds):", side=TOP, fill=X)
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), alpha1(), alpha2(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot2_gui(Root4)
                Root4.mainloop()

            def plot3():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels2(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{Alice Diff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot3_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot3_gui(Root4)
                Root4.mainloop()

            def plot4():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels3(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{BobDiff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot4_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot4_gui(Root4)
                Root4.mainloop()
                    

            def q10_gui(parent):
                B1=Button(parent, text='Plot I', bg="blue", command=plot1)
                B1.pack(side=TOP, fill=X)
                B2=Button(parent, text='Plot II', bg="red", command=plot2)
                B2.pack(side=TOP, fill=X)
                B3=Button(parent, text='Plot III', bg="green", command=plot3)
                B3.pack(side=TOP, fill=X)
                B4=Button(parent, text='Plot IV', bg="orange", command=plot4)
                B4.pack(side=TOP, fill=X)

            q10_gui(Root3)
            Root3.mainloop()
        
        def q11():
            Root3=Tk()
            Psi=Q11
            
            def plot1():
                Root4=Tk()

                def simulate(n, a, b, beta1, beta2):

                    x=np.linspace(-pi, pi, 20)
                    y=np.linspace(-pi, pi, 20)

                    X, Y=np.meshgrid(x, y)

                    Z1=[]
                    Z2=[]
                    z1=[]
                    for i in x:
                        for j in y:
                            z1+=[gqduels(Psi, n, a, b, i, j, beta1, beta2), ]
                        Z1+=[z1, ]
                        Z2+=[1-np.array(z1), ]
                        z1=[]
                    Z1=np.array(Z1)
                    Z2=np.array(Z2)
                    X, Y=np.meshgrid(x, y)
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 2, 1, projection='3d')
                    ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Alice}}$")
                    ax = fig.add_subplot(1, 2, 2, projection='3d')
                    ax.plot_surface(X, Y, Z2, rstride=1, cmap='Reds', cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Bob}}$")
                    plt.show()

                def plot1_gui(parent):
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    n=entryfieldInt(parent, "n:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot1_gui(Root4)
                Root4.mainloop()

            def plot2():
                Root4=Tk()

                def simulate(n, a, b, alpha1, alpha2, beta1, beta2):

                    

                    X=np.arange(1, n+1)
                    Y1=[gqduels(Psi, i, a, b, alpha1, alpha2, beta1, beta2) for i in X]
                    Y2=1-np.array(Y1)
                    
                    pylab.subplot(1, 2, 1)
                    pylab.plot(X, Y1)
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Alice}}$")

                    pylab.subplot(1, 2, 2)
                    pylab.plot(X, Y2, 'r-')
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Bob}}$")
                    
                    pylab.show()

                def plot2_gui(parent):
                    n=entryfieldInt(parent, "n(no. of rounds):", side=TOP, fill=X)
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), alpha1(), alpha2(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot2_gui(Root4)
                Root4.mainloop()

            def plot3():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels2(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{AliceDiff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot3_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot3_gui(Root4)
                Root4.mainloop()

            def plot4():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels3(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{BobDiff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot4_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot4_gui(Root4)
                Root4.mainloop()
                    

            def q11_gui(parent):
                B1=Button(parent, text='Plot I', bg="blue", command=plot1)
                B1.pack(side=TOP, fill=X)
                B2=Button(parent, text='Plot II', bg="red", command=plot2)
                B2.pack(side=TOP, fill=X)
                B3=Button(parent, text='Plot III', bg="green", command=plot3)
                B3.pack(side=TOP, fill=X)
                B4=Button(parent, text='Plot IV', bg="orange", command=plot4)
                B4.pack(side=TOP, fill=X)

            q11_gui(Root3)
            Root3.mainloop()

        def q_superpose1():
            Root3=Tk()
            qs=(Q0 + Q1)/np.sqrt(2)
            Psi=np.concatenate(np.outer(Q1, qs))
            
            def plot1():
                Root4=Tk()

                def simulate(n, a, b, beta1, beta2):

                    x=np.linspace(-pi, pi, 20)
                    y=np.linspace(-pi, pi, 20)

                    X, Y=np.meshgrid(x, y)

                    Z1=[]
                    Z2=[]
                    z1=[]
                    for i in x:
                        for j in y:
                            z1+=[gqduels(Psi, n, a, b, i, j, beta1, beta2), ]
                        Z1+=[z1, ]
                        Z2+=[1-np.array(z1), ]
                        z1=[]
                    Z1=np.array(Z1)
                    Z2=np.array(Z2)
                    X, Y=np.meshgrid(x, y)
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 2, 1, projection='3d')
                    ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Alice}}$")
                    ax = fig.add_subplot(1, 2, 2, projection='3d')
                    ax.plot_surface(X, Y, Z2, rstride=1, cmap='Reds', cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Bob}}$")
                    plt.show()

                def plot1_gui(parent):
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    n=entryfieldInt(parent, "n:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot1_gui(Root4)
                Root4.mainloop()

            def plot2():
                Root4=Tk()

                def simulate(n, a, b, alpha1, alpha2, beta1, beta2):

                    

                    X=np.arange(1, n+1)
                    Y1=[gqduels(Psi, i, a, b, alpha1, alpha2, beta1, beta2) for i in X]
                    Y2=1-np.array(Y1)
                    
                    pylab.subplot(1, 2, 1)
                    pylab.plot(X, Y1)
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Alice}}$")

                    pylab.subplot(1, 2, 2)
                    pylab.plot(X, Y2, 'r-')
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Bob}}$")
                    
                    pylab.show()

                def plot2_gui(parent):
                    n=entryfieldInt(parent, "n(no. of rounds):", side=TOP, fill=X)
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), alpha1(), alpha2(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot2_gui(Root4)
                Root4.mainloop()

            def plot3():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels2(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{AliceDiff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot3_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot3_gui(Root4)
                Root4.mainloop()

            def plot4():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels3(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{BobDiff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot4_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot4_gui(Root4)
                Root4.mainloop()
                    

            def q_superpose1_gui(parent):
                B1=Button(parent, text='Plot I', bg="blue", command=plot1)
                B1.pack(side=TOP, fill=X)
                B2=Button(parent, text='Plot II', bg="red", command=plot2)
                B2.pack(side=TOP, fill=X)
                B3=Button(parent, text='Plot III', bg="green", command=plot3)
                B3.pack(side=TOP, fill=X)
                B4=Button(parent, text='Plot IV', bg="orange", command=plot4)
                B4.pack(side=TOP, fill=X)

            q_superpose1_gui(Root3)
            Root3.mainloop()

        def q_superpose2():
            Root3=Tk()
            qs=(Q0 + Q1)/np.sqrt(2)
            Psi=np.concatenate(np.outer(qs, Q1))
            
            def plot1():
                Root4=Tk()

                def simulate(n, a, b, beta1, beta2):

                    x=np.linspace(-pi, pi, 20)
                    y=np.linspace(-pi, pi, 20)

                    X, Y=np.meshgrid(x, y)

                    Z1=[]
                    Z2=[]
                    z1=[]
                    for i in x:
                        for j in y:
                            z1+=[gqduels(Psi, n, a, b, i, j, beta1, beta2), ]
                        Z1+=[z1, ]
                        Z2+=[1-np.array(z1), ]
                        z1=[]
                    Z1=np.array(Z1)
                    Z2=np.array(Z2)
                    X, Y=np.meshgrid(x, y)
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 2, 1, projection='3d')
                    ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Alice}}$")
                    ax = fig.add_subplot(1, 2, 2, projection='3d')
                    ax.plot_surface(X, Y, Z2, rstride=1, cmap='Reds', cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Bob}}$")
                    plt.show()

                def plot1_gui(parent):
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    n=entryfieldInt(parent, "n:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot1_gui(Root4)
                Root4.mainloop()

            def plot2():
                Root4=Tk()

                def simulate(n, a, b, alpha1, alpha2, beta1, beta2):

                    

                    X=np.arange(1, n+1)
                    Y1=[gqduels(Psi, i, a, b, alpha1, alpha2, beta1, beta2) for i in X]
                    Y2=1-np.array(Y1)
                    
                    pylab.subplot(1, 2, 1)
                    pylab.plot(X, Y1)
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Alice}}$")

                    pylab.subplot(1, 2, 2)
                    pylab.plot(X, Y2, 'r-')
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Bob}}$")
                    
                    pylab.show()

                def plot2_gui(parent):
                    n=entryfieldInt(parent, "n(no. of rounds):", side=TOP, fill=X)
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), alpha1(), alpha2(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot2_gui(Root4)
                Root4.mainloop()

            def plot3():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels2(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{AliceDiff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot3_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot3_gui(Root4)
                Root4.mainloop()

            def plot4():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels3(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{BobDiff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot4_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot4_gui(Root4)
                Root4.mainloop()
                    

            def q_superpose2_gui(parent):
                B1=Button(parent, text='Plot I', bg="blue", command=plot1)
                B1.pack(side=TOP, fill=X)
                B2=Button(parent, text='Plot II', bg="red", command=plot2)
                B2.pack(side=TOP, fill=X)
                B3=Button(parent, text='Plot III', bg="green", command=plot3)
                B3.pack(side=TOP, fill=X)
                B4=Button(parent, text='Plot IV', bg="orange", command=plot4)
                B4.pack(side=TOP, fill=X)

            q_superpose2_gui(Root3)
            Root3.mainloop()

        def q_superpose3():
            Root3=Tk()
            qs=(Q0 + Q1)/np.sqrt(2)
            Psi=np.concatenate(np.outer(qs, qs))
            
            def plot1():
                Root4=Tk()

                def simulate(n, a, b, beta1, beta2):

                    x=np.linspace(-pi, pi, 20)
                    y=np.linspace(-pi, pi, 20)

                    X, Y=np.meshgrid(x, y)

                    Z1=[]
                    Z2=[]
                    z1=[]
                    for i in x:
                        for j in y:
                            z1+=[gqduels(Psi, n, a, b, i, j, beta1, beta2), ]
                        Z1+=[z1, ]
                        Z2+=[1-np.array(z1), ]
                        z1=[]
                    Z1=np.array(Z1)
                    Z2=np.array(Z2)
                    X, Y=np.meshgrid(x, y)
                    fig = plt.figure()
                    ax = fig.add_subplot(1, 2, 1, projection='3d')
                    ax.plot_surface(X, Y, Z1, rstride=1, cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Alice}}$")
                    ax = fig.add_subplot(1, 2, 2, projection='3d')
                    ax.plot_surface(X, Y, Z2, rstride=1, cmap='Reds', cstride=1, edgecolor='none')
                    ax.set_xlabel(r"$\alpha_1$")
                    ax.set_ylabel(r"$\alpha_2$")
                    ax.set_zlabel(r"$\overline{\pi_{Bob}}$")
                    plt.show()

                def plot1_gui(parent):
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    n=entryfieldInt(parent, "n:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot1_gui(Root4)
                Root4.mainloop()

            def plot2():
                Root4=Tk()

                def simulate(n, a, b, alpha1, alpha2, beta1, beta2):

                    

                    X=np.arange(1, n+1)
                    Y1=[gqduels(Psi, i, a, b, alpha1, alpha2, beta1, beta2) for i in X]
                    Y2=1-np.array(Y1)
                    
                    pylab.subplot(1, 2, 1)
                    pylab.plot(X, Y1)
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Alice}}$")

                    pylab.subplot(1, 2, 2)
                    pylab.plot(X, Y2, 'r-')
                    pylab.xlabel("rounds")
                    pylab.ylabel(r"$\overline{\pi_{Bob}}$")
                    
                    pylab.show()

                def plot2_gui(parent):
                    n=entryfieldInt(parent, "n(no. of rounds):", side=TOP, fill=X)
                    a=entryfield(parent, "a:", side=TOP, fill=X)
                    b=entryfield(parent, "b:", side=TOP, fill=X)
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    beta1=entryfield(parent, "beta1:", side=TOP, fill=X)
                    beta2=entryfield(parent, "beta2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="blue", command=lambda: simulate(n(), a(), b(), alpha1(), alpha2(), beta1(), beta2()))
                    B1.pack(side=TOP, fill=X)

                plot2_gui(Root4)
                Root4.mainloop()

            def plot3():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels2(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{AliceDiff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot3_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot3_gui(Root4)
                Root4.mainloop()

            def plot4():
                Root4=Tk()

                def simulate(alpha1, alpha2):
                    a=np.linspace(0, 1, 20)
                    b=np.linspace(0, 1, 20)


                    Z1=[]
                    z1=[]
                    z2=[]
                    z3=[]
                    for i in a:
                        for j in b:
                            z1+=[gqduels(Psi, 2, i, j, alpha1, alpha2, 0, 0), ]
                            z2+=[qduels3(Psi, i, j, alpha1, alpha2, 0, 0), ]
                            z3=np.array(z2)-np.array(z1)
                        Z1+=[z3, ]
                        z1=[]
                        z2=[]
                    Z1=np.array(Z1)
                    A, B=np.meshgrid(a, b)
                    fig = plt.figure()
                    ax = plt.axes(projection='3d')
                    ax.plot_surface(A, B, Z1, rstride=1, cstride=1,
                                    cmap='gnuplot', edgecolor='none')
                    ax.set_xlabel('a')
                    ax.set_ylabel('b')
                    ax.set_zlabel(r"$\overline{\pi_{BobDiff}}$")
                    plt.gca().invert_yaxis()
                    plt.show()

                def plot4_gui(parent):
                    alpha1=entryfield(parent, "alpha1:", side=TOP, fill=X)
                    alpha2=entryfield(parent, "alpha2:", side=TOP, fill=X)
                    B1=Button(parent, text='PLOT', bg="yellow", command=lambda: simulate(alpha1(), alpha2()))
                    B1.pack(side=TOP, fill=X)

                plot4_gui(Root4)
                Root4.mainloop()
                    

            def q_superpose3_gui(parent):
                B1=Button(parent, text='Plot I', bg="blue", command=plot1)
                B1.pack(side=TOP, fill=X)
                B2=Button(parent, text='Plot II', bg="red", command=plot2)
                B2.pack(side=TOP, fill=X)
                B3=Button(parent, text='Plot III', bg="green", command=plot3)
                B3.pack(side=TOP, fill=X)
                B4=Button(parent, text='Plot IV', bg="orange", command=plot4)
                B4.pack(side=TOP, fill=X)

            q_superpose3_gui(Root3)
            Root3.mainloop()

        def initialState_gui(parent):
            B1=Button(parent, text="|01>", bg="blue", command=q01)
            B1.pack(side=TOP, fill=X)
            B2=Button(parent, text="|10>", bg="red", command=q10)
            B2.pack(side=TOP, fill=X)
            B3=Button(parent, text="|11>", bg="green", command=q11)
            B3.pack(side=TOP, fill=X)
            B4=Button(parent, text="Alice:alive, Bob:superposed", bg="yellow", command=q_superpose1)
            B4.pack(side=TOP, fill=X)
            B5=Button(parent, text="Alice:superposed, Bob:alive", bg="magenta", command=q_superpose2)
            B5.pack(side=TOP, fill=X)
            B6=Button(parent, text="Alice:superposed, Bob:superposed", bg="orange", command=q_superpose3)
            B6.pack(side=TOP, fill=X)

        initialState_gui(Root2)
        Root2.mainloop()


    def twoPerson_gui(parent):
        B1=Button(parent, text="Initialise Qubits", bg="blue", command=initialState)
        B1.pack(side=TOP, fill=X, pady=4, padx=4)

    twoPerson_gui(Root1)
    Root1.mainloop()
    
master=Tk()

B1=Button(master,text="Two Person Duel", bg="red", command=twoPerson)
B1.pack(side=TOP, fill=X, pady=4, padx=4)
master.mainloop()       
