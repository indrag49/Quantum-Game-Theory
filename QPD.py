from QSimulator import *

def qpd(U_Alice, U_Bob, w, x, y, z):
    sigma_x=PauliX(I2)
    U=(np.kron(I2, I2)+1j*np.kron(sigma_x, sigma_x))/np.sqrt(2)
    U_dag=np.conj(U.T)

    initial=U.dot(Q00)
    
    PsiS=np.kron(U_Alice, U_Bob).dot(initial)
    PsiF=U_dag.dot(PsiS)
##    plot_measure(measure(PsiF))
    
    cpsif=np.conj(PsiF.T)

    def pi(w, x, y, z): return (w*np.abs(cpsif.dot(Q00))**2+
    y*np.abs(cpsif.dot(Q01))**2+z*np.abs(cpsif.dot(Q10))**2+
    x*np.abs(cpsif.dot(Q11))**2, w*np.abs(cpsif.dot(Q00))**2+
    z*np.abs(cpsif.dot(Q01))**2+y*np.abs(cpsif.dot(Q10))**2+
    x*np.abs(cpsif.dot(Q11))**2) 

    return pi(w, x, y, z)



