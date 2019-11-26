%%time

import numpy as np
from projectq import MainEngine  
from projectq.ops import H, All, Measure
from projectq.ops import QubitOperator, TimeEvolution

order4=[0.20724538589718786,
 0.4144907717943757,
 0.4144907717943757,
 0.4144907717943757,
 -0.12173615769156357,
 -0.6579630871775028,
 -0.12173615769156357,
 0.4144907717943757,
 0.4144907717943757,
 0.4144907717943757,
 0.20724538589718786]
order3=[7.0/24,2.0/3,3.0/4,-2.0/3,-1.0/24,1.0]
order2=[0.5,1,0.5]
order1=[1,1]


n_qubits = 13


def pauli_x(qubit, coeff):
    qubit = str(qubit) 
    return coeff*QubitOperator('X'+qubit)

def pauli_z(qubit, coeff):
    qubit = str(qubit) 
    return coeff*QubitOperator('Z'+qubit)

def product_pauli_z(q1, q2, coeff):
    q1,q2 = str(q1), str(q2)
    return coeff*QubitOperator(str('Z'+q1+' '+'Z'+q2))

def product_pauli_x(q1, q2, coeff):
    q1,q2 = str(q1), str(q2)
    return coeff*QubitOperator(str('X'+q1+' '+'X'+q2))
    
def evolution(steps=21,order=order4,T=3):
    eng = MainEngine()
    qureg = eng.allocate_qureg(n_qubits)
    dt = T/(steps-1)
    nc = len(order)

    for i in range(n_qubits):
        H | qureg[i] 
    
    for j in range(steps):
            s=(i+1)/steps
            for j in range(nc):
                if j%2==1:
                    TimeEvolution(dt*order[nc-j-1],Ham(s)[1]) | qureg
                else:
                    TimeEvolution(dt*order[nc-j-1],Ham(s)[0]) | qureg

    eng.flush()
    energy=eng.backend.get_expectation_value(H_c,qureg)
        
    All(Measure) | qureg
    return energy   


identity = pauli_x(0, 0)

H_0 = sum([pauli_x(i, -1/2) for i in range(n_qubits)], identity)

J = 2
alpha = 0.5
beta = 1.0
# J*T = 6????
# T = 3

H_c = sum([pauli_x(i, -alpha) for i in range(n_qubits)], identity)+\
    sum([pauli_z(i, -beta) for i in range(n_qubits)], identity)+\
    sum([product_pauli_x(i,i+1,-J) for i in range(n_qubits-1)], identity)+\
    sum([product_pauli_z(i,i+1,-J) for i in range(n_qubits-1)], identity)

def Ham(s):
    Ham = [0,0]
    Ham[0] =(1-s)*H_0+s*sum([pauli_x(i, -alpha) for i in range(n_qubits)], identity)+\
        s*sum([product_pauli_x(i,i+1,-J) for i in range(n_qubits-1)], identity)
    Ham[1] = s*(sum([pauli_z(i, -beta) for i in range(n_qubits)], identity)+\
        sum([product_pauli_z(i,i+1,-J) for i in range(n_qubits-1)], identity))
    return Ham

if __name__ == "__main__":
    gs_energy_slow=evolution(steps=210,order=order3,T=3)
    gs_energy_quick=evolution(steps=21,order=order4,T=3)
    
    print('Ground state energy by a slow evolution is {:.5f}, and a quick evolution is {:.5f}'.format(gs_energy_slow,gs_energy_quick))
    print('Error: {:.3%} '.format(abs(gs_energy_quick-gs_energy_slow)/gs_energy_slow))
