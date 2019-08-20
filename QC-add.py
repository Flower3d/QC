import numpy as np
import math as m

from qiskit import ClassicalRegister,QuantumRegister,QuantumCircuit,execute,Aer
from qiskit.tools.visualization import circuit_drawer


S_simulator = Aer.backends(name='statevector_simulator')[0]
M_simulator = Aer.backends(name='qasm_simulator')[0]


def QFT(qc,q,qubits):
    '''
    Assigns all the gate operations for a Ouantum Fourier Transformation
    '''
    R_phis = [0]
    for i in np.arange(2,int(qubits+1)):
        R_phis.append(2/(2**(i))*m.pi)
    for j in np.arange(int(qubits)):
        qc.h(q[int(j+qubits)])
        for k in np.arange(j+1,int(qubits)):
            qc.cu1(R_phis[k],q[int(k)+qubits],q[int(j)+qubits])


def QFT_dgr(qc,q,qubits):
    '''
    Assigns all the gate operations for an inverse Quantum Fourier Transfornation
    '''
    R_phis = [0]
    for i in np.arange(2,int(qubits+1)):
        R_phis.append(-2/(2**(i))*m.pi)
    for j in np.arange(int(qubits)):
        for k in np.arange(int(j)):
            qc.cu1(R_phis[int(qubits-(k+1))],q[int(qubits-(k+1))+qubits],q[int(qubits-(j+1))+qubits])
        qc.h(q[int(qubits-(j+1)+qubits)])


print('b>a, b&a are all decimal')
a = bin(int(input('a=')))
b = bin(int(input('b=')))
# b > a
n = len(b)-2

ab = QuantumRegister(2*n,name='ab')
c = ClassicalRegister(1,name='c')
abc = QuantumCircuit(ab,c,name='ac')

alst = list(a)[2:]
blst = list(b)[2:]
while(len(alst)<n):
    alst.insert(0,'0')


#Initialize a,b
bgates = [k for k in range(len(blst)) if int(blst[k]) != 0]
agates = [k for k in range(len(alst)) if int(alst[k]) != 0]
for k in range(len(bgates)):
    abc.x(ab[bgates[k]])
for k in range(len(agates)):
    abc.x(ab[agates[k]+n])


QFT(abc,ab,n)

R_phis = [0]
for i in np.arange(n):
    R_phis.append(1/(2**(i))*m.pi)
for k in np.arange(1,n+1):
    for j,l in zip(np.arange(k-1,n),np.arange(1,n-k+2)):
        abc.cu1(R_phis[l],ab[int(n+k-1)],ab[int(j)])
        
QFT_dgr(abc,ab,n)

circuit_drawer(abc)

job = execute(abc, S_simulator)
result = job.result()
result.get_statevector()

for i in range(4**n):
    if result.get_statevector()[i] != 0j:
        measure_result = str(bin(i))[2:]
        #print(i)
while len(measure_result)<2*n:
    measure_result = '0'+ measure_result
#Reverse
result_a,result_b = int(measure_result[n-1::-1],2),int(measure_result[2*n:n-1:-1],2)
#result_a,result_b --> QFT(a+b), b --> a+b, b

print('a+b='+ str(result_a))


