import itertools
import numpy as np


from qiskit import BasicAer as Aer
from qiskit import QuantumRegister, execute
from qiskit.quantum_info import Pauli
from qiskit import QuantumCircuit, ClassicalRegister, QuantumRegister
from qiskit.aqua import Operator, get_aer_backend
from qiskit.aqua.components.initial_states import Custom


n_qubits = 13


np.set_printoptions(precision=3, suppress=True)
backend = Aer.get_backend('statevector_simulator')
qr = QuantumRegister(n_qubits)
cr = ClassicalRegister(n_qubits)
circuit = QuantumCircuit(q, c)


def pauli_x(qubit, coeff):
    eye = np.eye((n_qubits)) 
    return Operator([[coeff, Pauli(np.zeros(n_qubits), eye[qubit])]])


def pauli_z(qubit, coeff):
    eye = np.eye((n_qubits))
    return Operator([[coeff, Pauli(eye[qubit], np.zeros(n_qubits))]])


def product_pauli_z(q1, q2, coeff):
    eye = np.eye((n_qubits))
    return Operator([[coeff, Pauli(eye[q1], np.zeros(n_qubits)) * Pauli(eye[q2], np.zeros(n_qubits))]])


def product_pauli_x(q1, q2, coeff):
    eye = np.eye((n_qubits))
    return Operator([[coeff, Pauli(np.zeros(n_qubits), eye[q1]) * Pauli(np.zeros(n_qubits), eye[q2])]])


#angle: the parameter t in $e^{iHt}$
def evolve(hamiltonian, angle, quantum_registers):
    return hamiltonian.evolve(None, angle, 'circuit', 1,
                              quantum_registers=quantum_registers,
                              expansion_mode='suzuki',
                              expansion_order=3)


def create_circuit(beta, gamma):
    circuit_evolv = sum([evolve(H_c, beta[i], qr) + evolve(H_c, gamma[i], qr)
                            for i in range(p)], evolve(identity, 0, qr))
    circuit = circuit_init + circuit_evolv
    return circuit


def evaluate_circuit(beta_gamma):
    n = len(beta_gamma)//2
    circuit = create_circuit(beta_gamma[:n], beta_gamma[n:])
    return np.real(H_c.eval("matrix", circuit, get_aer_backend('statevector_simulator'))[0])


identity = pauli_x(0, 0)


#transverse field Ising model Hamiltonian, B=1/2
H_0 = sum([pauli_x(i, -1/2) for i in range(n_qubits)], identity)
H_0.to_matrix()

J = 2
alpha = 0.5
beta = 1.0
# J*T = 6????
H_c = sum([pauli_x(i, -alpha) for i in range(n_qubits)], identity)+\
    sum([pauli_z(i, -beta) for i in range(n_qubits)], identity)+\
    sum([product_pauli_z(i,j,-J) for i,j in itertools.product(range(n_qubits), repeat=2)], identity)+\
    sum([product_pauli_x(i,j,-J) for i,j in itertools.product(range(n_qubits), repeat=2)], identity)
H_c.to_matrix()

#E
λ, v = np.linalg.eigh(H_0)

p = 2

#
init_state_vect = [1 for i in range(2**n_qubits)]
init_state = Custom(n_qubits, state_vector=init_state_vect)
circuit_init = init_state.construct_circuit('circuit', qr)




